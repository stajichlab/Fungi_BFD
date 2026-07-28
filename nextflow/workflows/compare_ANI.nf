//
// compare_ANI — symmetric all-vs-all ANI within taxonomic groups.
//
// Groups samples.csv by --compare rank, then sketches and compares every genome
// in each group with --ani_method (skani | mash | sourmash | fastani).
//
// When --run_ani_reuse=true, a final representative-selection step runs after
// ANI compute: CONCAT_ANI_TSVS merges per-group TSVs into all_pairs_merged.tsv,
// then PICK_REPRESENTATIVE_STRAIN writes genome_annotation/_reuse_assignments/
// abinitio_reuse_assignments.csv for consumption by --pipeline funannotate.
//
// The funannotate pipeline has two independent pre-requisites that use only
// raw genome sequence (no annotation required):
//   1. compare_ani --run_ani_reuse true  → abinitio_reuse_assignments.csv
//   2. busco_genome pipeline             → results/genome_stats/BUSCO_genome/*.BUSCO_summary.txt
// Both must complete before running --pipeline funannotate --run_ani_reuse true.
//
// See workflows/query_ANI.nf for the asymmetric orphan-placement counterpart.
//
// Usage (from project root):
//   # Pre-requisites
//   nextflow run nextflow/main.nf --pipeline compare_ani -profile ani \
//       --compare SPECIES --run_ani_reuse true -resume
//   nextflow run nextflow/main.nf --pipeline busco_genome -profile BFD -resume
//
//   # Prediction (requires both pre-requisites above)
//   nextflow run nextflow/main.nf --pipeline funannotate -profile funannotate \
//       --run_ani_reuse true -resume
//
// Parameter defaults live in conf/profile_ANI.config.
//

include { assertRank; writeNamesTsv; gatedGlobIn; toManifest } from '../modules/common/utils.nf'
include { ANI_SAMPLES }               from '../subworkflows/local/ANI_SAMPLES.nf'
include { ANI_COMPARE_METHOD }        from '../subworkflows/local/ANI_COMPARE_METHOD.nf'
include { REPORT_ANI }                from '../modules/ani/report/REPORT_ANI/main.nf'
include { COMBINE_ANI_TABLE }         from '../modules/ani/report/COMBINE_ANI_TABLE/main.nf'
include { CONCAT_ANI_TSVS }           from '../modules/ani/report/CONCAT_ANI_TSVS/main.nf'
include { PICK_REPRESENTATIVE_STRAIN } from '../modules/ani/report/PICK_REPRESENTATIVE_STRAIN/main.nf'

workflow COMPARE_ANI {

    // ── Validate params ───────────────────────────────────────────────────────
    def compareRank = assertRank(params.compare as String, 'compare')

    def method = (params.ani_method as String).toLowerCase()
    if (!(method in ['skani','mash','sourmash','fastani'])) {
        error "--ani_method must be one of: skani, mash, sourmash, fastani"
    }

    log.info "ANI method: ${method}"
    log.info "Genome dir: ${params.genome_dir}"
    if (method == 'fastani' && (params.fastani_prefilter as boolean)) {
        log.info "fastANI mash-prefilter cascade ENABLED (prefilter_ani=${params.prefilter_ani}%)"
    }

    // ── Samples → groups ──────────────────────────────────────────────────────
    // query_rank is '' here: compare_ANI treats every genome as a reference.
    ANI_SAMPLES(params.samples, compareRank, '')

    // n_test limits *groups* (applied after groupTuple so --n_test 3 = 3 groups).
    def grouped_ch = ANI_SAMPLES.out.samples
        .groupTuple()
        .filter { _gname, metas -> metas.size() >= params.min_group_size as int }
        .take(params.n_test > 0 ? params.n_test as int : -1)

    // ── Per-group genome list + names TSV ─────────────────────────────────────
    def prepared_ch = grouped_ch.map { group_name, metas ->
        tuple(group_name, metas.collect { m -> m.genome }, writeNamesTsv(group_name, metas))
    }

    def prepared_split = prepared_ch.multiMap { group_name, genomes, nameFile ->
        ani:   tuple(group_name, genomes)
        names: tuple(group_name, nameFile)
    }
    def genome_ch = prepared_split.ani
    def names_map = prepared_split.names

    // ── Sketch + compare ──────────────────────────────────────────────────────
    ANI_COMPARE_METHOD(genome_ch, method)
    def ani_tsv_ch = ANI_COMPARE_METHOD.out.ani_tsv

    // ── Reports ───────────────────────────────────────────────────────────────
    REPORT_ANI(ani_tsv_ch.join(names_map, by: 0))

    // ── Combine every group's pairs + labels into one queryable table ─────────
    def ani_sync   = ani_tsv_ch.collect().map { true }.ifEmpty(true)
    def names_sync = REPORT_ANI.out.names.collect().map { true }.ifEmpty(true)

    def ani_all = gatedGlobIn(ani_sync, params.outdir, "${params.ani_method}/${params.compare}/**/*.ani.tsv")
    def names_all = gatedGlobIn(names_sync, params.outdir, "${params.ani_method}/${params.compare}/**/*_genome_names.tsv")

    COMBINE_ANI_TABLE(
        toManifest(ani_all,   'ani_tsv.manifest.txt'),
        toManifest(names_all, 'names_tsv.manifest.txt')
    )

    // ── Representative strain selection (ANI-driven ab-initio reuse) ──────────
    // Runs after ANI compute when --run_ani_reuse=true.
    // Produces abinitio_reuse_assignments.csv consumed by --pipeline funannotate.
    // BUSCO_genome results (results/genome_stats/BUSCO_genome/) must exist from
    // a prior --pipeline busco_genome run; pick_representative_strain.py reads
    // them from disk and does not require a separate Nextflow channel.
    if (params.run_ani_reuse.toBoolean()) {
        log.info "ani_reuse: computing representative strain assignments"

        // Build predict_input channel from ANI_SAMPLES metadata (same tuple shape
        // as FUNANNOTATE_RNASEQ.out.predict_input).
        def predict_input_tuples = ANI_SAMPLES.out.samples
            .map { group_name, meta ->
                def out = meta.id.replaceAll(/\.(fa\.gz|masked\.fasta\.gz|fasta\.gz)$/, '')
                tuple(
                    out, meta.asmid, meta.species, meta.strain, meta.locustag,
                    'fungi_odb12',   // busco lineage (not used at ANI stage)
                    24,              // header_length
                    '1',             // transl_table
                    meta.genome.toString()
                )
            }
            .collect()

        // Wait for all groups to finish, then glob their TSVs into a manifest.
        def merge_done = ani_tsv_ch.collect().map { true }.ifEmpty(false)
        def tsv_manifest = merge_done
            .map { dummy -> 
                def dirPath = "${params.outdir}/${params.ani_method}/${params.compare}"
                def dir = new File(dirPath)
                def allFiles = []
                // Search root directory
                allFiles.addAll(dir.listFiles()?.findAll { f -> 
                    f.isFile() && f.name.endsWith('.ani.tsv')
                } ?: [])
                // Search subdirectories (1 level deep)
                def subdirs = dir.listFiles()?.findAll { it.isDirectory() } ?: []
                subdirs.each { subdir ->
                    allFiles.addAll(subdir.listFiles()?.findAll { f -> 
                        f.isFile() && f.name.endsWith('.ani.tsv')
                    } ?: [])
                    // Also search 2 levels deep (e.g., batches/)
                    subdir.listFiles()?.findAll { it.isDirectory() }?.each { subsubdir ->
                        allFiles.addAll(subsubdir.listFiles()?.findAll { f -> 
                            f.isFile() && f.name.endsWith('.ani.tsv')
                        } ?: [])
                    }
                }
                def paths = allFiles.collect { it.absolutePath }.join('\n') + '\n'
                log.info "Found ${allFiles.size()} ANI files in ${dirPath}: ${allFiles.collect { it.name }}"
                paths
            }
            .collectFile(name: 'tsv_manifest.txt')

        CONCAT_ANI_TSVS(tsv_manifest)

        // Write predict_input to a single file using the process
        WRITE_ANI_PREDICT_INPUT(predict_input_tuples)

        PICK_REPRESENTATIVE_STRAIN(
            CONCAT_ANI_TSVS.out.out.ifEmpty(file('/dev/null')),
            WRITE_ANI_PREDICT_INPUT.out.tsv,
            file(params.samples),
            file(params.busco_genome_dir as String)
        )
    }
}

// Write all predict_input tuples to a single TSV for PICK_REPRESENTATIVE_STRAIN.
process WRITE_ANI_PREDICT_INPUT {
    tag   "WRITE_ANI_PREDICT_INPUT"
    label 'report'

    publishDir "${params.target}/_reuse_assignments", mode: 'copy'

    input:
        val rows  // list of [out, asmid, species, strain, locustag, busco, hlen, table, genome_fa]

    output:
        path("predict_input_for_ani.tsv"), emit: tsv

    script:
    def json_rows = new groovy.json.JsonBuilder(rows).toString()
    """
    python3 - << 'PYEOF'
import json

header = ['out', 'asmid', 'species', 'strain', 'locustag', 'busco', 'hlen', 'table', 'genome_fa']
with open('predict_input_for_ani.tsv', 'w') as f:
    f.write('\\t'.join(header) + '\\n')
    data = json.loads('${json_rows}')
    # Handle both nested [[...], [...]] and flat [...] structures
    if data and isinstance(data[0], list):
        rows = data
    else:
        # Flat list - chunk into groups of 9
        rows = [data[i:i+9] for i in range(0, len(data), 9)]
    for row in rows:
        f.write('\\t'.join(str(x) for x in row) + '\\n')
PYEOF
    """
}


