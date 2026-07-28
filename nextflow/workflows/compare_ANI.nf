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

include { assertRank; writeNamesTsv; gatedGlobIn; toManifest; asmidRowFilter } from '../modules/common/utils.nf'
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
    // --asmid restricts the run the same way funannotate.nf/earlgrey_mask.nf do,
    // so representative-picking off this run's ANI data can't be influenced by
    // a strain that --asmid excluded (see .github issue #3).
    ANI_SAMPLES(params.samples, compareRank, '', asmidRowFilter())

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

        // Collect ASMIDs from all samples so CONCAT_ANI_TSVS can write the
        // companion asmid manifest (used for staleness checking on re-runs).
        def all_asmids = ANI_SAMPLES.out.samples
            .map { ignored, meta -> meta.asmid }
            .collect()
            .toSet()
        WRITE_ASMID_SET_COMPARE(all_asmids)

        // Build a manifest of every group's ANI TSV directly from the channel.
        // Two suffixes are possible:  skani/mash/sourmash emit
        //   <group>.full.ani.tsv   (SKANI_COMPARE writes this to batches/ subdir)
        //   <group>.ani.tsv        (fastani's MERGE_ANI writes this after batching)
        // Using the live channel avoids both the Channel.fromPath() deadlock inside
        // operator closures and a race against publishDir copying.
        def globbed = ani_tsv_ch
            .map { _group, tsv -> tsv }
            .collect()
            .map { tsv_list ->
                def allFiles = tsv_list.flatten().findAll { it.size() > 0 }.unique()
                log.info "Found ${allFiles.size()} ANI TSVs for concat"
                allFiles
            }
            .collectFile(name: 'tsv_manifest.txt', newLine: true)

        CONCAT_ANI_TSVS(globbed, WRITE_ASMID_SET_COMPARE.out.asmids)

        // Write predict_input to a single file using the process
        WRITE_ANI_PREDICT_INPUT(predict_input_tuples)

        PICK_REPRESENTATIVE_STRAIN(
            CONCAT_ANI_TSVS.out.out.ifEmpty(file('/dev/null')),
            WRITE_ANI_PREDICT_INPUT.out.tsv,
            file(params.samples),
            // Derive BUSCO dir from genome_stats_outdir rather than requiring an
            // explicit busco_genome_dir param.  This is consistent with how
            // --pipeline busco_genome names its output and prevents silent
            // misalignment when the two pipelines use different outdir names.
            file(params.genome_stats_outdir as String + "/BUSCO_genome")
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

// Write the run's asmid set to a one-per-line file so CONCAT_ANI_TSVS can
// read it as a val channel and write the companion asmid manifest alongside
// the merged TSV (used for staleness checking when re-running with new strains).
process WRITE_ASMID_SET_COMPARE {
    tag   "WRITE_ASMID_SET_COMPARE"
    label 'report'

    input:
        val asmid_set  // Set of ASMIDs in this run

    output:
        path("asmid_set.txt"), emit: asmids

    script:
    def asmid_list = asmid_set.collect { "\"${it}\"" }.join(' ')
    """
    printf '%s\\n' ${asmid_list} | tr ' ' '\\n' | sort -u > asmid_set.txt
    """

    stub:
    """
    touch asmid_set.txt
    """
}
