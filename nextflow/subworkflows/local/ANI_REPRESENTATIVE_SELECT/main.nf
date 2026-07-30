//
// ANI_REPRESENTATIVE_SELECT — compute ANI to drive species-level ab-initio reuse
// in the funannotate predict step. Fully wired into the pipeline.
//
// Flow:
//   1. Conditionally run ANI_COMPARE_METHOD (skani/mash/sourmash/fastani) to
//      produce per-species ANI TSVs if the merged TSV doesn't exist yet.
//   2. CONCAT_ANI_TSVS merges all group TSVs into one all-pairs TSV.
//   3. PICK_REPRESENTATIVE_STRAIN reads the merged TSV + predict_input,
//      picks one representative per species, writes abinitio_reuse_assignments.csv
//      (same format as species_reuse_clusters.py), and backfills the shared
//      ab-initio parameter store from any already-predicted representatives.
//
// emit:
//   assignments_csv  — path to abinitio_reuse_assignments.csv
//                      (FUNANNOTATE_PREDICTION reads this via loadAbinitioReuseMap)
//
// Required params:
//   ani_method, compare, outdir, samples, genome_dir, target,
//   gene_prediction_shared_abinitio, busco_genome_dir, ani_reuse_threshold,
//   run_ani_reuse
//

include { ANI_SAMPLES }                from '../../../subworkflows/local/ANI_SAMPLES.nf'
include { ANI_COMPARE_METHOD }         from '../../../subworkflows/local/ANI_COMPARE_METHOD.nf'
include { CONCAT_ANI_TSVS }            from '../../../modules/ani/report/CONCAT_ANI_TSVS/main.nf'
include { PICK_REPRESENTATIVE_STRAIN } from '../../../modules/ani/report/PICK_REPRESENTATIVE_STRAIN/main.nf'

workflow ANI_REPRESENTATIVE_SELECT {
    take:
    predict_input  // tuple(out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa)

    main:
    def method     = (params.ani_method as String).toLowerCase()
    // Representative selection is always species-level; default to SPECIES so
    // callers don't have to pass --compare just for this step.
    def compare    = (params.compare ?: 'SPECIES') as String
    def ani_out    = "${params.outdir}/${params.ani_method}/${compare}"
    def merged_tsv = file("${ani_out}/all_pairs_merged.tsv")
    def asmid_manifest_file = file("${ani_out}/all_pairs_merged.asmid_manifest.txt")

    def skip_ani   = !params.run_ani_reuse.toBoolean()

    // Staleness check:  require BOTH the merged TSV and its companion asmid
    // manifest to exist, and the manifest must contain every asmid in this run.
    // If a new strain was added to samples.csv the manifest will be stale and
    // ANI will be recomputed — fixing the params-only check bug.
    def run_asmid_set = predict_input.map { it[1] }.collect().map { it.toSet() }
    def has_merged    = merged_tsv.exists() && merged_tsv.size() > 0 &&
                        asmid_manifest_file.exists() && asmid_manifest_file.size() > 0 &&
                        isAsmidManifestCurrent(asmid_manifest_file, run_asmid_set.val)

    if (skip_ani) {
        log.info "ani_reuse: --run_ani_reuse=false — all strains will train independently"
    } else {
        log.info "ani_reuse: ${has_merged ? "merged TSV + asmid manifest current — using existing ANI data" : "computing ANI for representative selection"}"
    }

    // ── 1. ANI compute ─────────────────────────────────────────────────────────
    if (!skip_ani && !has_merged) {
        // Filtering to this run's asmids happens below via run_asmid_set, so no
        // extra --asmid filter is needed at the ANI_SAMPLES level here.
        def samples_ch = ANI_SAMPLES(params.samples, 'SPECIES', '', { _row -> true })
            .samples
            .groupTuple()
            .map { sp, metas ->
                def run_genomes = metas.findAll { m -> run_asmid_set.val.contains(m.asmid) }
                run_genomes ? tuple(sp, run_genomes.collect { m -> m.genome }) : null
            }
            .filter { it != null }

        ANI_COMPARE_METHOD(samples_ch, method)

        // Write the asmid set to a file so CONCAT_ANI_TSVS can read it as a
        // channel input and write the companion manifest next to the merged TSV.
        WRITE_ASMID_SET(run_asmid_set)

        // Build a manifest of every group's ANI TSV directly from the channel.
        // Using the live channel avoids both the Channel.fromPath() deadlock inside
        // operator closures and a race against publishDir copying.
        //
        // IMPORTANT: feed collectFile() a channel of *string* path emissions, one
        // per group — do NOT .collect() them into a single List first. See the
        // matching comment in workflows/compare_ANI.nf (same bug, same fix).
        def ani_tsv_paths = ANI_COMPARE_METHOD.out.ani_tsv
            .map { _group, tsv -> tsv }
            .filter { it.size() > 0 }

        ani_tsv_paths.collect().subscribe { log.info "Found ${it.size()} ANI TSVs for concat" }

        def tsv_manifest = ani_tsv_paths
            .map { it.toString() }
            .unique()
            .collectFile(name: 'tsv_manifest.txt', newLine: true, sort: true)

        CONCAT_ANI_TSVS(tsv_manifest, WRITE_ASMID_SET.out.asmids)
    }

    def ani_tsv = (skip_ani || has_merged)
        ? merged_tsv
        : CONCAT_ANI_TSVS.out.out

    // ── 2. Representative selection + CSV + backfill ──────────────────────────
    // Write predict_input to a process-published TSV so PICK_REPRESENTATIVE_STRAIN
    // can read it as a file input (avoids collect() + GString join pitfall).
    WRITE_PREDICT_INPUT(predict_input)

    PICK_REPRESENTATIVE_STRAIN(
        ani_tsv.ifEmpty(file('/dev/null')),
        WRITE_PREDICT_INPUT.out.tsv,
        file(params.samples),
        // Derive BUSCO dir from genome_stats_outdir rather than requiring an
        // explicit busco_genome_dir param.  Consistent with --pipeline
        // busco_genome's output naming; no silent misalignment possible.
        file(params.genome_stats_outdir as String + "/BUSCO_genome")
    )

    emit:
    assignments_csv = PICK_REPRESENTATIVE_STRAIN.out.outCSV
}

// Return true when every asmid in `current_set` appears in `manifest_file`.
// Missing entries mean a new strain was added → cached ANI is stale.
def isAsmidManifestCurrent(File manifestFile, Set currentSet) {
    def stored = manifestFile.readLines()*.trim().findAll { it } as Set
    currentSet.every { stored.contains(it) }
}

// Write the run's asmid set to a one-per-line file so CONCAT_ANI_TSVS can
// read it as a val channel and write the companion asmid manifest.
process WRITE_ASMID_SET {
    tag   "WRITE_ASMID_SET"
    label 'report'

    input:
        val asmid_set  // Set of ASMIDs for this run

    output:
        path("run_asmid_set.txt"), emit: asmids

    script:
    def asmid_list = asmid_set.collect { it }.join('\n')
    """
    printf '%s\\n' '${asmid_list}' | grep -v '^\${\$}' | sort -u > run_asmid_set.txt
    """

    stub:
    """
    touch run_asmid_set.txt
    """
}

// Write predict_input channel to a TSV file for PICK_REPRESENTATIVE_STRAIN.
// publishDir makes Nextflow stage the file so the downstream process can read it.
process WRITE_PREDICT_INPUT {
    tag   "WRITE_PREDICT_INPUT"
    label 'report'

    publishDir "${workflow.workDir}", mode: 'move', overwrite: true

    input:
        tuple val(out), val(asmid), val(sp), val(st), val(lt), val(bl), val(hl), val(tt), val(gfa)

    output:
        path("predict_input_for_ani.tsv"), emit: tsv

    script:
    def line = "${out}\t${asmid}\t${sp}\t${st}\t${lt}\t${bl}\t${hl}\t${tt}\t${gfa}"
    """
    printf 'out\\tasmid\\tspecies\\tstrain\\tlocustag\\tbusco\\thlen\\ttable\\tgenome_fa\\n' > predict_input_for_ani.tsv
    printf '${line}\\n' >> predict_input_for_ani.tsv
    """
}
