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
    def compare    = params.compare as String
    def ani_out    = "${params.outdir}/${params.ani_method}/${params.compare}"
    def merged_tsv = file("${ani_out}/all_pairs_merged.tsv")

    def skip_ani   = !params.run_ani_reuse.toBoolean()
    def has_merged = merged_tsv.exists() && merged_tsv.size() > 0

    if (skip_ani) {
        log.info "ani_reuse: --run_ani_reuse=false — all strains will train independently"
    } else {
        log.info "ani_reuse: ${has_merged ? "merged TSV found — using existing ANI data" : "computing ANI for representative selection"}"
    }

    // ── 1. ANI compute ─────────────────────────────────────────────────────────
    if (!skip_ani && !has_merged) {
        def run_asmid_set = predict_input.map { it[1] }.collect().map { it.toSet() }

        def samples_ch = ANI_SAMPLES(params.samples, 'SPECIES', '')
            .samples
            .groupTuple()
            .map { sp, metas ->
                def run_genomes = metas.findAll { m -> run_asmid_set.val.contains(m.asmid) }
                run_genomes ? tuple(sp, run_genomes.collect { m -> m.genome }) : null
            }
            .filter { it != null }

        ANI_COMPARE_METHOD(samples_ch, method)

        // Wait for all groups to finish, then glob their MERGE_ANI storeDir outputs.
        def merge_done = ANI_COMPARE_METHOD.out.ani_tsv
            .collect()
            .map { true }
            .ifEmpty(false)

        def tsv_glob = channel.fromPath("${ani_out}/*/*.ani.tsv")
            .filter { it.name.endsWith('.ani.tsv') && it.size() > 0 }
            .collectFile(name: 'tsv_manifest.txt', newLine: true)

        def merged_input = merge_done
            .combine(tsv_glob)
            .map { _done, manifest -> manifest }

        CONCAT_ANI_TSVS(merged_input)
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
        file(params.busco_genome_dir as String)
    )

    emit:
    assignments_csv = PICK_REPRESENTATIVE_STRAIN.out.outCSV
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
