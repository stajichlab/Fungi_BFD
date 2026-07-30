//
// BACKFILL_ABINITIO_PARAMS — copy a representative strain's trained AUGUSTUS/
// GeneMark/SNAP ab-initio parameters into the shared per-species parameter
// store (params.gene_prediction_shared_abinitio), so sibling strains can
// reuse them via `funannotate predict -p parameters.json`.
//
// One process, two callers (nextflow/bin/backfill_abinitio_params.py has the
// actual logic — staging dir + content-hash idempotency + atomic swap):
//   - FUNANNOTATE_PREDICTION.nf: wired with a real channel dependency right
//     after a representative's own FUNANNOTATE_PREDICT, in --predict_scope
//     all mode.
//   - workflows/backfill_abinitio.nf: the standalone sweep pipeline, for
//     representatives predicted before this feature existed (or by a
//     representative_only run in a separate invocation).
//
// "Option B" persistence model (same as FUNANNOTATE_PREDICT): the real output
// is written directly into the persistent gene_prediction_shared_abinitio
// store, not into the Nextflow work dir. There is no publishDir copy. The
// emitted marker file only keeps the DAG edge alive for downstream channel
// joins; callers re-derive the actual shared-params path from disk via
// sharedParamsJsonFor() (funannotate/utils.nf) once this task completes,
// since by then the write is guaranteed complete (it happens synchronously,
// before the marker is touched, within the same script block).
//
process BACKFILL_ABINITIO_PARAMS {
    tag   "$species"
    label 'report'

    input:
        tuple val(species), val(rep_out)

    output:
        tuple val(species), val(rep_out), path("*.backfill.done"), emit: done

    script:
    def species_tag = (species as String).replaceAll(/\s+/, '_')
    def shared_root  = params.gene_prediction_shared_abinitio
    def target       = params.target
    def threshold    = params.ani_reuse_threshold ?: 99.0
    def aug_cfg      = params.augustus_config ?: ''
    """
    python "${projectDir}/bin/backfill_abinitio_params.py" \
        --species "${species}" \
        --representative-out "${rep_out}" \
        --target "${target}" \
        --shared-root "${shared_root}" \
        --ani-threshold "${threshold}" \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''}
    touch "${species_tag}.backfill.done"
    """

    stub:
    def species_tag = (species as String).replaceAll(/\s+/, '_')
    """
    touch "${species_tag}.backfill.done"
    """
}
