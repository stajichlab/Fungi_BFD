//
// BACKFILL_ABINITIO_PARAMS — copy representative strains' trained AUGUSTUS/
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
// Batched, not one task per representative: the per-species work here is
// seconds of I/O (copy + content-hash compare), so with hundreds of
// representatives, one-SLURM-job-per-species was dominated by submission/
// queue overhead rather than actual work. Callers group candidates with
// .collate(100) (see workflows/backfill_abinitio.nf and
// FUNANNOTATE_PREDICTION.nf) before calling this process, so one task loops
// over up to ~100 species inside a single job via backfill_abinitio_params.py
// --manifest. Each line is independent and idempotent (content-hash
// short-circuit in backfill_species_store()), so re-running a whole batch on
// retry/-resume just re-confirms already-backfilled entries as up to date.
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
    tag   "batch_$batch_id"
    label 'report'

    input:
        // items: List of [species, rep_out, genemark_mod] -- 3rd field is ''
        // when GENEMARK_RUN is off (run_genemark=false) or the caller is the
        // standalone backfill_abinitio sweep (representatives predicted before
        // GENEMARK_RUN existed; backfill_abinitio_params.py falls back to
        // deriving genemark.mod from predict_misc/ when this field is empty).
        tuple val(batch_id), val(items)

    output:
        tuple val(items), path("*.backfill.done"), emit: done

    script:
    def shared_root  = params.gene_prediction_shared_abinitio
    def target       = params.target
    def threshold    = params.ani_reuse_threshold ?: 99.0
    def aug_cfg      = params.augustus_config ?: ''
    // Manifest content is fully resolved in Groovy before the shell ever sees
    // it, so there's no bare `$` in this string for Groovy to misinterpolate.
    def manifest = items.collect { sp, out, mod -> "${sp}\t${out}\t${mod ?: ''}" }.join('\n')
    """
    cat > manifest.tsv <<'BACKFILL_MANIFEST_EOF'
${manifest}
BACKFILL_MANIFEST_EOF
    python "${projectDir}/bin/backfill_abinitio_params.py" \
        --manifest manifest.tsv \
        --target "${target}" \
        --shared-root "${shared_root}" \
        --ani-threshold "${threshold}" \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''}
    touch "batch_${batch_id}.backfill.done"
    """

    stub:
    """
    touch "batch_${batch_id}.backfill.done"
    """
}
