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
// Symlink-resolved absolute path; falls back to the absolute path when the target
// does not exist yet (the shared store is created by beforeScript).
def backfillRealPath(p) {
    def f = file(p.toString())
    return f.exists() ? f.toRealPath().toString() : f.toAbsolutePath().toString()
}

process BACKFILL_ABINITIO_PARAMS {
    tag   "batch_$batch_id"
    label 'report'

    // The 'report' container only binds projectDir, so params.target (the
    // representatives' predict_results), the shared store and the Augustus config
    // were invisible inside it: backfill_abinitio_params.py reported every
    // representative "has no usable prediction on disk", exited 1, and the
    // following `touch` hid the failure (wave 0, 2026-09-27: 153 ms, empty store,
    // 21 siblings silently trained independently). Bind them explicitly, the same
    // set FUNANNOTATE_PREDICT binds; workDir covers the GENEMARK_RUN .mod paths.
    // beforeScript runs on the host, so the store exists before it is bound.
    // set -euo pipefail below catches shell errors; a per-species backfill failure is
    // reported as a warning instead, because the batch's other species did succeed.
    beforeScript "mkdir -p '${params.gene_prediction_shared_abinitio}'"
    // Real (symlink-resolved) paths, used for both the binds and the script
    // arguments: a symlinked path (e.g. <launchDir>/lib -> ../lib) cannot be a
    // bind destination on this GPFS mount ("destination ... doesn't exist in
    // container", wave 0 2026-09-27), and would not resolve inside the container.
    containerOptions {
        def binds = [params.target, params.gene_prediction_shared_abinitio, workflow.workDir.toString()]
        if (params.augustus_config) binds << params.augustus_config
        '--bind ' + binds.collect { b -> backfillRealPath(b) }.unique().collect { b -> "${b}:${b}" }.join(',')
    }

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
    // projectDir is usually reached through a symlink (<launchDir>/nextflow ->
    // Fungi_BFD/nextflow); apptainer.runOptions binds only its real path, so the
    // symlinked script path did not exist in the container ("can't open file",
    // the 153 ms failure in wave 0).
    def shared_root  = backfillRealPath(params.gene_prediction_shared_abinitio)
    def target       = backfillRealPath(params.target)
    def threshold    = params.ani_reuse_threshold ?: 99.0
    def aug_cfg      = params.augustus_config ? backfillRealPath(params.augustus_config) : ''
    // Manifest content is fully resolved in Groovy before the shell ever sees
    // it, so there's no bare `$` in this string for Groovy to misinterpolate.
    def manifest = items.collect { sp, out, mod -> "${sp}\t${out}\t${mod ?: ''}" }.join('\n')
    """
    set -euo pipefail
    cat > manifest.tsv <<'BACKFILL_MANIFEST_EOF'
${manifest}
BACKFILL_MANIFEST_EOF
    python "${backfillRealPath(projectDir)}/bin/backfill_abinitio_params.py" \
        --manifest manifest.tsv \
        --target "${target}" \
        --shared-root "${shared_root}" \
        --ani-threshold "${threshold}" \
        ${aug_cfg ? "--augustus-config ${aug_cfg}" : ''} \\
        || echo "[WARN] backfill_abinitio_params.py reported failures (see [ERROR] lines above); siblings of those species stay blocked -- FUNANNOTATE_PREDICTION checks each store on disk" >&2
    touch "batch_${batch_id}.backfill.done"
    """

    stub:
    """
    touch "batch_${batch_id}.backfill.done"
    """
}
