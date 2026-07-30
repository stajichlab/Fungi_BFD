//
// backfill_abinitio — standalone sweep: backfill the shared ab-initio
// parameter store for every representative strain that already has a
// prediction on disk, but whose store is missing or stale.
//
// Complements the wired backfill inside FUNANNOTATE_PREDICTION.nf's
// --predict_scope all mode (which only backfills representatives predicted
// in that same run). This pipeline catches every other case: representatives
// predicted before this feature existed, or via a `representative_only` run
// in a separate invocation. Reuses the exact same BACKFILL_ABINITIO_PARAMS
// process and backfill_abinitio_params.py script as that wired path -- one
// shared implementation, two callers.
//
// Discovery is CSV-driven only: every candidate comes from
// abinitio_reuse_assignments.csv's is_representative rows (computed by
// compare_ani/PICK_REPRESENTATIVE_STRAIN). This pipeline never re-derives
// who is representative -- it only reacts to what compare_ani already
// decided, so there is exactly one source of truth for that question.
//
// Usage:
//   nextflow run nextflow/main.nf --pipeline backfill_abinitio \
//       -c nextflow/nextflow.config -profile funannotate -resume
//

include { loadAbinitioReuseMap; gbkResult } from '../modules/funannotate/utils.nf'
include { BACKFILL_ABINITIO_PARAMS }        from '../modules/funannotate/predict/BACKFILL_ABINITIO_PARAMS/main.nf'

workflow BACKFILL_ABINITIO {
    // See the matching normalization + comment in workflows/funannotate.nf --
    // params.target gets passed on to BACKFILL_ABINITIO_PARAMS's task script,
    // where a relative path would resolve against that task's own work
    // directory instead of the persistent project tree.
    params.target = file(params.target as String).toAbsolutePath().toString()

    def abinitioReuseMap = loadAbinitioReuseMap()

    if (abinitioReuseMap.isEmpty()) {
        log.warn "backfill_abinitio: no ab-initio reuse assignments found (or " +
            "share_abinitio_params=false) -- nothing to backfill. Run compare_ani " +
            "(--pipeline compare_ani --compare SPECIES --run_ani_reuse true) first."
        return
    }

    // Plain Groovy collection filtering, not channel operators -- this whole
    // computation is synchronous (loadAbinitioReuseMap() already read the CSV
    // eagerly), so there's no combine()/collectFile() footgun to worry about
    // here; only the final Channel.fromList() hand-off is a real channel.
    def candidates = abinitioReuseMap
        .findAll { _out, v -> v.is_representative }
        .collect { out, v -> tuple(v.species, out) }
        .findAll { species, out ->
            def has_gbk = gbkResult("${params.target}/${out}/predict_results", out as String) != null
            if (!has_gbk) {
                log.info "backfill_abinitio: skipping ${out} (${species}) -- not yet predicted"
            }
            has_gbk
        }

    if (candidates.isEmpty()) {
        log.info "backfill_abinitio: no predicted representatives found to backfill"
        return
    }

    log.info "backfill_abinitio: backfilling ${candidates.size()} representative(s)"
    BACKFILL_ABINITIO_PARAMS(Channel.fromList(candidates))
}
