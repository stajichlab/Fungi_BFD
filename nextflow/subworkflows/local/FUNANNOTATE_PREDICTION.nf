//
// FUNANNOTATE_PREDICTION — ab-initio reuse resolution, predict.
//
// ANI-driven species-level reuse (--run_ani_reuse=true):
//   - Requires abinitio_reuse_assignments.csv to exist; the pipeline fails early
//     if it's missing so runs don't silently fall back to independent training.
//   - Representatives predict first (no -p), producing parameters.json on disk.
//   - Non-representatives predict with -p shared_params_json, reusing rep's training.
//   - All strains in a species share the same GBK freshness check via staleSharedParams.
//
// Pure independent training (--run_ani_reuse=false):
//   - abinitioReuseMap is empty; all strains train without shared params.
//   - No CSV is required.
//

include { FUNANNOTATE_PREDICT }  from '../../modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf'

include { gbkResult; staleRnaseq; sharedParamsJsonFor; staleSharedParams } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_PREDICTION {
    take:
    predict_input_ch       // tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
    abinitioReuseMap       // out -> [species, reuse_eligible, is_representative], loaded offline
    forceIndependentSet    // species that always train independently

    main:
    def predict_todo = predict_input_ch
        .map { out, asmid, sp, st, lt, bl, hl, tt, gfa ->
            def assignment   = abinitioReuseMap[out]
            def eligible     = assignment?.reuse_eligible && !forceIndependentSet.contains(sp)
            def is_rep       = assignment?.is_representative ?: false
            def sharedJson   = (eligible && !is_rep) ? sharedParamsJsonFor(sp) : null
            tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, sharedJson?.toString() ?: '')
        }
        .filter { out, _asmid, sp, _st, _lt, _bl, _hl, _tt, _gfa, shared_params_json ->
            gbkResult("${params.target}/${out}/predict_results", out as String) == null ||
                staleRnaseq(out as String, sp as String) ||
                staleSharedParams(out as String, shared_params_json ? file(shared_params_json as String) : null)
        }

    FUNANNOTATE_PREDICT(predict_todo)

    emit:
    metadata = FUNANNOTATE_PREDICT.out.metadata
}
