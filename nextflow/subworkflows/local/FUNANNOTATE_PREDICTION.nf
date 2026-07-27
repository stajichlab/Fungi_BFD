//
// FUNANNOTATE_PREDICTION — ab-initio reuse resolution, predict.
//
// Resolves per-row whether a genome should reuse shared ab-initio
// parameters (species-level reuse, todo/species_level_abinitio_reuse.md)
// or train independently, then runs FUNANNOTATE_PREDICT on every genome
// whose prediction is missing or stale.
//
// Representative-first ordering:
//   1. Predict on representatives (no shared params) → produces parameters.json
//   2. Predict on non-representatives with -p parameters.json (reuses rep's training)
//   3. All strains in a species share the same GBK freshness check via staleSharedParams.
//
// When ANI data is absent (ani.db not computed), all strains train independently.
// The abinitioReuseMap is loaded from abinitio_reuse_csv; if the CSV was written
// by ANI_REPRESENTATIVE_SELECT it includes the is_representative column.
//

include { FUNANNOTATE_PREDICT }  from '../../modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf'

include { gbkResult; staleRnaseq; sharedParamsJsonFor; staleSharedParams } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_PREDICTION {
    take:
    predict_input_ch       // tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
    abinitioReuseMap       // out -> [species, reuse_eligible, is_representative], loaded offline
    forceIndependentSet    // species that always train independently

    main:
    // Attach shared_params_json (empty string = train independently) and
    // is_representative per genome.
    def tagged = predict_input_ch
        .map { out, asmid, sp, st, lt, bl, hl, tt, gfa ->
            def assignment   = abinitioReuseMap[out]
            def eligible     = assignment?.reuse_eligible && !forceIndependentSet.contains(sp)
            def is_rep       = assignment?.is_representative ?: false
            def sharedJson   = (eligible && !is_rep) ? sharedParamsJsonFor(sp) : null
            tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, is_rep, sharedJson?.toString() ?: '')
        }

    // Which strains need a fresh predict run (GBK missing or stale)?
    def predict_todo = tagged
        .filter { out, _asmid, sp, _st, _lt, _bl, _hl, _tt, _gfa, _is_rep, shared_params_json ->
            gbkResult("${params.target}/${out}/predict_results", out as String) == null ||
                staleRnaseq(out as String, sp as String) ||
                staleSharedParams(out as String, shared_params_json ? file(shared_params_json as String) : null)
        }

    // ── Representative-first: run reps first so their params are on disk ──────
    // before non-representatives look for them via sharedParamsJsonFor().
    def repr_todo    = predict_todo.filter { it[9] == true  }  // is_representative
    def nonrep_todo  = predict_todo.filter { it[9] == false }  // is_representative

    FUNANNOTATE_PREDICT(repr_todo)    // no shared_params_json (it[10] = '')
    FUNANNOTATE_PREDICT(nonrep_todo)  // shared_params_json set for reuse_eligible

    emit:
    metadata = FUNANNOTATE_PREDICT.out.metadata
}
