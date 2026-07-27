//
// FUNANNOTATE_PREDICTION — ab-initio reuse resolution, predict.
//
// Resolves per-row whether a genome should reuse shared ab-initio
// parameters (species-level reuse, todo/species_level_abinitio_reuse.md)
// or train independently, then runs FUNANNOTATE_PREDICT on every genome
// whose prediction is missing or stale.
//
// Emits the metadata channel (tuple of out, asmid, species, strain,
// locustag, busco, header_length, transl_table) for genomes predicted
// in THIS run. Prior-run genomes with current predictions are handled
// by the caller's postpredict scan.
//

include { FUNANNOTATE_PREDICT  } from '../../modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf'

include { gbkResult; staleRnaseq; sharedParamsJsonFor; staleSharedParams } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_PREDICTION {
    take:
    predict_input_ch       // tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
    abinitioReuseMap       // out -> [species, reuse_eligible], loaded offline
    forceIndependentSet    // species that always train independently

    main:
    // Attach shared_params_json (empty string = train independently, today's
    // behavior). Resolved per-row from the offline abinitioReuseMap +
    // sharedParamsJsonFor(), not trusted blindly: reuse_eligible=true in the CSV
    // is a prediction made at species_reuse_clusters.py run time, not a guarantee
    // -- the representative's own PREDICT could have since failed/been skipped,
    // so this re-checks the actual store file exists on disk right now (plan
    // S4.3 "runtime re-verification" and "representative failure/skip fallback").
    def predict_ch = predict_input_ch
        .map { out, asmid, sp, st, lt, bl, hl, tt, genome_fa ->
            def assignment = abinitioReuseMap[out]
            def eligible = assignment?.reuse_eligible && !forceIndependentSet.contains(sp)
            def sharedJson = eligible ? sharedParamsJsonFor(sp) : null
            tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa, sharedJson?.toString() ?: '')
        }
        .filter { out, _asmid, sp, _st, _lt, _bl, _hl, _tt, _gfa, shared_params_json ->
            gbkResult("${params.target}/${out}/predict_results", out as String) == null ||
                staleRnaseq(out as String, sp as String) ||
                staleSharedParams(out as String, shared_params_json ? file(shared_params_json as String) : null)
        }

    FUNANNOTATE_PREDICT(predict_ch)

    emit:
    metadata = FUNANNOTATE_PREDICT.out.metadata
}
