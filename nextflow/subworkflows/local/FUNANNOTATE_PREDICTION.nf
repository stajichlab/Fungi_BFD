//
// FUNANNOTATE_PREDICTION — ab-initio reuse resolution, predict.
//
// ANI-driven species-level reuse (--run_ani_reuse=true):
//   - Requires abinitio_reuse_assignments.csv to exist; the pipeline fails early
//     if it's missing so runs don't silently fall back to independent training.
//
//   - --predict_scope representative_only: predict ONLY representative strains,
//     auto-derived from abinitioReuseMap's is_representative flag (no manual
//     --asmid hunting needed -- composes normally with --taxon/--asmid/--n_test).
//     Wires each representative's own FUNANNOTATE_PREDICT output straight into
//     BACKFILL_ABINITIO_PARAMS, so one invocation trains the rep AND populates
//     the shared store -- no separate compare_ani rerun needed for that leg.
//
//   - --predict_scope all (default): everyone. Representatives predict first (no
//     -p) and get backfilled the same way, in the same run. An eligible non-
//     representative sibling only proceeds if its species' shared params are
//     available -- either pre-existing on disk (from a prior representative_only
//     run, or the standalone --pipeline backfill_abinitio sweep), or freshly
//     backfilled by this same run's own representative. A species with neither
//     has its eligible siblings excluded (not trained independently) and
//     reported to predict_blocked_awaiting_representative.tsv -- per-species,
//     the rest of the run proceeds normally -- unless --allow_independent_fallback
//     is set, in which case they train independently instead (still logged).
//
// Pure independent training (--run_ani_reuse=false, or a species in
// forceIndependentSet): abinitioReuseMap has no eligible/representative
// assignment for those strains, --predict_scope is a no-op for them, they
// train without shared params exactly as before this feature existed.
//

include { FUNANNOTATE_PREDICT }                           from '../../modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf'
include { FUNANNOTATE_PREDICT as FUNANNOTATE_PREDICT_SIB } from '../../modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf'
include { BACKFILL_ABINITIO_PARAMS }                       from '../../modules/funannotate/predict/BACKFILL_ABINITIO_PARAMS/main.nf'

include { gbkResult; staleRnaseq; staleGenome; sharedParamsJsonFor; staleSharedParams } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_PREDICTION {
    take:
    predict_input_ch       // tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
    abinitioReuseMap       // out -> [species, reuse_eligible, is_representative], loaded offline
    forceIndependentSet    // species that always train independently

    main:
    def predictScope  = ((params.predict_scope ?: 'all') as String).toLowerCase()
    def allowFallback = (params.allow_independent_fallback ?: false).toString().toBoolean()

    // Classify once, up front. is_rep/eligible come from the offline-loaded CSV
    // map (a plain synchronous read, not a channel), so this is cheap per item.
    def classified = predict_input_ch
        .map { out, asmid, sp, st, lt, bl, hl, tt, gfa ->
            def assignment = abinitioReuseMap[out]
            def is_rep     = assignment?.is_representative ?: false
            def eligible   = (assignment?.reuse_eligible ?: false) && !forceIndependentSet.contains(sp)
            tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, is_rep, eligible)
        }

    def branched = classified.branch {
        representative:   it[9]
        eligible_sibling:  it[10]
        independent:       true
    }

    // Representatives never use -p on themselves -- they're what gets shared.
    def rep_todo = branched.representative
        .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _is_rep, _elig ->
            tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, '')
        }
        .filter { out, a, sp, _st, _lt, _bl, _hl, _tt, _gfa, _shared_json ->
            gbkResult("${params.target}/${out}/predict_results", out as String) == null ||
                staleRnaseq(out as String, sp as String) ||
                staleGenome(out as String, a as String)
        }

    def metadataOut
    if (predictScope == 'representative_only') {
        FUNANNOTATE_PREDICT(rep_todo)

        // Backfill every representative predicted (or already up to date) in this
        // run. Batched into groups of ~100 (see BACKFILL_ABINITIO_PARAMS/main.nf) --
        // one SLURM job per batch instead of one per representative.
        def freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
            .map { out, _a, sp, _st, _lt, _bl, _hl, _tt -> tuple(sp.toString(), out.toString()) }
            .collate(100)
            .map { batch -> tuple(batch.hashCode(), batch) }
        BACKFILL_ABINITIO_PARAMS(freshBackfillInput)

        metadataOut = FUNANNOTATE_PREDICT.out.metadata
    } else {
        // Species with no eligible/representative assignment at all (singleton
        // species, ANI below threshold, or run_ani_reuse off) train independently,
        // completely ungated -- they were never going to touch the shared store.
        def indep_todo = branched.independent
            .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _is_rep, _elig ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, '')
            }
            .filter { out, a, sp, _st, _lt, _bl, _hl, _tt, _gfa, _shared_json ->
                gbkResult("${params.target}/${out}/predict_results", out as String) == null ||
                    staleRnaseq(out as String, sp as String) ||
                    staleGenome(out as String, a as String)
            }

        FUNANNOTATE_PREDICT(rep_todo.mix(indep_todo))

        // Backfill every representative predicted in THIS run. Batched into
        // groups of ~100 (see BACKFILL_ABINITIO_PARAMS/main.nf) -- one SLURM
        // job per batch instead of one per representative.
        def freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
            .map { out, _a, sp, _st, _lt, _bl, _hl, _tt -> tuple(out.toString(), sp.toString()) }
            .filter { out, _sp -> abinitioReuseMap[out]?.is_representative ?: false }
            .map { out, sp -> tuple(sp, out) }
            .collate(100)
            .map { batch -> tuple(batch.hashCode(), batch) }
        BACKFILL_ABINITIO_PARAMS(freshBackfillInput)

        // .out.done emits one (items, marker) pair per BATCH, not per species --
        // flatten items (each a [species, rep_out] pair) back to one species per
        // emission so downstream stays the same shape it was pre-batching.
        def freshSpeciesCh = BACKFILL_ABINITIO_PARAMS.out.done
            .flatMap { items, _marker -> items.collect { sp, _rep_out -> sp.toString() } }

        // Species whose shared store already existed BEFORE this run started --
        // known synchronously from the CSV plus one filesystem check per species,
        // no channel needed for this half.
        def preexistingSpecies = abinitioReuseMap.values()
            .findAll { it.reuse_eligible }
            .collect { it.species }
            .unique()
            .findAll { sp -> sharedParamsJsonFor(sp) != null }

        // Species with shared params available *by the time this run needs them*:
        // pre-existing ones (known now) mixed with this run's own freshly-
        // backfilled ones (a real channel dependency on this run's own
        // representative predict+backfill), collected into a single Set once
        // that stream closes.
        //
        // NOT join(by:0, remainder:true): probed directly (see
        // .living/learnings.md) and confirmed it only satisfies ONE left-side
        // item per key -- with two eligible siblings of the same species, one
        // matched and the other was silently treated as a genuine non-match
        // and routed to the blocked report, which is wrong (it just hadn't
        // been the first one Nextflow happened to pair). A species key here
        // is inherently one-to-many (many eligible siblings can share it), so
        // join()'s implicit one-match-per-key assumption is the wrong tool.
        //
        // combine() with a *collected* Set as the broadcast value instead:
        // probed directly too, confirmed it broadcasts the whole Set as one
        // atomic tuple element to every left-side item, duplicate keys
        // included -- unlike the collectFile()-on-a-list and combine()-of-two-
        // collected-lists footguns already on file, this is combine()'s
        // documented safe shape (one real per-item channel + one genuinely
        // scalar-emission value on the other side).
        def availableSpeciesSet = Channel.fromList(preexistingSpecies)
            .mix(freshSpeciesCh)
            .unique()
            .collect()
            .map { it as Set }
            .ifEmpty([] as Set)

        def eligibleSiblingKeyed = branched.eligible_sibling
            .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _is_rep, _elig ->
                tuple(sp.toString(), tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa))
            }

        def gated = eligibleSiblingKeyed.combine(availableSpeciesSet)

        def readyRows = gated
            .filter { sp, _row, availSet -> availSet.contains(sp) }
            .map { sp, row, _availSet ->
                tuple(row[0], row[1], sp, row[3], row[4], row[5], row[6], row[7], row[8],
                      sharedParamsJsonFor(sp)?.toString() ?: '')
            }

        def blockedRows = gated
            .filter { sp, _row, availSet -> !availSet.contains(sp) }
            .map { _sp, row, _availSet -> row }

        def sibling_todo
        if (allowFallback) {
            def fallbackRows = blockedRows
                .map { out, asmid, sp, st, lt, bl, hl, tt, gfa ->
                    log.warn "predict: ${out} (${sp}) — shared ab-initio parameters not " +
                        "available; training independently (--allow_independent_fallback)"
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, '')
                }
            sibling_todo = readyRows.mix(fallbackRows)
        } else {
            blockedRows
                .map { out, asmid, sp, _st, _lt, _bl, _hl, _tt, _gfa ->
                    "${sp}\t${out}\t${asmid}\tshared_ab_initio_params_not_available"
                }
                .collectFile(name: 'predict_blocked_awaiting_representative.tsv',
                             storeDir: params.target, newLine: true, sort: true,
                             seed: "species\tout\tasmid\treason")
            sibling_todo = readyRows
        }

        def sibling_predict_todo = sibling_todo
            .filter { out, a, sp, _st, _lt, _bl, _hl, _tt, _gfa, shared_params_json ->
                gbkResult("${params.target}/${out}/predict_results", out as String) == null ||
                    staleRnaseq(out as String, sp as String) ||
                    staleGenome(out as String, a as String) ||
                    staleSharedParams(out as String, shared_params_json ? file(shared_params_json as String) : null)
            }

        FUNANNOTATE_PREDICT_SIB(sibling_predict_todo)

        metadataOut = FUNANNOTATE_PREDICT.out.metadata.mix(FUNANNOTATE_PREDICT_SIB.out.metadata)
    }

    emit:
    metadata = metadataOut
}
