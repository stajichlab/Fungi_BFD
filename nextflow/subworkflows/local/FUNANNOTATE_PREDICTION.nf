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
// GENEMARK_RUN (nextflow/docs/GENEMARK_RUN_DESIGN.md): standalone GeneMark
// step ahead of predict, so predict can move onto the rust container. Two
// invocations needed -- rep/indep (mutually exclusive with the SIB one
// within a single execution, since the rep/indep call happens once per
// predictScope branch and only one branch's Groovy code path ever actually
// executes) share the plain `GENEMARK_RUN` name; the eligible-sibling
// invocation, which genuinely coexists with it in --predict_scope all, needs
// its own alias -- same reasoning as FUNANNOTATE_PREDICT/_SIB above.
include { GENEMARK_RUN }                                  from '../../modules/funannotate/predict/GENEMARK_RUN/main.nf'
include { GENEMARK_RUN as GENEMARK_RUN_SIB }               from '../../modules/funannotate/predict/GENEMARK_RUN/main.nf'

include { gbkResult; staleRnaseq; staleGenome; sharedParamsJsonFor; staleSharedParams; sharedGenemarkModFor; trainingTranscriptBamFor } from '../../modules/funannotate/utils.nf'

workflow FUNANNOTATE_PREDICTION {
    take:
    predict_input_ch             // tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
    abinitioReuseMap             // out -> [species, reuse_eligible, is_representative], loaded offline
    forceIndependentSet          // species that always train independently (whole ab-initio bundle)
    forceIndependentGenemarkSet  // out values that always run GeneMark independently, finer-grained
                                  // than forceIndependentSet -- AUGUSTUS/SNAP reuse can stay on for a
                                  // strain while forcing GeneMark specifically to retrain for it

    main:
    def predictScope  = ((params.predict_scope ?: 'all') as String).toLowerCase()
    def allowFallback = (params.allow_independent_fallback ?: false).toString().toBoolean()
    def runGenemark   = (params.run_genemark ?: true).toString().toBoolean()
    // ET mode (nextflow/docs/GENEMARK_RUN_DESIGN.md's "ET mode" section):
    // recipe validated real end-to-end (10,776 genes), wired here. GENEMARK_RUN
    // itself falls back to ES per-genome when training_bam is empty (no RNA-seq
    // for that genome) -- this param picks which mode to REQUEST, not a
    // guarantee every genome actually gets it.
    def genemarkMode  = ((params.genemark_mode ?: 'ES') as String).toUpperCase()

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
        // GENEMARK_RUN ahead of predict (see include block comment for why this
        // single call, un-aliased, is safe here: only one predictScope branch's
        // Groovy code path ever actually executes). Representatives never reuse
        // a shared .mod -- they're what DEFINES the store -- so shared_mod is
        // always '' and GENEMARK_RUN always takes the fresh --ES path,
        // guaranteeing every row here gets both a gtf AND a mod. If
        // run_genemark=false, skip the process entirely and fall back to
        // predict's own --auto-skip-genemark degradation (empty genemark_gtf).
        def rep_with_gtf
        if (runGenemark) {
            def rep_genemark_in = rep_todo.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _shared_json ->
                def forceIndep = forceIndependentGenemarkSet.contains(out as String) ? 'true' : 'false'
                tuple(out, asmid, sp, st, gfa, tt, genemarkMode, trainingTranscriptBamFor(out as String), forceIndep, '')
            }
            GENEMARK_RUN(rep_genemark_in)
            rep_with_gtf = rep_todo.join(GENEMARK_RUN.out.gtf)
                .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf)
                }
        } else {
            rep_with_gtf = rep_todo.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, '')
            }
        }

        FUNANNOTATE_PREDICT(rep_with_gtf)

        // Backfill every representative predicted (or already up to date) in this
        // run. Batched into groups of ~100 (see BACKFILL_ABINITIO_PARAMS/main.nf) --
        // one SLURM job per batch instead of one per representative. Join in
        // GENEMARK_RUN's fresh .mod by `out` (1:1, safe -- not the species-keyed
        // many-to-one shape the availableSpeciesSet join below has to work
        // around) when run_genemark is on; empty genemark_mod field otherwise so
        // backfill_abinitio_params.py's augustus/snap backfill still runs.
        def freshBackfillInput
        if (runGenemark) {
            freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
                .map { out, _a, sp, _st, _lt, _bl, _hl, _tt -> tuple(out.toString(), sp.toString()) }
                .join(GENEMARK_RUN.out.mod.map { out, _sp, mod -> tuple(out.toString(), mod.toString()) })
                .map { out, sp, mod -> tuple(sp, out, mod) }
                .collate(100)
                .map { batch -> tuple(batch.hashCode(), batch) }
        } else {
            freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
                .map { out, _a, sp, _st, _lt, _bl, _hl, _tt -> tuple(sp.toString(), out.toString(), '') }
                .collate(100)
                .map { batch -> tuple(batch.hashCode(), batch) }
        }
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

        // rep_todo and indep_todo both always pass shared_mod='' to GENEMARK_RUN
        // (representatives define the store, independents were never eligible
        // to reuse it) -- functionally identical fresh-ES paths, so one shared
        // GENEMARK_RUN call on the mixed channel is correct, not just
        // convenient (avoids needing a 3rd process alias). Join back onto the
        // SAME mixed channel before predict; join-by-`out` doesn't care which
        // original branch a row came from.
        def rep_and_indep = rep_todo.mix(indep_todo)
        def rep_and_indep_with_gtf
        if (runGenemark) {
            def genemark_in = rep_and_indep.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _shared_json ->
                def forceIndep = forceIndependentGenemarkSet.contains(out as String) ? 'true' : 'false'
                tuple(out, asmid, sp, st, gfa, tt, genemarkMode, trainingTranscriptBamFor(out as String), forceIndep, '')
            }
            GENEMARK_RUN(genemark_in)
            rep_and_indep_with_gtf = rep_and_indep.join(GENEMARK_RUN.out.gtf)
                .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf)
                }
        } else {
            rep_and_indep_with_gtf = rep_and_indep.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, '')
            }
        }

        FUNANNOTATE_PREDICT(rep_and_indep_with_gtf)

        // Backfill every representative predicted in THIS run. Batched into
        // groups of ~100 (see BACKFILL_ABINITIO_PARAMS/main.nf) -- one SLURM
        // job per batch instead of one per representative. Representatives
        // always get a fresh .mod from GENEMARK_RUN above (never reuse), so
        // this join is 1:1/safe, same argument as the representative_only
        // branch's equivalent join.
        def freshBackfillInput
        if (runGenemark) {
            freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
                .map { out, _a, sp, _st, _lt, _bl, _hl, _tt -> tuple(out.toString(), sp.toString()) }
                .filter { out, _sp -> abinitioReuseMap[out]?.is_representative ?: false }
                .join(GENEMARK_RUN.out.mod.map { out, _sp, mod -> tuple(out.toString(), mod.toString()) })
                .map { out, sp, mod -> tuple(sp, out, mod) }
                .collate(100)
                .map { batch -> tuple(batch.hashCode(), batch) }
        } else {
            freshBackfillInput = FUNANNOTATE_PREDICT.out.metadata
                .map { out, _a, sp, _st, _lt, _bl, _hl, _tt -> tuple(out.toString(), sp.toString()) }
                .filter { out, _sp -> abinitioReuseMap[out]?.is_representative ?: false }
                .map { out, sp -> tuple(sp, out, '') }
                .collate(100)
                .map { batch -> tuple(batch.hashCode(), batch) }
        }
        BACKFILL_ABINITIO_PARAMS(freshBackfillInput)

        // .out.done emits one (items, marker) pair per BATCH, not per species --
        // flatten items (each a [species, rep_out, genemark_mod] triple) back to
        // one species per emission so downstream stays the same shape it was
        // pre-batching.
        def freshSpeciesCh = BACKFILL_ABINITIO_PARAMS.out.done
            .flatMap { items, _marker -> items.collect { sp, _rep_out, _mod -> sp.toString() } }

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

        // Eligible siblings: GENEMARK_RUN must read from sibling_predict_todo
        // (post availability-gating/readyRows), not branched.eligible_sibling
        // directly -- shared_mod is only resolvable once availableSpeciesSet has
        // been computed, and resolving it earlier would mean duplicating
        // readyRows'/blockedRows'/allowFallback's gating logic a second time.
        // No-op in --predict_scope representative_only (siblings are never
        // predicted there at all -- branched.eligible_sibling has no consumer
        // in that branch), which is why this lives only inside the `else`.
        def sibling_predict_with_gtf
        if (runGenemark) {
            def sib_genemark_in = sibling_predict_todo.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _shared_json ->
                def forceIndep = forceIndependentGenemarkSet.contains(out as String) ? 'true' : 'false'
                def sharedMod  = sharedGenemarkModFor(sp as String)?.toString() ?: ''
                tuple(out, asmid, sp, st, gfa, tt, genemarkMode, trainingTranscriptBamFor(out as String), forceIndep, sharedMod)
            }
            GENEMARK_RUN_SIB(sib_genemark_in)
            sibling_predict_with_gtf = sibling_predict_todo.join(GENEMARK_RUN_SIB.out.gtf)
                .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf)
                }
        } else {
            sibling_predict_with_gtf = sibling_predict_todo.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json ->
                tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, '')
            }
        }
        // A sibling's own GENEMARK_RUN_SIB.out.mod (only emitted when
        // force_independent=true forced a fresh train, or the very rare case
        // a sibling ran before its species had any shared mod at all) is
        // deliberately never wired into backfill -- only representatives
        // define the shared per-species store; a forced-independent sibling's
        // one-off GeneMark model must not silently become that.

        FUNANNOTATE_PREDICT_SIB(sibling_predict_with_gtf)

        metadataOut = FUNANNOTATE_PREDICT.out.metadata.mix(FUNANNOTATE_PREDICT_SIB.out.metadata)
    }

    emit:
    metadata = metadataOut
}
