# Consolidate FUNANNOTATE_PREDICTION / FUNANNOTATE_ANNOTATION stale-check duplication

**Status**: scoped, shelved (2026-09-09) — not blocking, revisit only on trigger conditions below.

## Background

`FUNANNOTATE_PREDICTION.nf` classifies every genome in `predict_input_ch` into
"needs `FUNANNOTATE_PREDICT`" (via `rep_todo`/`indep_todo`/`sibling_predict_todo`,
each filtering on `gbkResult()==null || staleRnaseq() || staleGenome()`, the
sibling one also `staleSharedParams()`) vs. implicitly "already resolved" — but
only emits the "needs predict" side (`metadata`). The "already resolved" side is
computed but discarded.

`FUNANNOTATE_ANNOTATION.nf`'s `postpredict` (lines ~38-65) independently
rebuilds that same "already resolved" set from scratch: re-reads the entire
`samples.csv` (23,683 rows), reapplies the same `taxonFilter`/`asmidFilter`/
`suppressFilter`/`n_test` closures (confirmed to be the literal same closure
instances, passed to both from `workflows/funannotate.nf`), then filters on
`gbkResult()!=null && !staleRnaseq() && !staleGenome()`.

This is real, confirmed redundant computation for the subset of genomes that
are simultaneously (a) already predicted from a prior run, (b) currently stale,
and (c) part of this run's `predict_input_ch` — such a genome gets its
staleness computed and logged in both places. The `FUNANNOTATE_ANNOTATION.nf`
header comment claiming the predict/postpredict split is "disjoint by
construction" is currently **false**: a timing gap (postpredict evaluates at
launch; PREDICTION evaluates after this run's own RNASEQ/training may have
changed the trinity/genome files) and a missing `staleSharedParams` check in
`postpredict` both let a genome land in both `metadata_ch` and `postpredict`
under real conditions.

**This was investigated as the suspected cause of a measured ~2,200-line
duplicate-log-line symptom — it was NOT the cause.** The actual cause (two
independent `.filter{}`s in `FUNANNOTATE_RNASEQ.nf`'s `train_todo`/`train_done`)
was found and fixed in commit `191ab3a` (see `.living/decisions.md` and
`.living/learnings.md`, both 2026-09-09). This item is the leftover, genuinely
separate architectural redundancy — real, but smaller and higher-risk to fix
than it first appeared.

## Proposed fix, if picked back up

1. **`FUNANNOTATE_PREDICTION.nf`**: convert each of the three "needs predict"
   filters into a `.branch{}` capturing both sides in one evaluation (mirroring
   the `train_branched` pattern from `FUNANNOTATE_RNASEQ.nf`, so this doesn't
   introduce a *new* double-eval bug):
   - `rep_split = rep_rows.branch { todo: P1(it); resolved: true }` where
     `P1 = gbkResult==null || staleRnaseq || staleGenome` (currently lines ~89-93)
   - `indep_split` — same `P1`, `all` scope only (currently ~156-160)
   - `sib_split = sibling_todo.branch { todo: P2(it); resolved: true }` where
     `P2 = P1 || staleSharedParams` — sits after the `availableSpeciesSet` gate
     (currently ~297-303), so the "resolved" signal for a sibling is delayed
     behind every representative's predict finishing in this run, same as the
     `todo` side already is today
   - `blocked_resolved = blockedRows.filter { !P1(it) }` — needed only when
     `!allow_independent_fallback`, otherwise a blocked sibling with a current
     GBK silently drops out of annotation entirely
   - `--predict_scope representative_only` mode needs its own explicit path,
     since `eligible_sibling`/`independent` currently have no consumer at all
     in that mode
   - Emit `already_resolved = mix of all "resolved" branches`, mapped to
     `(out, asmid, sp, st, lt, bl, hl, tt)`
2. **`workflows/funannotate.nf`**: wire `FUNANNOTATE_PREDICTION.out.already_resolved`
   into `FUNANNOTATE_ANNOTATION`.
3. **`FUNANNOTATE_ANNOTATION.nf`**: replace the `samples.csv` rescan
   (`postpredict`, lines ~38-65) with consuming `already_resolved` directly.
   **Requires an explicit policy decision**: `predict_input_ch`'s genome scope
   is a strict subset of what the current rescan covers — additional drop
   points exist upstream (missing source genome, failed cleaning, `--only_clean`,
   inner-join channel combines in RNASEQ, and `FUNANNOTATE_TRAIN` failures/skips).
   One of these — a genome with a current GBK but unresolved training state —
   isn't logged anywhere, so its real-world frequency is unknown. A naive
   substitution could silently narrow what gets annotated. Either accept the
   narrowing (and monitor before/after annotate row counts on a real run) or
   add a supplementary check for that specific case.
4. Fix the stale "disjoint by construction" comments in both files once this
   actually makes the claim true.

## Why shelved (2026-09-09)

- Not the cause of the measured symptom (already fixed separately).
- Unmeasured performance benefit — likely modest (tens of thousands of cheap
  `stat()` calls per launch; no evidence this costs more than low single-digit
  seconds of wall-clock time).
- Real implementation risk: scope-narrowing hazard (item 3 above),
  `representative_only`-mode gap, and a sibling-annotation delay if built the
  straightforward way.
- The correctness angle (today's false "disjoint by construction" claim) is
  the more compelling reason to eventually do this, independent of performance.

## Revisit triggers

- (a) Actual slow DAG-construction/pipeline-launch time is observed and traced
  specifically to this rescan (not just assumed).
- (b) A concrete double-annotation or missing-annotation case is observed in
  practice (a genome processed by `FUNANNOTATE_ANNOTATE`/etc. twice, or a
  genome with a current GBK that never gets annotated).
