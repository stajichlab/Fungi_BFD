---
topic: funannotate-train-performance
description: Where funannotate train wall-clock time actually goes, and which upstream (funannotate-live/PASApipeline) optimization work is worth prioritizing.
created: 2026-07-23
last_updated: 2026-07-23
status: active
---

# funannotate train Performance

> Opened 2026-07-23 while triaging whether to prioritize Trinity-GG threading fixes
> (`funannotate-live`, upstream issues
> [#1178](https://github.com/nextgenusfs/funannotate/issues/1178)/
> [#1179](https://github.com/nextgenusfs/funannotate/issues/1179)) or PASA
> parallelization work (`PASApipeline`) for pipeline throughput.

## Findings

### F-004 — PASA (`Launch_PASA_pipeline.pl`), not Trinity-GG, dominates `funannotate train` wall-clock

**Status:** established (2026-07-23, n=400 random sample of 5,471 `funannotate-train.log` files)

**Claim:** For a typical per-strain `funannotate train` run, PASA's own
`Launch_PASA_pipeline.pl` invocation is a median **84.8%** (p90 94.4%) of the
*entire* run's wall time — median PASA time 73.2 min out of a median 73.8 min
total run. Trinity-GG genome-guided assembly (the stage the upstream
Butterfly-pool-sizing fix, #1178/PR #1180, targets) **does not run at all** in
336/336 of the PASA-dominant sampled logs, because per-strain training calls
pass `--trinity <shared_fasta>` (built once per species by `RNASEQ_PREPARE`),
which makes `funannotate train` skip its internal hisat2-align +
Trinity-GG-assembly path entirely. Only 4/400 sampled runs (1%) use an
older/different hisat2+StringTie code path without shared Trinity, and even
there PASA still consumes 37–254 min in absolute terms. Second-largest cost,
a distant second: TransDecoder training-set extraction
(`pasa_asmbls_to_training_set.dbi`) at 6.4% of sampled wall-clock.

**Evidence:** `analysis/funannotate_train_stage_timing/` —
`FUNANNOTATE_TRAIN_STAGE_TIMING.md`, `outputs/stage_summary.csv`,
`outputs/per_run_summary.csv`, `scripts/profile_train_stage_timing.py`. PASA
timing is directly measured (both its start and end lines are always
explicitly logged), independent of the stage-attribution methodology used for
other stages.

**Implications:** Prioritize the PASA parallelization gaps already identified
by code review (unthreaded `assign_clusters_by_gene_intergene_overlap.dbi` /
`assign_clusters_by_stringent_alignment_overlap.dbi`, unbatched
`subcluster_builder.dbi`, un-threaded TransDecoder call in
`pasa_asmbls_to_training_set.dbi`) over further Trinity-GG work — they sit
inside the ~80-85%-of-wall-time stage. The already-opened Trinity-GG fix
(#1178/PR #1180) remains a correct bugfix but has negligible real-world
throughput impact on this training corpus, since that code path essentially
never executes for per-strain runs.

**Related**: [[nextflow_memory_profile]] (FUNANNOTATE_TRAIN cpu/memory sizing);
decision log entry "Bump FUNANNOTATE_TRAIN attempt-1 cpus 2→4" (`.living/decisions.md`,
2026-07-23) — that cpu bump still helps the (now known to be tiny) Trinity-GG-running
minority of runs and the aligner sub-step inside PASA, but is not the primary lever
for overall training throughput.

**Tags**: funannotate, pasa, trinity, throughput, log-parsing, upstream-bugfix, funannotate-train
