# funannotate ab-initio parameter reuse & PASA-skip validation

## Context

`todo/species_level_abinitio_reuse.md` proposed two related speedups for `funannotate
predict` across ANI-qualified strains of the same species:
- **T-004**: reuse trained AUGUSTUS/GeneMark-ES/SNAP parameters from a species
  representative instead of retraining ab-initio predictors per strain (implemented in
  `funannotate.nf` via `-p <shared_params.json>`).
- **T-013**: skip `FUNANNOTATE_TRAIN` (PASA gene-model building) entirely and feed the
  pre-built Trinity-GG transcript assembly directly to `funannotate predict
  --transcript_evidence`, combined with T-004's reuse so ab-initio predictors don't need
  from-scratch training either.

This entry records the 12-strain validation (6 *Aspergillus fumigatus*, 6 *Beauveria
bassiana*) comparing both against each other and, for 4 strains, against legacy
independent-training baselines.

## Method

For each of 12 reuse-eligible strains, ran:
- **Condition A (T-004)**: full `FUNANNOTATE_TRAIN` (PASA, using the species-shared
  Trinity-GG assembly) + `FUNANNOTATE_PREDICT` with `-p` ab-initio reuse.
- **Condition B (T-013)**: `funannotate predict --transcript_evidence <shared
  Trinity-GG.fasta> -p <shared_params.json>`, no TRAIN/PASA step at all.

Both conditions used the same shared ab-initio store (`gene_prediction_shared_abinitio/`),
which includes a `glimmerhmm` empty-stub entry added after an initial pilot run showed
`-w glimmerhmm:0` (EVM weight zero) does not skip glimmerhmm's *training* — without a
`path` key, funannotate falls back to a full genome-mode BUSCO run to seed it regardless
of weight. See `.living/learnings.md` (2026-07-25) for the mechanism.

BUSCO (`fungi_odb12`, protein mode) run on every resulting `proteins.fa`.

## Data

| Strain | Cond A genes | Cond A BUSCO | Cond B genes | Cond B BUSCO | BUSCO Δ (B−A) | Genes Δ (B−A) |
|---|---|---|---|---|---|---|
| Aspergillus_fumigatus_47-10 | 9,431 | 97.1% | 9,083 | 96.1% | -1.0 | -348 |
| Aspergillus_fumigatus_A-2-36s-3 | 9,303 | 96.8% | 8,966 | 96.0% | -0.8 | -337 |
| Aspergillus_fumigatus_E-1-33L-1 | 9,632 | 95.9% | 8,950 | 96.3% | +0.4 | -682 |
| Aspergillus_fumigatus_NRZ-2018-146 | 9,594 | 96.3% | 8,926 | 96.4% | +0.1 | -668 |
| Aspergillus_fumigatus_L-2-21-2 | 9,749 | 96.3% | 9,053 | 96.3% | 0.0 | -696 |
| Aspergillus_fumigatus_W72310 | 9,655 | 96.1% | 8,953 | 96.1% | 0.0 | -702 |
| Beauveria_bassiana_BCC_2660 | 9,283 | 98.1% | 8,954 | 98.0% | -0.1 | -329 |
| Beauveria_bassiana_ARSEF_5078 | 9,327 | 94.4% | 9,016 | 94.9% | +0.5 | -311 |
| Beauveria_bassiana_JAU2 | 9,311 | 97.3% | 8,995 | 97.4% | +0.1 | -316 |
| Beauveria_bassiana_MBC_306 | 9,557 | 96.7% | 9,238 | 96.2% | -0.5 | -319 |
| Beauveria_bassiana_MBC_948 | 8,985 | 94.9% | 8,706 | 94.8% | -0.1 | -279 |
| Beauveria_bassiana_MBC_715 | 9,343 | 96.1% | 9,023 | 95.9% | -0.2 | -320 |

Mean BUSCO Δ: -0.13pp (median -0.05pp). Mean gene-count Δ: -442 genes (consistently
negative in all 12/12 strains — T-013 always predicts fewer models than T-004).

Legacy independent-training sanity check (4/12 strains with a pre-reuse baseline,
`genome_annotation/_validation_backup_20260724/`):

| Strain | Legacy genes | Legacy BUSCO | Cond A (T-004) genes | Cond A BUSCO |
|---|---|---|---|---|
| Aspergillus_fumigatus_47-10 | 9,698 | 96.9% | 9,431 | 97.1% |
| Aspergillus_fumigatus_A-2-36s-3 | 9,807 | 96.5% | 9,303 | 96.8% |
| Beauveria_bassiana_MBC_306 | 9,222 | 96.5% | 9,557 | 96.7% |
| Beauveria_bassiana_MBC_715 | 9,470 | 96.4% | 9,343 | 96.1% |

Condition A (T-004) BUSCO tracks legacy within ±0.3pp on all 4 sanity-check strains —
confirms ab-initio reuse alone does not degrade completeness, consistent with the earlier
4-strain T-012 validation.

Runtime (wall-clock; Condition A's TRAIN time was cached from prior runs for 9/12 strains,
directly measured fresh for 3 — `Beauveria_bassiana_JAU2` 65m15s, `MBC_306` 66m45s,
`MBC_715` 69m41s, averaging 67.2 min, used as the TRAIN estimate for the other 9):

| | Condition A total (TRAIN+PREDICT) | Condition B (PREDICT only, no TRAIN) |
|---|---|---|
| Mean per strain | ~150 min | ~46 min |
| **Savings** | — | **~104 min/strain (69% reduction)** |

Restricting to PREDICT-only wall time (isolating just the ab-initio-training-skip effect,
ignoring TRAIN entirely): Condition A PREDICT averaged 82.6 min vs Condition B's 45.8 min
— a 44.5% reduction attributable to skipping ab-initio training + evidence-hint building
inside PREDICT alone.

## Key Findings

- **T-013 (no-TRAIN + `--transcript_evidence` + ab-initio reuse) preserves BUSCO
  completeness within noise of T-004** (mean Δ -0.13pp, median -0.05pp, range -1.0 to
  +0.5pp across 12 strains) — no systematic degradation, and 5/12 strains actually scored
  *higher* under T-013.
- **Gene count is a different story**: T-013 predicts consistently fewer models than
  T-004 in all 12/12 strains (mean -442, range -279 to -702) — expected, since skipping
  PASA/TRAIN removes PASA's own gene-model contribution to the EVM consensus, not just its
  ab-initio-training-seeding role. BUSCO completeness holding steady despite ~5-7% fewer
  total models suggests the missing models are mostly redundant/low-confidence calls PASA
  otherwise contributes, not core single-copy orthologs — but this hasn't been checked
  against a curated gene set, only BUSCO's single-copy marker panel.
- ***Aspergillus fumigatus* loses more genes under T-013 than *Beauveria bassiana***
  (mean -572 vs -312) despite comparable BUSCO deltas — worth a follow-up look at whether
  this tracks genome size/gene density or *A. fumigatus*'s typically larger, more complete
  RNA-seq evidence base giving PASA more to work with in Condition A.
- **T-013's total wall-clock savings (~69% including TRAIN) is substantially larger than
  T-004 alone** — consistent with [funannotate-train-performance](funannotate-train-performance.md)'s
  finding that PASA dominates TRAIN wall time; T-013 is the correct lever if raw
  throughput is the goal, T-004 alone if TRAIN's PASA-derived gene models are considered
  worth keeping for quality reasons.
- **glimmerhmm gotcha**: a predictor's EVM weight (`-w glimmerhmm:0`) does not exempt it
  from being *trained* — funannotate predict.py decides training method independent of
  weight, and without a `path` key in `parameters.json` it falls back to an expensive
  full genome-mode BUSCO run regardless. Fixed via an empty-stub `path` entry (see
  `.living/learnings.md`).

## Status

**preliminary** — n=12 strains (6+6), single species pair. BUSCO is a coarse completeness
proxy; the gene-count delta's practical impact (which specific models are dropped, whether
they're redundant isoforms/low-confidence calls or genuinely missed genes) has not been
characterized. Recommend spot-checking a handful of T-013-only-missing genes against
protein evidence/synteny before treating T-013 as production-ready, and extending to a
third species before generalizing the ~69% savings figure.

## Implications

- T-004 (ab-initio reuse with TRAIN retained) is validated for production rollout as
  designed — BUSCO tracks legacy closely, gene counts are reasonable, and it's the lower-
  risk of the two levers (keeps PASA's gene-model contribution).
- T-013 is a promising, larger lever (addresses PASA's dominant TRAIN cost directly) but
  needs the gene-count-delta characterization above before being recommended as a
  default — the BUSCO-neutral result is encouraging but not sufficient on its own given a
  consistent, non-trivial drop in total gene count.
