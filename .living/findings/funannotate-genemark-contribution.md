---
topic: funannotate-genemark-contribution
description: How much GeneMark-ES's ab-initio gene models actually change the final EVM-consensus gene set in funannotate predict, and whether that varies by genome.
created: 2026-08-12
last_updated: 2026-08-12
status: active
---

# funannotate predict: GeneMark-ES contribution to the final gene set

> Opened 2026-08-12 while validating the `ghcr.io/nextgenusfs/funannotate` rust
> container (see `.living/learnings.md`, 2026-08-12 container smoke-test entry)
> — GeneMark is deliberately absent from that image (license restrictions
> forbid redistribution), motivating a possible standalone `GENEMARK_RUN`
> Nextflow task. Before investing in that split, measured whether GeneMark is
> worth keeping at all.

## Findings

### F-008 — GeneMark-ES's contribution to the final EVM-consensus gene set is highly genome-dependent, from negligible to dominant

**Status:** preliminary (2026-08-12, n=3 genomes — deterministic paired rerun, not a sampled measurement)

**Claim:** Reran `funannotate predict` for 3 already-trained genomes, once at
production weights (baseline) and once with `-w genemark:0` (nogenemark),
against the identical existing PASA training data, then compared final gene
sets by genomic coordinate (`bedtools intersect`, not just counts — EVM can
substitute a GeneMark call for what would otherwise be an Augustus/snap call
at the *same* locus, which a count-only diff would miss):

| Genome | Assembly | Baseline final genes | Nogenemark final genes | Genes only in baseline (GeneMark-attributable) | Genes only in nogenemark |
|---|---|---|---|---|---|
| Saccharomyces_kudriavzevii_IFO10991 | 9.7 Mb, 1,145 contigs (fragmented draft) | 882 | 328 | **603 (68%)** | 47 |
| Kluyveromyces_marxianus_YG-4 | — | 2,632 | 2,256 | **381 (14%)** | 3 |
| Penicillium_citrinum_NRRL_1841 | — | 11,198 | 11,363 | 44 (0.4%) | 152 |

GeneMark's contribution ranges from **dominant** (68% of *S. kudriavzevii*'s
final genes exist only because GeneMark supplied them — nogenemark loses
nearly two-thirds of the annotation) to **substantial** (14% for
*K. marxianus*) to **negligible/net-neutral** (Penicillium actually gained a
few more genes without GeneMark, net +165, likely EVM redistributing
consensus weight onto Augustus/PASA/snap rather than a real improvement).

The dominant case (*S. kudriavzevii*) is the most fragmented, smallest
assembly of the three — consistent with PASA/RNA-seq training evidence being
sparser on a fragmented draft, leaving GeneMark's self-trained ab-initio
model as the largest single evidence source feeding EVM. This baseline rerun
reproduced 882 final genes vs. the original production run's 878 (same
genome, same training data, different day) — near-identical, validating the
A/B setup's fidelity to production.

**Evidence:** `analysis/genemark_es_contribution/` —
`GENEMARK_ES_CONTRIBUTION.md`, `outputs/genemark_contribution_summary.csv`,
`scripts/compare_results.py`, `scripts/run_predict_variant.sbatch`. Raw
per-source model counts parsed from each run's `Summary of gene models: {...}`
log line; final gene sets diffed by `bedtools intersect -f 0.9`.

**Implications:** GeneMark is not a low-value/droppable ab-initio source in
general — its value is genome-dependent, and for at least one of three
sampled genomes it was carrying the majority of the final annotation. This
argues **for** building the standalone `GENEMARK_RUN` Nextflow task (rather
than deprioritizing GeneMark or accepting `--auto-skip-genemark` silently
degrading small/fragmented genomes) so container-based Trinity/PASA/EVM
adoption doesn't come at the cost of losing GeneMark's biggest wins on
exactly the genomes (small, fragmented, RNA-seq-poor) where it matters most.
n=3 is not enough to generalize the *magnitude* pattern (e.g. "fragmented →
GeneMark dominant") beyond this sample — worth revisiting with more genomes
spanning the small/fragmented ↔ large/complete axis before treating that
correlation as established.

**Related**: [[funannotate-train-performance]] (F-004, same predict/train
pipeline); `.living/learnings.md` 2026-08-12 container smoke-test entry
(motivating context); `.living/learnings.md` 2026-08-12 `AUGUSTUS_CONFIG_PATH`
basename gotcha (hit while building this analysis).

**Tags**: funannotate, genemark, evm, predict, ab-initio, ab-test, container-migration, fragmented-genome
