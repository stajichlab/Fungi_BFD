---
topic: annotation-failure-modes
description: Why genomes fail or stall in the funannotate annotation pipeline (too-small/fragmented assemblies, training-model shortfall, oversized-genome timeouts) and how to detect/gate them.
created: 2026-06-25
last_updated: 2026-06-25
status: active
---

# Annotation Failure Modes

> Opened 2026-06-24/25 from the funannotate predict-failure screen
> (`analysis/funannotate_model_failures/`). Tracks the empirical failure modes of
> the genome-annotation pipeline and the detection/gating rules that address them.

## Findings

### F-003 — funannotate predict failures split into two disjoint causes; only the "too-small" one is assembly-stat-detectable

**Status:** established (2026-06-24, n=8,075 genomes)

**Claim:** Of 8,075 genomes scanned (deduped from 8,082 `funannotate-predict.log`
files), **8,054 completed**. Non-completions split into two physically distinct
classes:
1. **Too small / fragmented (unfixable property of the assembly):** 9 genomes hit the
   authoritative line `Not enough gene models N to train Augustus (30 required),
   exiting` (+ 2 `AUGUSTUS training failed`, 1 GeneMark `CMD ERROR`, 2
   `0 valid BUSCO predictions`). Every one is small **and** fragmented
   (127 kb–15.7 Mb; rusts that should be 100 Mb+ sitting at 2–7 Mb). A screen rule of
   **small AND fragmented** (assembled bp < ~16 Mb AND (N50 < ~10–20 kb OR contigs >
   ~1000)) flags **9/9** of these with **zero** false-skips of complete small genomes.
2. **Oversized / slow (resource problem):** 6 `INCOMPLETE_unknown` are the opposite —
   large genomes (0.5–1.3 Gb rusts; good 18–34 Mb assemblies with N50 up to 47 Mb)
   that time out / OOM / are still running. The size screen correctly does **not**
   flag them; they need walltime/memory, not a size filter.

**Evidence:** `analysis/funannotate_model_failures/` — `crosstab_summary.txt`,
`failed_genomes_annotated.tsv`, `parse_predict_failures.py`. Metrics for the failing
genomes were computed on the fly with `seqkit` because they are absent from the (stale)
asm_stats table; the portable N50 calc was verified identical to seqkit.

**Implications:**
- A too-small **pre-screen** is valid and high-precision, but **requires a current
  asm_stats table** — the existing `results/genome_stats/asm_stats/` (4,094 rows) was
  built on an earlier genome set and contained none of the 21 failing genomes.
- A flat genome-size cutoff is wrong: *Malassezia* (~7–9 Mb, complete) and Ashbya
  yeasts (~9 Mb, chromosome-level) annotate fine; the AND-with-fragmentation rule is
  what makes the screen safe.
- Implemented as an in-pipeline guard (pre-flight + post-predict catch) in
  `funannotate.nf`; the 9 genomes feed the manual `suppress.txt` curation (7 newly
  added, 2 pre-existing).

**Safety check (false-skip audit):** Ashbya/Eremothecium, complete Microsporidia,
Rozella (11.3 Mb), Paramicrosporidium (7.2 Mb / N50 70 kb) all pass the gate (frag=0)
→ never skipped.

### Evidence Ledger
| Date | Run/Session | Dataset | Result | Direction |
|------|-------------|---------|--------|-----------|
| 2026-06-24 | predict-failure screen | 8,082 funannotate-predict.log; samples.csv; input_clean_genomes/ | 8,054 completed; 9 too-few-models (all small+fragmented); 6 oversized-incomplete | too-small failures fully separable by asm-stats; oversized failures are not |

### Open Questions
- After `CALC_ASM_STATS` is re-run, does the merged `asm_stats.tsv.gz` flag the same 9
  (and no complete genomes)? (validates the pre-screen end-to-end)
- Should the oversized-genome class get its own handling (masking/longer walltime) —
  e.g. `GCA_025201825.1_Phapa1` is already in `suppress.txt` as "Too big - need to mask"?
