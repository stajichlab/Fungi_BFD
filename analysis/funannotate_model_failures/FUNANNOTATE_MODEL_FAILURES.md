# Funannotate "Not enough training models" failures — detection & guard

**Question:** identify (esp. basidiomycota) genomes too small to annotate, that fail
`funannotate predict` because it cannot assemble the 30 gene models needed to train
Augustus. Detect both (a) reactively from logs and (b) proactively from assembly stats,
and add an in-pipeline guard.

## Ground-truth failure signature
`funannotate predict` aborts with the verbatim line:
```
Not enough gene models N to train Augustus (30 required), exiting
```
This is the authoritative "too few models" signal. Related non-completions:
`AUGUSTUS training failed, check logfiles`; `0 valid BUSCO predictions found`;
`CMD ERROR ... genemark`.

## What we did
1. **`parse_predict_failures.py`** — scans every persisted `funannotate-predict.log`
   (8,082 under `genome_annotation/` + `do_annotation/genome_annotation/`), classifies
   each genome's outcome, and joins to taxonomy + assembly stats.
   - Join key: log line 1 carries `--name <LOCUSTAG>`; LOCUSTAG is **unique** in
     `samples.csv` (0 dups in 22,929 rows) → ASMID, PHYLUM. (Do NOT parse ASMID from the
     `-i` path: gzipped genomes are inflated to `genome_input.fa`, losing the ASMID.)
   - Outputs: `predict_outcomes.tsv` (all genomes), `failed_genomes.tsv`,
     `crosstab_summary.txt`, and (added) `failed_genomes_annotated.tsv` with freshly
     computed seqkit metrics + screen flag.

## Results (run 2026-06-24)
- 8,075 genomes (deduped). Outcomes: **8,054 COMPLETED**, 9 `FAILED_too_few_models`,
  2 `FAILED_augustus_training`, 1 `FAILED_cmd_error`, 9 INCOMPLETE.
- **All 9 `FAILED_too_few_models` are tiny/fragmented** (127 kb – 15.7 Mb). Rusts
  (Uromyces, Puccinia) that should be 100 Mb+ sit at 2–7 Mb = grossly partial assemblies.
- Applying the screen rule (small **AND** fragmented) to fresh seqkit metrics:
  **100% of the too-few-models failures are flagged** (8 HIGH + 1 MED). The
  Augustus-training and cmd-error junk are also flagged HIGH.
- The 6 `INCOMPLETE_unknown` are the **opposite** problem: large genomes (0.5–1.3 Gb
  rusts; or good 18–34 Mb assemblies, N50 up to 47 Mb) timing out / OOM / still running.
  The screen correctly does **not** flag them — they need resources, not a size filter.

### Critical caveat discovered
The pre-existing `results/genome_stats/asm_stats/` table (4,094 genomes) was built on an
**earlier genome set and never recomputed** — *none* of the 21 failing genomes had an
asm_stats entry. So the asm_stats pre-screen can only be used after `CALC_ASM_STATS` is
re-run over the current `samples.csv`. The parser therefore computes metrics on the fly
from `input_clean_genomes/` for genomes missing from the table.

### Relationship to `suppress.txt`
Manual curation already lists several of these ("Too small (800kb)", "Too small for Yeast
genome", "All contigs too small"). 3 of the 21 detected failures are already suppressed.
**The guard automates what is currently manual.**

## In-pipeline guard (FUNANNOTATE_PREDICT)
Two layers, both in `nextflow/funannotate.nf`, tunable in `conf/profile_funannotate.config`:
1. **Pre-flight gate** (before predict): cheap contig stats from the input FASTA; if the
   assembly is **small AND fragmented** it is flagged to
   `${target}/predict_skipped_too_small.tsv` and skipped (touch `.predict.done`, exit 0)
   instead of burning hours then crashing. Requires **both** gates so complete small
   genomes are never false-skipped. Params:
   `predict_min_asm_bp=8000000`, `predict_frag_max_n50=10000`,
   `predict_frag_max_contigs=1000` (set `predict_min_asm_bp=0` to disable).
2. **Post-predict catch** (backstop): if predict produced no GBK, grep the log for
   `Not enough gene models .* to train Augustus`; if matched, flag + skip cleanly;
   otherwise still hard-fail so real errors surface.

### Safety against legitimately small genomes (verified 2026-06-24)
The **BOTH-gates** requirement is what protects compact genomes:
- **Ashbya / Eremothecium yeasts**: 8.9–9.7 Mb, 7–8 contigs, N50 ~1.2–1.5 Mb →
  small=0, frag=0 → **never skipped**.
- **Microsporidia** (already skipped upstream by the project): *complete* ones (e.g.
  Antonospora 3.2 Mb / 17 c / N50 184 kb) → frag=0 → **not skipped**; only fragmented
  sub-Mb drafts (975 kb / 305 c / N50 3 kb) trip both gates — and those genuinely fail
  funannotate anyway. A complete genome (high N50, few contigs) cannot trip the
  fragmentation gate, so the guard cannot false-skip it.

## Reproduce
```bash
bash analysis/funannotate_model_failures/run.sh
```
