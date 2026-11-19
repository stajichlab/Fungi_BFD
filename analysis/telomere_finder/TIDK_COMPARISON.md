# tidk vs. BFD telomere scripts (2026-11-18)

Question: does `tolkit/telomeric-identifier` (tidk 0.2.1, module `tidk/0.2.1`) do better than
`nextflow/bin/find_telomeres.py` + `denovo_telomere_repeats.py`?

## Method
Five genomes from `input/dna`, production parameters from `nextflow.config`
(`--both-ends`, exact mode). tidk: `search -s TTAGGG -w 1000` (an end is called when the
first/last 1 kb window has >=5 repeat copies) and `explore -l 6 -m 5 -x 12 --distance 10000`.
Scaffolds >=100 kb only. Script: `compare_tidk.py`. Raw output was in the session scratchpad
(not kept). No curated ground truth exists except *F. graminearum* CML3066 (4 chromosomes + mito).

## Ends called (scaffolds >=100 kb)
| Genome | Ends | find_telomeres | denovo | tidk search | tidk explore |
|---|---|---|---|---|---|
| F. graminearum (5 scaffolds) | 10 (8 real) | 7 | 7 | 7 | 4 |
| A. nidulans 1680 (8 scaffolds) | 16 | 10 | 14 | 14 | 7 |
| F. graminearum 03132 (16 scaffolds) | 32 | 0 | 4 | 0 | 0 |
| C. albicans 101 (10 scaffolds) | 20 | 0 | 15 | 0 | 0 |
| N. crassa FGSC_2489 (210 scaffolds, fragmented) | 228 | 10 | 19 | 11 | 34 |

## Findings
- tidk `search` needs a motif you supply. It then equals `find_telomeres.py` on clean cases and
  beats it on degraded ends (A. nidulans 14 vs 10). `denovo` matches tidk there (14).
- tidk `find` has no fungal clade (only Agaricales). `explore` is capped at 12-mers, so it
  misses the 23 bp *C. albicans* repeat, and it lists rotations of one repeat as separate hits.
- tidk gives windowed counts. It gives no tract coordinates, `terminal` flag, `distance_to_end`
  or flank sequence, which `telomeres.parquet` needs.
- `find_telomeres.py` misses on A. nidulans came from (a) the 10 bp `--terminal-tolerance`
  (e.g. JBEUPP010000004.1 5': tract starts 14 bp in; JBEUPP010000002.1 3': ~14 bp) and (b) exact
  matching broken by a single indel ~20 bp from the end (JBEUPP010000003.1 3').
  `--fuzzy --max-indel 1` raised A. nidulans tracts 10 -> 12; `--terminal-tolerance 30` -> 14.
  Fuzzy mode on *C. albicans* did not finish in 120 s.
- Not tested: tidk `plot`, other tidk versions, a large genome set.

## Conclusion
tidk does not replace the existing scripts. `denovo_telomere_repeats.py` already covers
tidk's `explore` role and does better. Tidk `search` is a useful independent cross-check.
Cheapest recall gain: raise `telomere_terminal_tolerance` (10 -> ~30) and test on a larger
set before changing the default (T-019).

## Addendum: proxy ground truth and terminal-tolerance test (2026-11-18)
Genomes: 40 random species (seed 1) with 5-40 contigs, N50 > 1.5 Mb, zero gaps
(`outputs/ground_truth_genomes.tsv`). 38 yielded a repeat unit, 2 did not.
Proxy truth (`score_ground_truth.py`): an end of a >=100 kb scaffold is a telomere when its
terminal 300 bp holds >= 8 copies (any rotation, either strand) of the genome's unit.
**Caveat:** the unit comes from `denovo_telomere_repeats.py`, so denovo is favoured in this
score. The truth is a proxy, not curated, and counts one unit per genome only.
Per-genome numbers: `outputs/ground_truth_scores.tsv`. Exact mode, `--both-ends`.

| Subset | Truth ends | find tol=10 | find tol=30 | denovo |
|---|---|---|---|---|
| All 38 genomes | 453 | 255 (56%) | 296 (65%) | 419 (93%) |
| Unit in default list (24 genomes) | 341 | 250 (73%) | 284 (83%) | 312 (91%) |
| Other units (14 genomes) | 112 | 5 (4%) | 12 (11%) | 107 (96%) |

- `--terminal-tolerance 30` adds 41 ends over tol=10 with 0 new false positives
  (2 in both runs). It is the cheapest gain.
- Genomes whose unit is not in `telomere_patterns` (long units such as the *Candida*/
  Saccharomycotina ones) are almost invisible to `find_telomeres.py`; denovo finds them.
- denovo's 21 calls outside the proxy truth mostly name a second unit at that end
  (e.g. `AGAGGGGGTC`), so they are not shown to be errors.
- Not investigated: why 57 of 341 default-unit ends are still missed at tol=30. Fuzzy mode was
  not scored (too slow on some genomes). tidk was not scored on this set.

## Addendum 2: the remaining misses, `--merge-gap`, and the pipeline change (2026-11-18)
**Why 57 default-unit ends were missed at tol=30.** 44 of the 57 sit within 30 bp of the edge.
The tracts hold variant repeats (`CCTAA`, `GGGG`) and 1 bp indels between exact copies, which break
an exact run below `--min-repeats`/`--min-length`. Tolerance was not the cause.

**Fix.** New `find_telomeres.py --merge-gap N` (exact mode, default 0): copies of a monomer
separated by up to N bp of any sequence form one tract. It is regex-based, so it is not frame-locked
like `--fuzzy`. Serial runtime on one genome: 0.93 s (N=0) vs 0.86 s (N=8).

Same 38 genomes and proxy truth, `--terminal-tolerance 30` throughout
(`outputs/ground_truth_scores_merge_gap.tsv`):

| Subset (truth ends) | tol30 only | gap 8 | gap 15 | gap 25 | denovo |
|---|---|---|---|---|---|
| Unit in default list (341) | 284 (2 outside truth) | 319 (4) | 319 (5) | 320 (5) | 312 |
| Other units (112) | 12 (0) | 45 (10) | 81 (12) | 91 (12) | 107 |
| All (453) | 296 | 364 | 400 | 411 | 419 |

- Gap 8 adds 35 default-unit ends (83% -> 94%) for 2 more calls outside truth.
- Gaps 15/25 add little for default units. Their extra "other unit" gains come with 12 calls
  outside truth that I did not inspect, so I did not adopt them.
- Still missed with gap 8: 22 of 341 default-unit ends. Not examined.
- Units not in `telomere_patterns` stay mostly invisible to `find_telomeres.py` (45/112). denovo remains the
  tool for those.

**Changed:** `telomere_terminal_tolerance` 10 -> 30, new `telomere_merge_gap = 8`
(`nextflow.config`, `conf/profile_BFD.config`, `nextflow_schema.json`, `BATCH_TELOMERES/main.nf`);
`find_telomeres.py` default tolerance 10 -> 30. Two tests added to `tests/test_find_telomeres.py`
(16 telomere tests pass). The denovo step is unchanged.

## How to re-run telomere identification
`BATCH_TELOMERES` uses `storeDir` = `results/genome_stats/telomeres_v2`, and `planFreq` skips
any genome that already has output there. A parameter change alone re-runs nothing.
1. Set the old results aside (same pattern as `telomeres.pre_both_ends_fix_20261112`):
   `mv results/genome_stats/telomeres_v2 results/genome_stats/telomeres_v2.pre_tol30_gap8_20261118`
   Do this only when no pipeline run is active. The path is 23,582 genomes in 256 buckets.
2. Keep the old merged tables until the new run succeeds:
   `cp tables/telomere_summary.parquet tables/telomere_summary.parquet.pre_tol30_gap8`
   (same for `telomere_tracts.parquet`).
3. Run only the BFD phase: `sbatch nextflow/run_unified.sh --only_phase1`.
   Other BFD steps (asm_stats, BUSCO, aa_freq, ...) have their own storeDirs and should not re-run.
4. `MERGE_TELOMERES` runs `check_table_coverage.py` against the old table (`coverage_min_frac` 0.9).
   If the run is partial, the merge fails instead of overwriting with a short table. Do not set
   `allow_shrink` to get past it. Re-run so the missing genomes are planned and finished.
5. Check: `telomere_summary.parquet` row count should be near the old count, and
   the share of genomes with a terminal call should rise. Stats from before the change are not
   comparable to stats after it.
6. After a good merge, delete the `.pre_tol30_gap8` directory to free `/bigdata`.

Task sizing: 23,582 genomes at `telomere_batch_size = 100` is about 236 tasks. I measured about
1 s (find) + 1 s (denovo) per genome on 30-40 Mb genomes, so a task may finish in minutes. That is
unconfirmed for the full set. If real task times are short, raise `telomere_batch_size` (the
global guidance is 1-1.5 h per task).
