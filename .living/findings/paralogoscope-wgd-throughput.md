# Paralogoscope/wgd throughput characteristics

## Observation
The two-step `wgd` dating pipeline (committed 2026-08-28) has a very uneven
per-genome cost profile, measured on the r3/resume4 real-data validation run
(2 *Aaosphaeria* genomes, local executor, 4 cpus/task):

- **`wgd dmd`** (DIAMOND all-vs-all + MCL family clustering) is cheap:
  ~1 min/genome for a 12–13 k-gene genome (61–64 s measured). Never the
  bottleneck.
- **`wgd ksd`** (K_a/K_s dating of paralog pairs) is the bottleneck:
  94.3 min for *A. pasadenensis* (12,244 genes, 1,184 multi-copy families,
  34,164 gene pairs → 10.9 MB `ks.tsv` + svg/pdf plots). The runtime is driven
  by the *multi-copy* family load, i.e. ≈ gene count (singleton families are
  skipped — 7,075 of 8,259 families for pasadenensis).
- Runtime variance is large and **not biological**: *A. arxii*, with a
  near-identical family-size spectrum (1,262 vs 1,184 multi-copy families,
  top families 123/101/56 vs 117/99/53 members), was SIGTERM'd at the 2 h test
  `time` limit after only 151 "Analysing family" log lines. Early-afternoon
  login-node contention and/or unflushed python log buffering in the killed
  task are the plausible causes — treated as worst-case contention noise.

## Validation
- Timings from `/scratch/jstajich/27843486/plg_real/work_r3/*/*/.command.begin`
  mtimes + `.command.log` completion times + `.exitcode` (ksd pasadenensis
  ".command.log" prints `Total run time: 94.31 minutes`); reproduced by
  `analysis/WGD_PERFORMANCE_ANALYSIS/scripts/collect_measurements.py`.
- Per-family "Analysing family" line counts in ksd logs match the published
  multi-copy family counts exactly (1,184 / 1,262), confirming the skips.
- The 4,365-genome input set (`function_neurospora/input/cds`) has median
  **9,105 genes** (n=198 random sample; mean 9,321, IQR 6,694–11,253, range
  996–28,891).

## Implications
- Extrapolated full-run cost (ksd at 8 cpus, assumed 1.7x the 4-cpu rate;
  normalized to median gene count): **~41 min per median genome → ~3,000
  task-h / ~24,000 CPU-h → 2.2–7.8 days wall clock** at 16–56 concurrent
  tasks on the preempt partition. dmd adds ~30 task-h and overlaps the ksd
  queue, so it is wall-clock-neutral. `wgd syn` (i-ADHoRe) is **unmeasured at
  scale** and must be benchmarked before enabling.
- The contention-driven ~10x spread means the 41-min number carries real
  uncertainty: run a ~20-genome SLURM **pilot wave** spanning the size
  spectrum (p10→p90) to verify the 8-cpu scaling factor and per-genome spread
  *before* committing the full run.
- Storage is a non-issue at this scale: raw buckets ≈ 11–12 MB/genome
  (**~35–50 GB total**), plus small derived parquet summary/density tables
  that become the analysis hub — the 45 GB of raw `ks.tsv` is never merged
  into one table. See
  `analysis/WGD_PERFORMANCE_ANALYSIS/WGD_PERFORMANCE_ANALYSIS.md` for the
  full storage + visualization plan (per-genome ksd.pdf gallery, Ks peak
  calling for WGD dating, class-overlaid densities, profile embeddings).

## Pilot wave (2026-08-28) — 20 genomes, SLURM preempt, 8 cpus, measured

The pilot (job 27927803, 18:32→~20:25 wall, queueSize 20, all 40 tasks exit 0)
now pins the earlier estimates with real SLURM data; trace
`logs/nextflow/paralogoscope_trace.2026-08-28_18_32_15.txt`, table
`analysis/WGD_PERFORMANCE_ANALYSIS/outputs/pilot_profile.tsv`:

- **Per 5,000 genes** (8 cpus): runtime slope **24.8 min** (R²=0.37; median
  observed rate 20.0 min, mean 37.3 — small-genome outliers inflate the mean);
  artifact size slope **12.7 MB** (R²=0.69, ks.tsv-dominated); peak RSS is
  **flat** (≈2.5 GB/task median, R²≈0) — memory does not scale with genome
  size, the 8 cpus/32 GB/task cap is comfortable. dmd remains trivial
  (11–15 s/genome).
- **Full 4,365-genome extrapolation**: median 9,578 genes → **~3,277
  task-hours** → ~102 h (4.3 d) at 32 concurrent, ~58 h at 56; **~92 GB**
  artifacts; worst-case genome seen 111 min ksd (Talaromyces_liani, 11,868
  genes) — all < 24 h preempt cap.
- **Runtime variance is biological, not just contention**: at equal gene
  counts (~11.8–12.0 k), Schizophyllum_commune ran 27 min while
  Talaromyces_liani ran 111 min, both on dedicated preempt nodes — ksd cost
  tracks multi-copy family structure, so raw gene count is a weak predictor
  (R²=0.37). Refinement: regress on multi-copy family count from the dmd
  families TSV, and re-fit after the first full-run wave.
- **Data QC flag**: 52 of 4,365 `input/cds/*.cds-transcripts.fa` symlinks
  dangle (46 at pilot time — the annotation tree is still being rewritten; 38
  empty `predict_results/`, 12 no `genome_annotation/<tag>/` dir, 2 no
  `predict_results/`; none re-pointable, 51/52 in samples.csv), presenting as
  zero-gene entries; excluded from the pilot, must be handled before the full
  run.

## wgd syn benchmark (2026-08-28) — dmd-tier, not a bottleneck

`wgd syn` (i-ADHoRe), run on the same Aaosphaeria pair (12,244 / 13,393 genes)
via job 27929339 (run_wgd_syn true, outdir `paralogoscope_synb`, trace
`logs/nextflow/paralogoscope_trace.2026-08-28_20_57_00.txt`):

- **~37 s/genome** (36.9 / 37.2 s) at 8 cpus, ~0.5 GB peak RSS — the same
  order as `wgd dmd`, **not** the ksd-class cost previously feared. Anchors
  exist on both genomes (zero-anchor guard not triggered).
- **~350 KB/genome** artifacts (`anchors.csv` + `*-vs-*.dot.pdf` only).
- At 4,365 genomes syn adds ≈45 task-h (<1.5% of the ~3,300 task-h total) —
  the full-run estimate is unchanged and `--run_wgd_syn true` is safe to
  enable for the full run.
- Confirm on a couple of small/medium genomes in the first full-run wave if
  i-ADHoRe cost scales with genome size; here both genomes were the
  mid-large tier and stayed trivial.
