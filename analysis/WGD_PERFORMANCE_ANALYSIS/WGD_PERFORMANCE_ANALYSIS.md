# WGD Per-Genome Performance Estimate & Output/Visualization Plan

## Purpose

Give a defensible, measurement-based estimate of how long the paralogoscope
(`--pipeline paralogoscope`) full run over the **4,365-genome**
`function_neurospora` input set will take, and lay out how the results should
be stored and visualized. The estimate is anchored on real wall times from the
r3/resume4 validation run (2 Aaosphaeria genomes), then extrapolated with
explicitly documented assumptions. Numbers are reproducible via the scripts in
this folder.

## 0. Workflow overview (what runs, in order)

Paralogoscope is a **per-species** dating pipeline: every genome flows through
`wgd dmd` then `wgd ksd` independently; `wgd syn` is optional (gated behind
`--run_wgd_syn`, **default `false`** — skip it in most cases; it needs the
GFF3 input and i-ADHoRe is heavy). No cross-genome/downstream-merging steps
exist — the workflow ends after the per-genome tasks publish.

```
 samples.csv ──┬── --group / --taxon / --ignore / --n_test ──► per-species selection
               │
               ▼
      input/cds/{id}.cds-transcripts.fa ──► (input/gff3/{id}.gff3 only if --run_wgd_syn)
               │
               ▼
   ┌───────────────────────┐
   │  STEP 1 — WGD_DMD     │  wgd dmd   (all-vs-all DIAMOND BLASTp + MCL)
   │  8 cpus / 32 GB / 24 h│  ~1 min per genome             (~61-64 s measured)
   │  SIF wgd-2.0.38       │
   └──────────┬────────────┘
              │
              │  {id}.cds-transcripts.fa.tsv   (MCL paralog families)
              │  └─► storeDir bucket: paralogoscope/wgd_dmd/{bucket256}/
              ▼
   ┌───────────────────────┐
   │  STEP 2 — WGD_KSD     │  wgd ksd    (per-family mafft → Ka/Ks dating)
   │  8 cpus / 32 GB / 24 h│  ~41 min per median genome @ 8 cpus (94.3 min @ 4)
   │  --nthreads = cpus    │  [the bottleneck step]
   └──────────┬────────────┘
              │
              │  {id}.cds-transcripts.fa.tsv.ks.tsv  (+ .ksd.pdf, .ksd.svg)
              │  └─► storeDir bucket: paralogoscope/wgd_ksd/{bucket256}/
              ▼
   ┌───────────────────────┐   (only if --run_wgd_syn true)
   │  STEP 3 — WGD_SYN     │  wgd syn     (i-ADHoRe multiplicon/synteny anchors)
   │  (opt-in, default off)│  NOT measured at scale — benchmark before enabling
   └──────────┬────────────┘
              │
              │  anchors.csv + *.dot.pdf
              │  └─► storeDir bucket: paralogoscope/wgd_syn/{bucket256}/
              ▼
         done (no merge step: results stay per-genome in the buckets)
```

Notes on the machinery, so the artifact/throughput review below is grounded in
what actually executes:

- **Execution model**: each step is a SLURM task on the `preempt` partition
  (`-N 1 -n 1`, `-A preempt -p preempt`), 8 cpus / 32 GB / 24 h, running the
  wgd SIF via inline `module load apptainer squashfuse fuse`. All analysis is
  local to `$SCRATCH`.
- **Container (no per-run pull)**: the profile default is a *materialized SIF
  in the shared apptainer cache* —
  `params.singularity_cache`/`hyphaltip_wgd2_complete-2.0.38--9e3a9744900b_0.sif`
  (sha256 `9e3a9744900b…` of the file; distinct from the biocontainers
  `wgd_2.0.38--pyhdfd78af_0.sif` in the same cache). `params.singularity_cache`
  resolves strictly from `NXF_APPTAINER_CACHEDIR` / `NXF_SINGULARITY_CACHEDIR`
  / `APPTAINER_CACHEDIR` / `SINGULARITY_CACHEDIR`, which `run_paralogoscope.sh`
  exports with the repo-standard backwards-compat default
  `/bigdata/stajichlab/shared/singularity_cache`. The image never pulls over
  the network in production; override with `--wgd_sif /path/to/other.sif`.
- **Pipeline-level concurrency**: dmd/ksd/syn for different genomes are fully
  independent, so a run is just a stream of per-genome tasks (ksd dominates
  the queue). There is no merge/aggregate process; the storeDir hash buckets
  *are* the final tree.
- **`storeDir` bucketing**: every step publishes through a SHA-1 hash bucket
  (`hashBucketForType(type, LOCUSTAG)`), so outputs live at
  `paralogoscope/wgd_{dmd,ksd,syn}/{bucket}/{file}` — no per-genome directory
  exceeds a few hundred files, and `-resume` reruns pick up completed genomes
  for free.

## 0a. Data artifacts produced

Each step's own `mv` moves its outputs to the workdir root, from where
`storeDir` copies them into the bucket. Filenames embed the source CDS id, so
a re-annotated sample yields a new filename and can never collide with a stale
storeDir hit.

### Step 1 — `wgd dmd` → families TSV (input to ksd/syn)

| artifact | content | size (real measurement) |
|----------|---------|------------------------|
| `{id}.cds-transcripts.fa.tsv` — "families" table | wgd MCL gene-family output. Header: `\t{cds filename}`; each row: `GF********\t{gene1},\ {gene2}, ...` (comma-space separated members, one family per row) | **0.33 MB** (A. pasadenensis: 327,528 B for 12,244 genes / 8,259 families) ⇒ ~0.24 MB at the 9,105-gene median |

No separate intermediate files are retained — only this one TSV escapes the
workdir.

### Step 2 — `wgd ksd` → per-pair Ks table + histogram plots

| artifact | content | size (real measurement) |
|----------|---------|------------------------|
| `{id}.cds-transcripts.fa.tsv.ks.tsv` — **the main data table** | one row per paralog gene pair with the Ka/Ks machinery columns: `pair, N, S, alignmentcoverage, alignmentidentity, alignmentlength, dN, dN/dS, dS, family, g1, g2, l, node, node_averaged_dS_outlierexcluded, node_averaged_dS_outlierincluded, strippedalignmentlength, t, weightoutlierexcluded, weightoutlierincluded, gene1, gene2`. dS (and the two node-averaged dS columns) are the WGD-dating inputs; `dN/dS` available per pair | **10.9 MB** (34,164 pair rows for A. pasadenensis) ⇒ ~8 MB at the median genome |
| `{id}.cds-transcripts.fa.tsv.ksd.pdf` | per-genome Ks histogram (density, colored/decorated wgd default plot) | **23 KB** |
| `{id}.cds-transcripts.fa.tsv.ksd.svg` | vector version of the same Ks histogram | **96 KB** |

**2026-08-28 review fix**: `WGD_KSD` previously moved/emitted only the
`.ks.tsv` + `.ksd.pdf`, silently dropping the `.svg` that wgd actually writes;
module now also captures `*.ksd.svg` (emit `ks_plot_svg`).

### Step 3 — `wgd syn` (opt-in) → synteny anchors

| artifact | content | size |
|----------|---------|------|
| `anchors.csv` | i-ADHoRe anchor points (gene-pair anchors across collinear/multiplicon regions) | *unmeasured* (see caveats: "No anchors found" is a legitimate rc=0 with no file) |
| `*.dot.pdf` | dot-plot images of collinear regions | *unmeasured* |
| (`iadhore-out/anchorpoints.txt` intermediate) | raw i-ADHoRe output retained in workdir | — |

### Full-generation size estimate (4,365 genomes)

| artifact | per genome (median) | × 4,365 | note |
|----------|---------------------|---------|------|
| dmd families TSV | ~0.24–0.33 MB | **~1.0–1.5 GB** | wgd_dmd buckets |
| ks.tsv | ~5–11 MB (genome-size dependent) | **~22–48 GB** | wgd_ksd buckets; dominates storage |
| ksd.pdf + ksd.svg | ~0.12 MB | ~0.5 GB | wgd_ksd buckets |
| syn anchors (+ dot pdfs) | TBD | TBD | opt-in |
| **raw total** | ~11–12 MB | **~35–50 GB** | comfortable on /bigdata |

Per-pair `ks.tsv` scales sub-linearly with genome size (pairs come from
multi-copy families only); the p90 (~13.8 k-gene) genome lands closer to the
11 MB end, the p10 (~5.2 k) to ~4–5 MB.

## 1. Empirical measurements (r3/resume4 run, local executor, login node, 4 cpus/task)

Source of truth: `collect_measurements.py` → `outputs/wgd_perf_measurements.tsv`
parsed from `/scratch/jstajich/27843486/plg_real/` (workdir `.command.begin`
mtimes, `.command.log` completion times, `.exitcode`, published dmd TSVs).

| step | genome | genes | multi-copy families | outcome | wall time |
|------|--------|-------|---------------------|---------|-----------|
| `wgd dmd` | Aaosphaeria pasadenensis | 12,244 | – | exit 0 | **~61 s** |
| `wgd dmd` | Aaosphaeria arxii | 13,393 | – | exit 0 | **~64 s** |
| `wgd ksd` | Aaosphaeria pasadenensis | 12,244 | 1,184 | exit 0 | **94.3 min** |
| `wgd ksd` | Aaosphaeria arxii | 13,393 | 1,262 | exit 143 (2 h `time` limit) | >120 min, incomplete |
| `wgd syn` | — | — | — | unmeasured at real scale | *TBD* |

Key observations:

- **dmd is cheap**: ~1 min/genome at 4 cpus for a 12–13 k-gene genome
  (DIAMOND all-vs-all + MCL). Never the bottleneck.
- **ksd is the bottleneck**: one genome (pasadenensis) fully completed:
  1,184 multi-copy paralog families → 34,164 gene pairs → `ks.tsv` + plots in
  94.3 min at 4 cpus.
- **ksd cost driver** = number/size of *multi-copy* families, not total genes.
  Singleton families (7,075 of 8,259 for pasadenensis) are skipped.
- **The two genomes have near-identical family-size spectra** (top families
  117/99/53 vs 123/101/56 members, 1,184 vs 1,262 multi-copy families), yet the
  arxii run appeared ~10x slower per family (151 "Analysing family" log lines
  in 2 h vs 1,184 in 94 min). The family compositions can't explain this — the
  most plausible causes are **login-node contention** (afternoon vs morning)
  and/or **unflushed python log buffering** in the killed task. Since the
  arxii log only shows 151 lines, the real progress may be undercounted, but
  the run still failed to finish within 2 h. This ~10x spread is treated as
  worst-case contention noise, **not** a biological signal, and argues for a
  SLURM pilot wave before locking estimates (see §4).

Raw size data per genome (drives storage estimate, §5):

- dmd families TSV: **~0.33 MB** (8,259 families / 12,244 genes)
- ks.tsv (34,164 pairs): **10.9 MB**; ksd.svg 96 KB + ksd.pdf 23 KB
- pasadenensis total published wgd outputs: ~11.4 MB

## 2. Extrapolation model

Normalization: ksd wall ∝ multi-copy-family load ≈ ∝ gene count. Measured
genome: 12,244 genes. Input-set median: **9,105 genes** (n=198 random sample;
mean 9,321, IQR 6,694–11,253, range 996–28,891 — wide spread, tail matters).

Assumptions (all stated, tunable in `estimate_throughput.py`):

1. 8 cpus ≈ **1.7x** the 4-cpu ksd throughput (`wgd ksd --nthreads ${task.cpus}`
   parallelizes per family across the threads; per-family mafft/codeml is
   single-threaded so scaling is sub-linear). **Unverified — pilot must check.**
2. Median-genome ksd time = pasadenensis 94.3 min × (9,105/12,244) / 1.7
   ≈ **41 min** at 8 cpus. Per-genome time varies roughly linearly with genes:
   p10 (~5.2 k) ≈ 24 min, p90 (~13.8 k) ≈ 62 min, the 28.9 k max ≈ 2.2 h.
3. dmd scales ~linearly with genes: ≈ 0.4 min at 8 cpus per median genome.
   dmd overlaps ksd of other genomes in the queue → adds ~nothing to wall.
4. `wgd syn` (i-ADHoRe) is **not yet measured at scale**; excluded from totals.
   Must be benchmarked on ≥2 real genomes before deciding to enable it
   (`--run_wgd_syn true`), since the SLURM profile gives each task 24 h.
5. Run overhead (resume, retries, scheduler idle) ≈ +10–15 % (not folded into
   the table below; the numbers are *compute* wall).

## 3. Full-run estimate (4,365 genomes)

Produced by `estimate_throughput.py` → `outputs/wgd_throughput_estimate.tsv`:

| step | cpus/task | per-genome (median) | total task-h | total CPU-h |
|------|-----------|---------------------|--------------|-------------|
| wgd dmd | 8 | ~0.4 min | ~30 | ~250 |
| wgd ksd | 8 | **~41 min** | **~3,000** | **~24,000** |
| wgd syn | 8 | TBD (opt-in) | TBD | TBD |

ksd wall clock by concurrency on the preempt partition (currently 7 nodes × 64
cpus "up" = 448 cpus, partly allocated; 8-cpu tasks):

| concurrent ksd tasks | 16 | 24 | 32 | 40 | 56 |
|----------------------|----|----|----|----|----|
| wall (days) | 7.8 | 5.2 | 3.9 | 3.1 | 2.2 |

**Bottom line: ~3,000 task-hours / ~24,000 CPU-h for ksd, i.e. roughly 2–8 days
of wall clock at 16–56 concurrent tasks.** At the high end (28.9 k-gene
genomes, or if 8-cpu scaling is worse than 1.7x) expect the tail to stretch:
plan for up to **~10 days** wall in the worst realistic case; a quick pilot
wave narrows this.

## 4. Recommended run plan

1. **Pilot wave (next step):** submit ~20 genomes spanning the size spectrum
   (p10→p90) via a small `samples.csv`/`--n_test` + `--ignore` slice on SLURM
   (`sbatch nextflow/run_paralogoscope.sh ...`), 8 cpus, preempt queue.
   - Measure real per-task wall from `logs/nextflow/paralogoscope_trace.*.txt`.
   - Verify the 8-vs-4-cpus scaling factor and the per-genome time spread.
   - Confirm ksd memory stays well under 32 GB.
2. **syn benchmark (parallel):** run `wgd syn` on the 2 pilot genomes
   (`--run_wgd_syn true`) to measure anchors.csv time/memory before committing.
3. **Full run:** launch in waves with `-resume`; the hash-bucket `storeDir`
   makes reruns pick up completed genomes for free. Waves also protect against
   preempt-partition volatility. Monitor via the trace file.
4. **Local-test caveat:** the test config's `time='2 h'` killed the arxii ksd;
   if ever re-run locally raise the limit, but SLURM (24 h) makes this moot.

## 5. Storage plan

Raw wgd outputs (per genome ≈ 11–12 MB) already land in BFD hash buckets:
`paralogoscope/wgd_{dmd,ksd,syn}/{bucket}/`, `storeDir` = cache + publish.

| artifact | per genome | × 4,365 | notes |
|----------|-----------|---------|-------|
| dmd families TSV | ~0.33 MB | ~1.5 GB | wgd_dmd buckets |
| ks.tsv (gene pairs) | ~7–11 MB | **~30–45 GB** | wgd_ksd buckets; largest |
| ks plots (svg+pdf) | ~0.12 MB | ~0.5 GB | wgd_ksd buckets |
| syn anchors.csv | TBD | TBD | opt-in |
| **raw total** | ~11–12 MB | **~35–50 GB** | comfortable on /bigdata |

Derived tables (recommended, small, become the analysis hub):

- **`tables/wgd.ks.parquet`** — **merged analysis hub** (final step of the
  pipeline, `MERGE_WGD_KSD` → `bin/merge_wgd_ks.py`): one row per
  (genome, wgd gene pair), pruned to the 14 analysis columns (string
  `pair/family/g1/g2/gene1/gene2/node` + float `N/S/dN/dN/dS/dS/
  alignmentlength/t`; `node` is the wgd tree-node id the pair was dated at
  and is required by wgd's own mix/peak node-averaging) and keyed by
  `species_prefix` (LOCUSTAG) + `genome` (sampletag). zstd; ≈ 19 MB per 1,000
  genomes of raw ks.tsv → **~0.8–1 GB for the full run** (vs ~45 GB raw).
  ~43 % of rows are NULL-dS placeholders (wgd emits one row per gene pair
  regardless of whether Ka/Ks dating succeeded; cheap in parquet).
- **`tables/wgd.ks.summary.parquet`** — per-genome counts: n_pairs,
  n_pairs_with_ds, n_families (4,365 rows).
- **`tables/wgd_ksd_summary.parquet`** — **implemented** (see "Derived-table
  builder" below): one row per genome (4,365): species_prefix, the summary
  counts, `n_ks` (pairs passing the wgd fit filter), `n_ks_peaks` (best-model
  component count) + `bic_best`/`aic_best`, and wide per-peak columns
  `peak1_ks_mean/sd/weight … peak4_...` (peaks ordered by increasing mean
  Ks; NaNs past the genome's component count). Joins samples.csv/BFD tables.
- **`tables/wgd_ksd_density.parquet`** — long form: genome × shared dS grid
  (200 bins on dS∈[0,5]) of the per-genome KDE (silverman) → ~0.9 M rows, tens
  of MB. This is the artifact the aggregate figures consume.
- Raw ks.tsv (35–50 GB) stays in the wgd_ksd buckets as the archive; per-clade
  subset extraction reads directly from `tables/wgd.ks.parquet` (column
  pruning + predicate pushdown make per-genome reads sub-second).

**Derived-table builder** — answers the per-genome question "how many
duplicates, and what is the mean Ks of each Ks peak?" using **wgd's own
mixture-model machinery** (there is no separate wgd script for this: `wgd
mix`/`wgd peak` only ship per-pair posterior TSVs + plots, never the
component means):

- `scripts/wgd_ks_mix_fit.py` — worker that runs **inside the wgd SIF** and
  calls `wgd.mix.filter_group_data` + `wgd.mix.fit_gmm`/`fit_bgmm` (the exact
  code behind `wgd mix`; node-averaged pairs, min-BIC model selection,
  seed 2352890 = wgd peak's default). Emits a JSON with per-model
  BIC/AIC and per-component mean-Ks/sd/weight plus a silverman KDE on
  dS∈[0,5].
- `scripts/build_wgd_ksd_summary.py` — host driver: slices
  `tables/wgd.ks.parquet` per genome (`pair/family/node/alignmentlength/dS`),
  runs the worker in parallel (`--jobs`), and assembles
  `wgd_ksd_summary.parquet` + `wgd_ksd_density.parquet` (zstd). Resolves
  apptainer from PATH or the HPCC module install; ~2.5 CPU-min/genome at
  n_init=200 ≈ 180 core-h / 4,365 genomes (~3 h at 64-way).

```bash
python3 analysis/WGD_PERFORMANCE_ANALYSIS/scripts/build_wgd_ksd_summary.py \
  --sif /bigdata/stajichlab/shared/singularity_cache/hyphaltip_wgd2_complete-2.0.38--9e3a9744900b_0.sif \
  --ks-parquet tables/wgd.ks.parquet --summary-parquet tables/wgd.ks.summary.parquet \
  --outdir analysis/WGD_PERFORMANCE_ANALYSIS/tables --jobs 8
```

**Validated on the 2-genome synb pair** (fits match a direct in-container
`fit_gmm` run): arxii `n_ks=893`, **3 peaks** at mean Ks 0.080 / 1.559 / 3.500
(weights 0.080 / 0.191 / 0.729, BIC 1388.57); pasadenensis `n_ks=813`,
**4 peaks** at 0.053 / 0.853 / 2.449 / 3.766 (weights 0.019 / 0.129 / 0.300 /
0.552, BIC 1052.96).

## 6. Visualization plan

Per-genome (wgd already emits): `ksd.pdf`/`ksd.svg` per genome ⇐ gallery index
(HTML listing 4,365 genomes with embedded per-genome Ks plots).

Aggregate (derive from `wgd_ksd_density.parquet` + `wgd_ksd_summary.parquet`):

1. **Ks peak calling** — **implemented** as model-based peaks via the wgd GMM
   (`wgd mix` machinery, min-BIC, components = peaks; see "Derived-table
   builder" in §5): `n_ks_peaks` + per-peak mean/sd/weight land directly in
   `wgd_ksd_summary.parquet`. Recent WGD events show a substantial low-Ks
   peak (≈0.1–0.5); slowly-decaying single-copy genomes show fewer/no peaks
   → classify each genome **WGD-candidate / ambig / none** from
   n_ks_peaks + peak1 position/weight. (KDE `find_peaks` kept as an optional
   cross-check.)
2. **Peak-position distribution** — histogram/rug of 1st-peak Ks, faceted by
   taxonomic class (`samples.csv` rank columns) → the WGD landscape/recent
   duplication pulse.
3. **Class-overlaid Ks density** — mean ± CI per class on the shared grid
   (small multiples per class).
4. **Profile-space embedding** — each genome's KDE vector → PCA/UMAP, colored
   by class, marker by WGD-call; exposes clades behaving differently.
5. **Interactive report** — plotly HTML combining 1–4 with click-through from
   genome points to the stored per-genome `ksd.pdf`.
6. **Run/throughput dashboard** (ops, not science): wall-time CDF from the
   trace once the full run starts, updating the per-genome model with real data.

## 7. Open questions / uncertainties

- 8-cpu ksd scaling (assumed 1.7x) — verify on pilot.
- Contention-driven spread (the arxii 10x) — quantify on dedicated nodes.
- Gene-count distribution for the *full* 4,365 (n=198 sample) — refinable from
  existing genome-stats tables if present.
- syn anchors volume and time — **resolved, see §8a** (dmd-tier, not a bottleneck).
- Peak-calling thresholds (Ks window, KDE bandwidth) — science decisions to
  make once the first real aggregate plots exist.

## 8. Pilot wave measurements (2026-08-28, SLURM preempt, 8 cpus)

Pilot: 20 genomes spanning the rank-sorted gene-count spectrum of the 4,365
input (354 → 18,307 genes; `samples.csv` slice in
`outputs/pilot_samples_20.csv`, selection table in
`outputs/pilot_genome_counts.tsv`). Submitted 18:32 via
`run_paralogoscope.sh` (job 27927803), `-process.queueSize 20`, image
`hyphaltip_wgd2_complete-2.0.38--9e3a9744900b_0.sif`; **run_wgd_syn = false**
(opt-in; syn still unbenchmarked at scale). All 40 tasks completed exit 0 in
**~1 h 50 m wall**. Trace: `logs/nextflow/paralogoscope_trace.2026-08-28_18_32_15.txt`;
per-genome table: `outputs/pilot_profile.tsv` (`scripts/profile_pilot.py`).

**Per-stage timing.** `wgd dmd` is trivially fast (11–15 s/genome). `wgd ksd`
dominates: 8 min (Metschnikowia_rubicola, 6,235 genes) to 111 min
(Talaromyces_liani_FKII-L2-CM-DRAB3, 11,868 genes) and 83 min
(Paramarasmius_palmivorus, 18,307 genes). Within-size run-to-run variance is
large (e.g. Schizophyllum 27 min vs Talaromyces_liani 111 min at ~11.9–11.8k
genes) — ksd cost tracks protein-family structure (multi-copy family count),
not raw gene count.

**Rates per 5,000 genes** (8 cpus, preempt):
| metric | linear slope (per 5k) | R² | median observed rate | interpretation |
|---|---|---|---|---|
| total runtime | **24.8 min / 5k** | 0.37 | 20.0 min / 5k (mean 37.3, small-genome outliers) | weak predictor; expect per median-genome (9,578 genes) ≈ **38–45 min** |
| peak RSS / task | −0.17 GB / 5k (≈flat) | 0.03 | 1.28 GB / 5k | **memory does not scale with genome size** — 1.8–3.3 GB/task, median ~2.5 GB |
| artifact size | **12.7 MB / 5k** | 0.69 | 10.6 MB / 5k | ks.tsv dominates; Polyploid 18k-gene case 61.5 MB |

**Extrapolation to the full 4,365 genomes** (41,155,977 total genes, median
9,578/genome):
- **Total compute ≈ 3,277 task-hours** → at concurrency 16: ~205 h; **32: ~102
  h (4.3 d)**; 56: ~58 h (2.4 d). Tasks are 0.1–2 h each, well under the
  preempt 24 h cap.
- **Total artifacts ≈ 92 GB** (avg ~21 MB/genome; ks.tsv dominates) — trivial
  for /bigdata, plus syn ≈350 KB/genome; budget the 52 dangling-symlink genomes
  (zero-gene) as no-op skips (or drop their links, §8).
- **Memory is not a constraint**: 32 GB/task cap is >10× the ~2.5 GB median
  peak. Keep the cap for the heavy multi-copy tail (cf. pasadenensis 11k-gene
  spike observed in §1).

**Refinement path:** regress ksd runtime/RSS against multi-copy family count
(available from each `wgd_ksd` input families TSV) instead of raw gene count
to cut the R²=0.37 noise; re-fit after the first full-run wave.

**Data-quality flag:** 52 `input/cds/*.cds-transcripts.fa` symlinks dangle
(46 at pilot time — the tree is still being rewritten by the upstream
annotation wave, e.g. Ogataea ran fine in the pilot then its dir vanished):
38 have an empty `predict_results/`, 12 have no `genome_annotation/<tag>/` dir
at all, 2 have no `predict_results/` subdir. None have a valid CDS anywhere
else (not re-pointable). 51/52 have rows in `samples.csv`, so a full-samples
run would read-fail them; exclude via an ignore list or drop the dangling
links at launch time. All 4,365 gene counts are unaffected (grep-based).

## 8a. wgd syn benchmark (synb, 2026-08-28, SLURM preempt, 8 cpus)

Two-genome validation of `wgd syn` (opt-in, previously unbenchmarked): the r3
genomes Aaosphaeria_arxii_CBS_175.79 (13,393 genes) and
Aaosphaeria_pasadenensis_FJI-L9-BK-P1 (12,244 genes). Submitted 20:52 via
`run_paralogoscope.sh` (job 27929339), `-process.queueSize 2`, `--run_wgd_syn
true`, outdir `paralogoscope_synb` (separate storeDir from the main run).
COMPLETED in 58:28. Trace: `logs/nextflow/paralogoscope_trace.2026-08-28_20_57_00.txt`.

| stage | arxii (13,393 g) | pasadenensis (12,244 g) |
|---|---|---|
| dmd | 39.7 s | 36.3 s |
| **syn** | **36.9 s** | **37.2 s** |
| ksd | 55 m 23 s | 42 m 31 s |
| peak RSS / task | 2.5 GB | 2.6 GB |

**Result: syn is dmd-tier, not ksd-tier — ~37 s/genome (~0.5 GB peak RSS).
It is not a bottleneck.** At 4,365 genomes syn adds ≈45 task-h (<1.5% of the
~3,300 task-h total); the full-run estimate stands unchanged. `wgd_syn`
publishes only `anchors.csv` + a `*-vs-*.dot.pdf` per genome (~175 KB each,
≈350 KB/genome); both genomes produced real anchors, so the zero-anchor guard
was not exercised here (no fix needed on these data). Ksd re-measured at
8 cpus (42.5 / 55.4 min) is consistent with the §8 pilot spread, confirming
ksd remains the sole cost driver.

## Files

- `scripts/collect_measurements.py` — re-derives the r3 measurement table.
- `scripts/estimate_throughput.py` — the extrapolation model.
- `scripts/profile_pilot.py` — pilot trace → per-genome runtime/RSS/artifact
  table + linear-fit extrapolation (per 5k genes, full-dataset projection).
- `outputs/wgd_perf_measurements.tsv` — measured step timings.
- `outputs/pilot_samples_20.csv` — 20-genome pilot samples.csv slice.
- `outputs/pilot_genome_counts.tsv` — rank / sampletag / gene count for the
  pilot slice.
- `outputs/pilot_profile.tsv` — per-genome dmd/ksd min, peak RSS, artifact MB.
- `outputs/wgd_throughput_estimate.tsv` — full-run estimate table.
- `nextflow/bin/merge_wgd_ks.py` — final merge of per-genome ks.tsv →
  `tables/wgd.ks.parquet` + `tables/wgd.ks.summary.parquet` (globs disk,
  prunes to analysis columns, streams file-by-file into one zstd parquet).
- `nextflow/modules/comparative/wgd/MERGE_WGD_KSD/main.nf` — terminal pipeline
  step (label 'merge', publishDir tablesDir()); validated end-to-end on the
  synb data (72,755 rows, 2 genomes) including a fully-cached `-resume`
  re-trigger of the merge.
