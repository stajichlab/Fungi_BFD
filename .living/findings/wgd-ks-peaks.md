# wgd Ks-peak structure of the first two genomes

## Observation
Fit per-genome Gaussian-mixture models to the paralog-pair Ks distribution
using **wgd's own machinery** (`wgd.mix.filter_group_data` + `wgd.mix.fit_gmm`,
the exact code behind `wgd mix`; seed 2352890, node-averaged pairs,
alignment length ≥ 300, dS ∈ (0, 5], min-BIC model selection over 1–4
components). On the two *Aaosphaeria* synb genomes (2026-08-30):

- **A. arxii** (38,591 pairs, 1,262 multi-copy families; `n_ks` = 893 after the
  wgd fit filter): **3 components** (BIC 1388.57) at mean Ks **0.080 / 1.559 /
  3.500**, weights **0.080 / 0.191 / 0.729**, sd 1.888 / 0.500 / 0.220.
- **A. pasadenensis** (34,164 pairs, 1,184 families; `n_ks` = 813): **4
  components** (BIC 1052.96) at mean Ks **0.053 / 0.853 / 2.449 / 3.766**,
  weights **0.019 / 0.129 / 0.300 / 0.552**.

Both genomes are dominated by the high-Ks (≈3.5–3.8, saturating) component
(55–73% of pairs); the low-Ks "recent duplication" mass is small (8% / 2% at
Ks < 0.1) — consistent with no recent whole-genome duplication pulse in either
genome, matching their single-copy-dominated gene content.

## Validation
- Framework reproduces a direct in-container `fit_gmm` run (relative error
  < 1% on component means; seed-fixed). Reproducible via
  `analysis/WGD_PERFORMANCE_ANALYSIS/scripts/build_wgd_ksd_summary.py`
  → `analysis/WGD_PERFORMANCE_ANALYSIS/tables/wgd_ksd_summary.parquet`.
- wgd's `mix`/`peak` CLIs were checked: they never export component means
  (only per-pair posterior TSVs + plots), so the worker captures them.

## Implications
- The low-Ks component (0.05–0.08, sd > 1.8 in these two) is the direct
  WGD/duplication-pulse signal; a per-genome heuristic (e.g. peak1 weight and
  n_ks_peaks) can classify genomes as **WGD-candidate / ambig / none** once
  the full 4,365-genome `wgd_ksd_summary.parquet` exists (T-032).
- `n_ks` (≈800–900 here vs 17–22 k pairs with dS) shows how aggressive the
  wgd default aln-length/dS filter is at the pair level — worth a documented
  sensitivity check if peak means are compared across runs.

## Status
preliminary (n = 2 genomes; full-dataset peak landscape pending the full
paralogoscope run, T-032).
