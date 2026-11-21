# Paralogoscope dN/dS estimator: yn00 (kaks) or codeml

Status: 2026-11-20. Branch `paralogoscope-giant-cap`.

## Decision
- Default: `wgd_dnds_method = yn00` (kaks, Yang and Nielsen 2000), per pair, on the gap-free codon alignment.
- Optional second pass: `wgd_refine_codeml_min` / `wgd_refine_codeml_max`. Pairs whose yn00 dS is inside the range are re-estimated with codeml on the same alignment. Default max = 0 (off). Suggested range: 0.5 to 2.
- In the refined rows, `dS`, `dN`, `dN/dS`, `t`, `S`, `N` hold the codeml values. `dS_yn00` holds the yn00 dS. `estimator` is `codeml`. If codeml fails, the yn00 values stay and `estimator` stays `yn00`.
- The merge step keeps its own column list, so `estimator` and `dS_yn00` stay in the per-genome files and are not in the merged table.
- Report dS above about 1.5 as rank order only.

## Evidence (all measured, one run each)
Speed (`eval/bench_excl`, isolated node): 2.8 ms per pair with yn00, 86 ms with codeml. Whole genome, 16 workers: Mycena 240 s against 1,382 s; Allomyces macrogynus 102 s against 382 s.

Agreement on 3,034 pairs: yn00 and codeml agree below dS about 0.5. Above that, the codeml/yn00 median ratio is 1.45 at yn00 dS 0.5 to 0.8 and about 2 above 0.8. yn00 saturates near 2.

Accuracy against known dS (`eval/kaks_browser_testset/manifest.tsv`, 384 simulated pairs, 50 per bin above 0.3, indel sets included; stratified, not a random sample). "Within" means |estimate - true| <= max(0.1, 0.5 x true).

| True dS | codeml median est/true | yn00 median est/true | codeml within | yn00 within |
|---|---|---|---|---|
| 0-0.3 (n=84) | 1.16 | 1.21 | 57% | 57% |
| 0.3-0.6 | 1.17 | 1.22 | 58% | 58% |
| 0.6-1.0 | 1.13 | 1.38 | 64% | 56% |
| 1.0-1.5 | 1.29 | 1.56 | 60% | 46% |
| 1.5-2.0 | 1.89 | 1.36 | 40% | 58% |
| 2-3 | 1.77 | 1.28 | 40% | 70% |
| 3-5 | 6.16 | 0.86 | 28% | 84% |

Reading: codeml is equal or better up to true dS 1.5. yn00 is better above 1.5. Neither is accurate above 2 in these sets. The refinement range of 0.5 to 2 follows the yn00-vs-codeml disagreement, not a measured accuracy gain, because codeml is not clearly better than yn00 in the simulation for true dS above 1.0. The refinement pass is therefore a consistency option with the earlier codeml results, not a proven accuracy fix.

## Not decided
Giant families. Five genomes were skipped in `nf_full/tables_yn00` (Ambispora leptoticha, Geosiphon pyriformis, Gigaspora margarita, Jimgerdemannia flammicorona, Jimgerdemannia lactiflua). Their largest families have up to 3,041 genes (4.6 million pairs). The protein-alignment first pass takes 99% of worker time (`eval/giant_scale`: 1,800 genes, 1.6 million pairs, 30 min on 16 CPUs). Capping family size loses close pairs (`eval/cap_recluster`). Options: process the first N genes, subsample pairs, or run these genomes alone with a long walltime.

## Test of the refinement pass (2026-11-20, `eval/refine_test`)
Allomyces macrogynus, node sampling, 16 CPUs, one run each. Without refinement: 132 s wall, 19,933 ok rows. With range 0.5 to 2: 219 s wall (1.7x), 20,513 ok rows, 8,183 rows refined with codeml (every refined row has yn00 dS from 0.5 to 2.0; none left unrefined inside the range). Median codeml/yn00 dS ratio in refined rows: 1.07. The row count differs because node sampling gates on the refined dS, so different pairs are drawn. `dS_yn00` is filled in all rows. Not tested: a full Nextflow run with the new parameters.

## Giant families: K random partners per gene (2026-11-20, `eval/giant_k200`, `giant_kvar`, `giant_scaled`)
Implemented as `--giant-genes N` and `--giant-k K` (params `wgd_giant_genes`, default 0 = off; `wgd_giant_k`, default 200). In a family with more than N distinct sequences the first pass aligns K random partners per gene, not every pair. The NJ tree needs a full matrix. Each unaligned distance is the shortest two-step path through common partners, multiplied by the median ratio (aligned distance / two-step path) measured on 20,000 aligned pairs (0.69 to 0.78 here). Node sampling then draws only among aligned pairs. Unaligned pairs are counted as status `giant_unaligned`. The family list is written to `<out>.giant.tsv.gz`.

Test: first 1,800 genes of the largest Jimgerdemannia lactiflua family (1,619,100 pairs), yn00, 16 CPUs, one run per setting. Reference is the full all-pairs run (1,808 s). Median and quantiles are of the node-mean dS (gated pairs, dS <= 2.5, >= 30 codons).

| Run | Wall | Nodes | 10% | 25% | 50% | 75% | 90% |
|---|---|---|---|---|---|---|---|
| full | 1,808 s | 1,797 | 0.157 | 0.294 | 0.488 | 0.762 | 1.055 |
| K=100, scaled fill | 136 s | 1,771 | 0.216 | 0.370 | 0.560 | 0.840 | 1.170 |
| K=200, scaled fill | 250 s | 1,794 | 0.192 | 0.347 | 0.553 | 0.841 | 1.190 |
| K=400, scaled fill | 449 s | 1,792 | 0.170 | 0.324 | 0.522 | 0.815 | 1.163 |
| K=200, no scaling (first version) | 312 s | 1,794 | 0.207 | 0.382 | 0.608 | 0.935 | 1.257 |

Result: the sampled runs are 7 to 13 times faster (K=200: 7x) but read high. At K=200 with scaling the median is 13% above the full run and the 10% quantile is 22% above. At K=400 the median is 7% above. Bias falls as K rises. The fill scaling removed about half of the unscaled bias. A 600-gene family at K=200 (unscaled) gave a median of 0.540 against 0.505 (+7%). I did not test the cause of the remaining bias (tree topology from the filled distances is my hypothesis). Not tested: other giant families, a full Nextflow run, the five skipped genomes.

Use: the feature is off by default. For the five skipped genomes, treat the dS distribution of their giant families as biased high by about 10 to 15% at K=200, and flag them in the results.
