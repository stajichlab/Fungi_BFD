# Paralogoscope pair-level dating: what it gives and what it can lose

Status: 2026-11-01. Written after the cap, re-cluster, pre-filter and Nextflow tests.
Evidence: `eval/sim_controls/report.md` (sections of 2026-11-01) and `eval/cap_recluster/`
in `paralogoscope_pilot/`. All numbers below are measured, one run each, unless the text says otherwise.

## 1. What the default method does now

`wgd_ks_method = pair` (module `WGD_PAIR_KS`, script `bin/wgd_pair_ks.py`).

1. `wgd dmd` makes the families (unchanged: e-value 1e-5, MCL inflation 1.5).
2. For every gene pair in a family, the script makes one global protein alignment (BLOSUM62). It does not build a family alignment.
3. It computes protein identity and NG86 dS from that alignment.
4. A pre-filter keeps pairs with identity >= 0.4 and >= 30 aligned codons.
5. One pairwise `codeml` runs on each kept pair.
6. The merge step keeps pairs with dS <= 2.5 and >= 30 codons used by `codeml`. It stores the number of codons as `pair_codons`.

Other methods remain available: `ksd_pairwise` and `ksd`. Each writes to its own folder.

## 2. What we gain (measured)

- Task time at 16 CPUs: 1 to 7 min for six of seven test genomes. Rhizopus took 19.5 min.
- Five genomes shared with the Oct 29 pilot: measured CPU use (realtime x %cpu from the Nextflow traces) is 5.7 CPU-hours now and 38.0 before, a 6.7x ratio. On allocated CPUs the figures are 6.8 and 71. The old run used mixed node types. Like-for-like offline runs (same genome, same pairs, same CPU model, 32 CPUs) show 5.7x (Candida), 8.5x (Zymoseptoria) and 2.1x (Rhizophagus families up to 50 genes) for all pairs. Most of the gain in small-family genomes comes from the pre-filter, which skips pairs. The Nextflow timings were measured on the `preempt` queue, the offline runs on `intel`.
- Default `wgd ksd` gives no usable dS for many close pairs, because it strips every alignment column that has a gap in any family member. In Rhizopus the old gate kept 13.1% of the NG86-close pairs. The pair-level gate keeps 99.9%.
- NG86 close-pair recall of the final table: 99.9 to 100% in all seven test genomes. This figure is mostly circular (see section 3.2) and is not independent evidence.

## 3. What can be lost or changed

### 3.1 Pairs that get no dS

- **Pre-filter.** Pairs below 40% identity or below 30 aligned codons are not sent to `codeml`. The share of pairs sent to `codeml` in the Nextflow run was 55% (Candida), 55% (Rhizopus), 38% (Zymoseptoria), 45% (Coprinopsis), 15% (N. crassa), 15% (Botrytis) and 12% (Hortaea). In an offline test on Rhizophagus families up to 50 genes the identity filter passed 63%.
- **dS > 2.5.** Pairs above the cut are dropped. Old duplications are not dated. NG86 saturates near dS 2. I did not test `codeml` accuracy in simulation for true dS above 2, so I cannot say how reliable dS 2 to 2.5 is.
- **Short alignments.** In simulation (true dS <= 2), accuracy rose with alignment length. With the lenient criterion |estimate - truth| <= max(0.1, 0.5 x truth), 56.6% of pairs are accurate at 20 codons, 65.0% at 30, 71.4% at 50 and 85.1% at 150. With the stricter criterion of within 1.5x of truth, the figures are 38.2%, 47.8%, 60.4% and 82.1%. At the 30-codon limit fewer than half of the pairs are within 1.5x of the true dS. There is no cliff. The 30-codon limit is a trade-off.
- **Genes with no family.** Genes that `wgd dmd` leaves as singletons are never examined. This is unchanged from the old method.

### 3.2 Pairs that differ from the old `--pairwise` result

Compared with `wgd ksd --pairwise` on the same families, at the 2.5 cut:

| Genome | Baseline gated pairs kept | Pairs only in the new result |
|---|---|---|
| Candida | 95.6% | 19 |
| Zymoseptoria | 95.1% | 98 |
| Rhizopus | 90.6% | 722 |
| Coprinopsis | 82.6% | 153 |
| Botrytis | 72.7% | 69 |
| N. crassa | 46.2% (93 baseline pairs) | 12 |

- The two methods give different dS for some pairs, because the alignments differ. Agreement within 1.5x or 0.02 was 97.7% (Candida), 89.0% (Zymoseptoria) and 90.8% (Rhizophagus families up to 50 genes). In Rhizophagus, 12.3% of baseline gated pairs were lost: their baseline dS median was 1.54 and their pair-level dS median was 3.3. I do not know which value is closer to the truth.
- NG86 truth classes partly define "close". The close-recall figures are therefore partly circular. The reviewer's count (2026-11-01): among baseline gated pairs that are in the "distant" class, 116 of 119 (Candida), 197 of 207 (Zymoseptoria) and 2,131 of 2,703 (Rhizopus) are distant only because identity < 0.4. That is the pre-filter rule itself. The independent evidence is the share of baseline gated pairs kept, the dS agreement with `--pairwise`, and the homeolog control. The identity >= 0.4 pre-filter also removes pairs that NG86 calls distant by definition, so the low distant share of the gated pairs (0.1 to 1.8%, N. crassa 10.9% on 55 pairs) is partly circular.

### 3.3 Families and structure

- **No family alignment, no tree, no node.** `wgd ksd` records a tree node for each pair. The pair-level output has no `node` column.
- **CRITICAL OPEN ISSUE. The peak-summary script needs `node`, and the choice of method changes the result.** On the Zymoseptoria pair-level table (gated pairs, 4,355 pairs in 168 families) the median dS is 0.039 with no averaging, 0.151 with UPGMA-on-dS node averaging and 1.55 with 1/(pairs in family) weights (reviewer calculation). In Candida the medians are 0.056, 0.020 and 0.014. S. pastorianus is insensitive, because 54% of its gated pairs are in 2-gene families. `merge_wgd_ks.py` writes an empty `node`, so the summary script would average each family to one value without an error. The pipeline default is now `pair`. Do not use pair-level output for peak statements until the method is chosen and validated. `analysis/WGD_PERFORMANCE_ANALYSIS/scripts/build_wgd_ksd_summary.py` node-averages pairs by `(family, node)` before it fits the mixture models, as `wgd mix` does. With no node, that step would average each family to one value. Do not run it on pair-level output until this is resolved. Options: skip averaging, weight each pair by 1/(pairs in family), or derive a node from a per-family tree. I have not tested any of them.
- **Redundancy.** A family of n genes gives n(n-1)/2 pairs. Without node averaging, large families contribute many more pairs than small ones. Histograms of dS then weight large families heavily.
- **Diamond hit cap.** The family-building step uses diamond defaults. Diamond returns at most 25 targets per query. In Rhizopus, 2,280 queries reached that limit. Families built from capped hits may miss some close relatives. This is unchanged from the old method.

### 3.4 Other limits

- Standard genetic code only. `TRANSL_TABLE` is not passed through. 1,154 genomes use table 12. The 1,154 count is for the full sample sheet of 23,680 genomes. The count among the 1,866 representatives is not known. `--table` and `--icode` exist in the script, but the module does not pass them. My check on Candida albicans showed standard-code behavior, but dS for other table-12 genomes is not validated.
- The pre-filter thresholds (identity 0.4, 30 codons) were chosen on the pilot genomes and not tested on separate genomes.
- One run per genome. No replicates. Task times were measured on the `intel` partition and include no contention estimate.
- `codeml` runs one pair at a time with default settings (CodonFreq 2, one omega). A `codeml` timeout of 60 s per pair marks the pair `timeout`. In the seven Nextflow-run genomes every pair has status `ok` or `prefiltered`. No pair timed out and none failed. Known silent-failure paths, none triggered in the tests: an empty families file raises an error in the pair generator, and the pipeline then drops the genome after two retries (default `ignore`); genes missing from the CDS file are counted only in the log; the `status` column is dropped at the merge step; one worker exception stops a genome; a `codeml` timeout is not recorded in the merged table.
- **dS = 0 pairs.** Identical proteins give dS = 0. Their share of the merged rows is 40.1% (Candida), 16.1% (Zymoseptoria), 13.2% (Botrytis), 2.1% (Rhizopus), 1.1% (Coprinopsis), 58.3% (Hortaea, 48 rows) and 0.0% (N. crassa, 55 rows). Report this share with any dS distribution.
- **`alignmentlength` meaning differs.** In pair-level output it is 3 x aligned columns of the pair. In `wgd ksd` it is the family alignment length. `wgd mix` removes pairs with `alignmentlength` below 300, which removes about 30% of the Candida gated pairs, 19% of Zymoseptoria and 5% of S. pastorianus (reviewer calculation).

## 4. What we did not test

- **Time for the full run.** Seven non-random genomes are not enough. A random stratified sample of about 30 representatives is needed first.
- **Task size.** Each task handles one genome and runs for 1 to 20 min. This is far below the 1 to 1.5 h target. Group about 20 genomes per task for production.

## 4a. Large families: Rhizophagus irregularis DAOM_197198 (measured 2026-11-01)

This genome has 95 baseline families over 50 genes (51 to 293 genes, 41 over 100). Together they hold 586,770 pairs. This is the only large-family test so far. It is one genome.

- **Pair-level `codeml` on all 586,770 pairs:** 5,951 s on 32 CPUs (0.010 s per pair). Every pair finished with status `ok`. No pair timed out.
- **Gate at dS <= 2.5 and >= 30 codons:** 251,561 pairs pass (42.9%). NG86 close-pair recall is 99.8%. The NG86-distant share of the gated pairs is 0.7% (partly circular).
- **Identity pre-filter (>= 0.4, >= 30 codons):** 73.7% of these pairs pass. The filter removes only 0.7% of the gated pairs (1,660). It reduces the pair count much less here than in the small-family genomes. I did not run the pre-filtered job, so the time saved by the filter in this genome is not measured.
- **`wgd ksd --pairwise` on a random sample of 8 of these families (67 to 162 genes, 49,056 pairs):** 16,543 s on 32 CPUs (0.337 s per pair). Its rate per pair is 33 times the pair-level rate above. This ratio compares two different pair sets and is indirect. At that rate, `--pairwise` on all 586,770 pairs would take about 55 h. Family-alignment cost grows faster than linearly with family size, so the true time may be higher.
- **Agreement on the 8-family sample (dS <= 2.5):** the pair-level run keeps 83.9% of the `--pairwise` gated pairs (85.2% at dS <= 2). It also passes 688 pairs that `--pairwise` does not. Among shared pairs, 94.5% agree within 1.5x or 0.02 in dS, with Spearman 0.905. These values are close to those from the Rhizophagus families up to 50 genes (87.6% kept, Spearman 0.943). The loss is therefore not worse in large families.
- **Not measured:** which method's dS is closer to the truth for the pairs that differ.

## 4b. Known-duplicate controls rerun with the pair-level method (measured 2026-11-01)

Same families as the old pilot (`wgd dmd`), same gate. Files: `eval/cap_recluster/controls_results.txt`.

**S. pastorianus homeologs (3,593 pairs, dS about 1.2, young set).** This is a hybrid of two parent species, not a WGD, and the reference is not independent.
- Pair-level method: 3,544 pairs (98.6%) pass dS <= 2.5 and >= 30 codons. The old default `ksd` passed 3,378 (94.0%).
- Homeolog dS quantiles (10/25/50/75/90%): 0.86 / 1.02 / 1.20 / 1.40 / 1.62. The earlier prototype gave 0.86 / 1.03 / 1.20 / 1.41 / 1.66. dS values for the same pairs are identical to the prototype (median difference 0.0000).
- Genome-wide gated dS (5,344 pairs, new; 4,884 pairs, old default `ksd`): median 1.16 in both, histogram mode 1.15 in both. This agreement is weak evidence, because the 3,544 homeolog pairs are 66% of the gated pairs. Spearman 0.961 on 4,657 shared pairs. Within 1.5x or 0.02: 98.1%.

**S. cerevisiae SX6 genome-wide gated dS (dS <= 2.5).** This is a small set: 347 pairs (new) and 286 pairs (old). Histogram mode 0.25 in both. Median 0.75 (new) against 0.90 (old). Spearman 0.969 on 228 shared pairs; within 1.5x or 0.02: 90.8%. 119 pairs appear only in the new result and 58 only in the old.

**YGOB ancient-WGD ohnologs in SX6 (369 mapped pairs, 334 in a family).** The ancient yeast WGD is NOT recovered by either method.
- Gate at dS <= 2.5: 3 pairs pass (0.8%) with the pair-level method, 4 (1.1%) with the old method. At dS <= 5: 49 (13.3%) and 58 (15.7%).
- Pair-level dS of the 293 ohnologs that were run: 10/25/50/75/90% = 4.37 / 6.26 / 10.25 / 64.61 / 71.43. These values are saturated. The pairs have low protein identity, often about 25%.
- The dS <= 2.5 cut therefore hides old duplications, as it did in the old pipeline. Only younger events (roughly dS below 2) are visible. This limit is shared with the old method and is not new.

Summary: for young duplicates, the pair-level method gives the same dS peak position as the old method. For ancient WGD it recovers nothing, and neither did the old method.

## 5. Methods tested and rejected

- **Family re-clustering (`cov50_I2`).** It saves 1.8x to 14x in time, but lost 12 to 35 points of NG86-close pairs: recall 88.0% (Candida), 87.4% (Zymoseptoria), 75.5% (Rhizopus), 88.4% (Botrytis), 68.8% (Coprinopsis). Not recommended.
- **Hard cap on family size.** Cap at 50 genes: close recall 56.5% (Zymoseptoria), 39.1% (Rhizopus), 68.8% (Coprinopsis). Not recommended.
- **Pre-filter on diamond hit statistics.** 25 to 26% of gated pairs in Zymoseptoria and Rhizopus have no direct diamond hit. Not usable.

## 6. How to look deeper for one genome

Run the same genome with `--wgd_ks_method ksd_pairwise` (family alignment, per-pair `codeml`), or with `ksd`. Results go to separate folders. Use `--wgd_max_ds 2` or a lower `wgd_prefilter_identity` to widen or narrow the set.

## 7. Changes after the second review (2026-11-02)
- Node: `wgd_pair_ks.py` writes a `node` column. It is the last common ancestor in a neighbour-joining tree built on protein identity distance (1 - identity), midpoint-rooted. On simulated families with a known tree the node-mean dS distribution is close to the truth (Wasserstein-1 0.045-0.057; 0.41-0.43 with no averaging; FastTree on the family alignment 0.007). Node labels agree only partly with the true tree (ARI 0.55-0.62). The histogram mode of node means is unstable. Source: `eval/sim_controls/report.md`.
- Failures are loud: empty families file, more than 1% missing genes, more than 5% failed codeml runs and worker errors stop `wgd_pair_ks.py`. A genome that fails inside a bundle is written to `bundle_<locustag>.failed.txt`. The merge writes `tables/wgd.ks.check.txt` and the workflow logs an error when it is not OK (missing genomes, failure lists). The `status` column is kept in the merged table.
- Genetic code: `TRANSL_TABLE` is passed (`--table`, codeml `icode 8` for table 12). Other non-standard tables are skipped and listed as failed unless `wgd_allow_untested_table=true`. On the 7-genome test, table 12 changed Candida: 3,477 rows (3,529 before), and 10.9% of common pairs differ in dS by more than 0.05. Zymoseptoria was identical.
- Bundling: `wgd_pair_bundle` genomes per task (default 20, a first guess; 4 in the test). Output is flat in `wgd_ksd_pp/`.
- Switch: `wgd_speedups=false` turns off the pre-filter, dS cut and pair-codon gate (output `wgd_ksd_pp_full/`). The `false` path is not yet tested.
- Bug found by the test run: the NJ code failed on families with identical genes (all-zero distances). Fixed.
