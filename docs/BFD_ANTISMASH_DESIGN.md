# antiSMASH-lite in BFD_FUNCTIONAL: data capture design

Status: DRAFT, 2026-11-22. `run_antismash` stays false by default. Not run at scale.
Source of numbers: interim report from the `ebc-dgf` session (32 genomes, antiSMASH 8.0.4,
MIBiG 4.0). Labels: MEASURED = counted on those 32 genomes. EXTRAPOLATION = scaled, not run.
PROPOSED = design choice, not yet tested. I did not re-measure these numbers.

## Goal
Keep cluster location, member genes, MIBiG hits, product types and cluster stats for about
24,000 genomes. Keep run time and disk use low. Never store raw antiSMASH output.

## Result layout (decided)
- Key: LOCUSTAG, not `out`. This is the same key as `results/function/cazy` and
  `results/function/pfam_hmmscan`. `out` is rebuilt from species and strain, so a rename
  would orphan results. Note: `README_DATA_LAYOUT.md` calls `out` the per-genome folder key.
  That applies to `genome_annotation/`, not to `results/function/`.
- Path: `results/function/antismash/<sha1(LOCUSTAG)[:2]>/<LOCUSTAG>.antismash.<table>.tsv.gz`.
  Bucket width 2 is set in `modules/common/utils.nf` and `bin/genome_stats_paths.py`.
- Files per genome: `regions`, `cds`, `known`, `known_pairs` (tsv.gz) and `run.json`.
  Size 21-46 KB gz per genome (MEASURED, 3 genomes).

## Tables (MEASURED on 32 genomes, sorted by LOCUSTAG, zstd parquet)
| Table | Grain | Rows (32 genomes) | Size |
|---|---|---|---|
| bgc_cluster | one per region | 2,066 | 51.6 KB |
| bgc_cluster_gene | one per gene in a region | 40,435 | 390 KB |
| bgc_mibig_hit | top-3 MIBiG hits per region | - | 25 KB |

- bgc_cluster: locustag, record, region, start, end, length, products, categories,
  contig_edge, core_start, core_end, candidate_kinds, n_cds, n_biosynthetic, n_additional,
  best_mibig, best_similarity, best_rank.
- bgc_cluster_gene: locustag, record, region, gene_id, start, end, is_minus, gene_kind,
  functions (smCOG text). 40,435 of 40,435 gene_id values matched `gene_info.parquet`.
- bgc_mibig_hit: similarity, core_gene_hits, gene_hits, n_entry_proteins, synteny_score,
  blast_score. `known_pairs` (our gene to MIBiG gene, identity, coverage) is kept as a TSV.
- SQL: `DeltaGain_Fungi_wt_realdata/analysis/bgc_subset/build_bgc_parquet.sql`
  (branch `pathway-census-realdata`, not pushed).
- Query speed: NRPS clusters with best similarity >= 75 took 1-3 ms (217 of 2,066).
  Synthetic test, 32 genomes repeated 750 times in memory: 9-62 ms. This tests speed only.
- Size at 24,000 genomes: about 350 MB (EXTRAPOLATION, unverified).

## Rules for users of the data
- AntiSMASH `rank` is not similarity order. Use the highest similarity of the top 3.
- 188 of 430 regions with best similarity >= 75 sit on single-protein MIBiG entries.
  Filter on `n_entry_proteins`.
- smCOG text is in `cds.functions`. Pfam domains are not captured (see below).

## Run options and cost
- Lite flags: `--taxon fungi --cb-general --cb-knownclusters`, plus the gene-source option
  below. Dropped from the old module: `--fullhmmer`, `--clusterhmmer`, `--pfam2go`.
- Af293: lite found the same 44 regions as full (same coordinates and classes), 786 s vs
  6,535 s. Other genomes not compared.
- 32 genomes, 8 CPUs (MEASURED): wall time median 965 s, mean 1,343 s, max 4,408 s.
  Peak memory 8.7-12.4 GB. Raw antiSMASH output median 264 MB, deleted after extraction.
- 24,000 genomes: about 71,600 CPU hours (EXTRAPOLATION). The 32 genomes are all
  Pezizomycotina. Yeasts, basidiomycetes and Mucoromycota are untested.
- Bundling (PROPOSED): about 4 genomes per job, 1-1.5 h per job. Wall time varies 8-73 min,
  so bundle by contig count or stop starting genomes after a time budget.
- Raw output stays in `$SCRATCH`. Only the extracted tables go to /bigdata.

## Decisions pending (2026-11-22)
1. Output trimming: add `--no-zip-output --no-summary-gbk --no-region-gbks`.
   Test on A1163 (job 32375093) is running. No time or size result yet.
2. Gene source: `--genefinding-gff3` (FASTA + GFF3, measured) vs the funannotate GBK
   (`predict_results/<out>.gbk.gz`, not tested). See the next section.
3. Extra JSON fields for the minimal set (PROPOSED, not yet checked against a real JSON):
   antiSMASH version and per-module timings, if present at top level. Check one JSON first.
4. Shared checkout: `ebc-dgf` committed antiSMASH changes to `unify-funannotate-sif`
   (75dc089, c5d2407). Check that no live Nextflow run reads a frozen copy.

## Gene source: GFF3 or GBK
- AntiSMASH for fungi takes either GenBank (uses its CDS annotation) or FASTA. With FASTA it
  needs a gene source: its own gene finder, or `--genefinding-gff3`. So the flag is needed
  only because the module feeds FASTA.
- GBK needs no flag and no GFF3. It is already in `predict_results/`. It is larger (about
  21 MB gz for a small genome) and must be inflated.
- GFF3 route side effect: the extractor reads the `gene` qualifier in place of
  `locus_tag` (`q1()` in `extract_antismash_subset.py`). The 40,435-of-40,435 join shows it works.
- Not measured: time and result differences between the two routes.

## Update 2026-11-22 (second report from `ebc-dgf`; full tables in
`DeltaGain_Fungi_wt_realdata/analysis/bgc_subset/FINDINGS.md`, branch pathway-census-realdata, d546551, not pushed)

### Where the time goes (MEASURED)
- Per-module time on the slowest genome (N. sphaerica CBS 166.26): clusterblast 3,660 s of
  3,859 s module time (94.8%). genefunctions 1.7%, hmm_detection 1.5%.
- A1163 without `--cb-general` (`--cb-knownclusters` only): 233 s vs 881 s (951 s on a repeat),
  3.8x-4.1x faster. Regions identical (43 = 43). `known` and `known_pairs` (top 3) identical.
  Only the `blast_top` table needs `--cb-general`, and BFD does not store it.
- Confirmed on 3 genomes (MEASURED, same flags, compared with the earlier run that had it):
  A1163 881 s -> 233 s (3.8x); A. niger CICC 2151 3,896 s -> 562 s (6.9x);
  N. sphaerica CBS 166.26 4,308 s -> 807 s (5.3x). Regions, `known` (top 3) and
  `known_pairs` (top 3) identical on all 3. Speed-up is larger for slow, fragmented genomes.
- Bound only (EXTRAPOLATION, not a result): if 3.8x-6.9x held on average, the 24,000-genome
  cost would fall from about 8,950 to roughly 1,300-2,400 8-CPU hours. Mean time without
  `--cb-general` is not measured.
- PROPOSED: delete `--cb-general` from the module. Not yet edited. Test one yeast, one
  basidiomycete and one Mucoromycete first.
- Not tested: non-Pezizomycotina genomes.

### Output-trimming flags (MEASURED, A1163, job 32375093)
`--no-zip-output --no-summary-gbk --no-region-gbks`: raw directory 102 MB vs 206 MB.
Wall time 951 s vs 881 s, no saving (node-to-node spread is 10-15%). Regions identical.
These flags save scratch space, not time.

### Contig filter (MEASURED, 33 genomes, 2,140 regions)
- Regions on contigs < 5 kb: 0. < 10 kb: 3 (0.1%), all on a contig edge. < 20 kb: 6 (0.3%),
  all on a contig edge. Shortest contig with a region: 5.9 kb.
- Wall time follows contig count (Spearman 0.80), not genome size (0.11). But the second
  slowest genome has only 141 contigs, so short contigs do not explain the time.
- Slowest genome, 4 runs: control 4,308 s (72 regions); 5 kb prefilter 4,742 s (72 regions);
  10 kb prefilter 4,061 s (71, lost one 5.9 kb terpene on a contig edge);
  `--minlength 10000` 3,708 s (71). Differences are inside the node noise.
- `antismash --minlength` exists in 8.0.4, default 1000.
- PROPOSED: no contig filter. `filter_contigs.py` exists in `bgc_subset` if wanted.

### Top 3 vs top 5 MIBiG hits (MEASURED, 12 genomes, 310 regions with a hit)
- Best-similarity hit outside ranks 1-3: 13 (4.2%). Outside ranks 1-5: 5 (1.6%).
- Rank 4-5 hit with similarity >= 75 beating ranks 1-3: 4 (1.3%); 3 of the 4 are on
  one-protein MIBiG entries.
- Size per genome (gz TSV): known 1,862 B (top 3), 2,247 B (top 5); known_pairs 7,445 B and
  8,809 B. Two tables at 24,000 genomes: about 165 MB (top 3), 191 MB (top 5)
  (EXTRAPOLATION). Mixed (known 5, pairs 3) saves little.
- Choice of TOPN is the owner's decision. The Fungi_BFD extractor is unchanged (TOPN = 3).

### Change made 2026-11-22 (user approved)
`modules/BFD/ANTISMASH/main.nf`, lite flags are now
`--cb-knownclusters --no-zip-output --no-summary-gbk --no-region-gbks` (no `--cb-general`).
Full mode is unchanged. Uncommitted. Not yet run through Nextflow. Still untested without
`--cb-general`: yeast, basidiomycete, Mucoromycete. TOPN stays 3 until the owner decides.

### Decision 2026-11-22: top 5 MIBiG hits (owner)
`TOPN = 5` in `nextflow/bin/extract_antismash_subset.py` (uncommitted). Applies to `known` and
`known_pairs`. Size at 24,000 genomes for the two tables: about 191 MB (EXTRAPOLATION).
Comments that say "top-3" in the extractor and module are now out of date.

### Genetic code (TRANSL_TABLE), added 2026-11-22 (uncommitted)
- Problem (MEASURED by `ebc-dgf`, D. hansenii CBS767, table 12): antiSMASH 8.0.4 reads the
  GFF3 as code 1. 3,923 of 5,951 CDS (65.9%) differ from the BFD protein; all changed
  residues are Ser in BFD and Leu in antiSMASH. 1,164 BFD genomes use table 12, 1 uses 26.
- Fix: `RUN_ANTISMASH` appends `;transl_table=<N>` to every CDS line of the GFF3 when
  `meta.transl_table` is not 1. With the tag, 5,951 of 5,951 translations matched BFD.
- Code: `workflows/BFD.nf` adds `transl_table` to `meta` (from `samples.csv` TRANSL_TABLE,
  default '1'); `modules/BFD/ANTISMASH/main.nf` runs the awk step.
- Checked: the Groovy string renders the expected awk command, and the awk output on a
  3-line GFF3 is correct. NOT run through Nextflow or antiSMASH yet.
- Table 26 is untested. The effect of the wrong code on regions and MIBiG hits is not
  measured (D. hansenii regions are the same either way; it has no MIBiG hits).
- Side effect: `meta` now has one more key for every tool that takes it. Most BFD function
  modules use `storeDir`, so existing results should not rerun, but a `-resume` of a
  non-storeDir step may. Not checked.

### Nextflow test and merge test, 2026-11-22
- Test run `.nf_launch/20261122_antismash_lite_nf_test` (job 32377249): 3 genomes through
  `BFD_FUNCTIONAL` (B. adeninivorans F634A7E4, D. hansenii F8C91A38 table 12, M. hiemalis
  F00CDB58). 3 of 3 completed. Task time 56-151 s, peak memory 1.6-3.6 GB (MEASURED).
  Regions 4/5/11, `known` 0/0/1, `known_pairs` 0/0/1: same counts as the hand runs.
  Buckets match sha1(LOCUSTAG)[:2]. Top-5 limit not exercised (no region had > 1 hit).
- `merge_antismash.py` run by hand on the 3 `regions` tables: 20 rows (4+5+11), parquet
  written with duckdb 1.1.3, counts match. Not run through the Nextflow `MERGE_ANTISMASH`.
- GAP: `MERGE_ANTISMASH` merges only the `regions` table, into `antismash.parquet`, with
  10 columns (species_prefix, record_id, cluster_num, product, contig_edge, start, end,
  categories, n_cds, n_biosynthetic). It drops length, core_start, core_end,
  candidate_kinds, n_additional and the best MIBiG fields. It does not merge `cds`, `known`
  or `known_pairs`. So `bgc_cluster_gene` and `bgc_mibig_hit` are not built by the pipeline
  yet. The SQL that builds them (`build_bgc_parquet.sql`) is on an unpushed DeltaGain branch.

### Merge rewritten to schema v2, 2026-11-22 (uncommitted) -- closes the GAP above
- Schema source: `DeltaGain_Fungi_wt_realdata/analysis/bgc_subset/build_bgc_parquet.sql`
  (commit dcc239a, local, NOT pushed by `ebc-dgf`; it needs its owner's approval).
- `bin/merge_antismash.py` now reads the regions, cds, known and known_pairs tables and
  writes `bgc_cluster`, `bgc_cluster_gene`, `bgc_mibig_hit`, `bgc_mibig_pair` (zstd parquet,
  same columns, types and sort order as the SQL). Join key: (locustag, record, region).
  `gene_id` joins to `gene_info.gene_id`. Hits above `--topn` (5) are dropped.
  Header-only files are fine. Empty gzip files (stub output) are skipped.
- `modules/BFD/MERGE_ANTISMASH/main.nf` emits the four parquets. The coverage check uses
  `--key locustag` (only when `merge_all` is on). `BFD_FUNCTIONAL.nf` passes all four
  table kinds to the merge (mix of the RUN_ANTISMASH outputs, or a brace glob in glob mode).
- The old `antismash.parquet` (10 columns, regions only) is no longer written.
- Test (MEASURED, 3 genomes): standalone and through Nextflow (job 32377318, 2 tasks, 0
  failed). Rows: bgc_cluster 20, bgc_cluster_gene 314 (= 61+78+175 cds rows),
  bgc_mibig_hit 1, bgc_mibig_pair 1. bgc_cluster from Nextflow equals the standalone
  result. Column types match the schema (e.g. region USMALLINT, contig_edge BOOLEAN).
- Not tested: glob mode (`merge_all`), the coverage check, the stub block, top-5 with
  several hits per region, more than 3 genomes, a genome with no regions.
- Expected size at 24,000 genomes (`ebc-dgf`, EXTRAPOLATION from 12 genomes): about 551 MB
  for the four tables.

## TODO: check back on regions with several MIBiG hits (opened 2026-11-22)
The 3-genome Nextflow test had only 1 MIBiG hit in total, so no path that handles several
hits per region has run through the pipeline. Check these places before a production run:
1. `bgc_cluster.best_mibig / best_similarity / best_rank` keep ONE hit per region
   (highest similarity among ranks 1-5; ties go to the lower rank). The other hits are only
   in `bgc_mibig_hit`. A query that reads `bgc_cluster` alone hides them.
2. Rank is antiSMASH order, not similarity order. The best-similarity hit was outside
   ranks 1-3 in 4.2% of regions and outside ranks 1-5 in 1.6% (12 genomes, `ebc-dgf`).
   Hits at rank 6 and above are never stored.
3. `bgc_mibig_pair` is filtered to rank <= 5 and keyed by (locustag, record, region, rank,
   gene_id). Check that joins to `bgc_mibig_hit` do not multiply rows.
4. `region` is unique only inside a record, and `record` only inside a genome. A join
   without `locustag` merges hits of different genomes.
5. One-protein MIBiG entries (`n_entry_proteins = 1`) give high similarity with weak
   evidence. A top-hit filter should also use `n_entry_proteins`.
How to check: run the merge on the 12 `ebc-dgf` genomes (mean 3.2 hits per region, 310
regions with a hit) and compare row counts of `bgc_mibig_hit` and `bgc_mibig_pair` with its
778 and 3,110 rows. Then test the joins in points 3 and 4 with a query.

### SQL and database registration, 2026-11-22
- `nextflow/sql/merge/antismash.sql` is now the primary definition of the four `bgc_*`
  tables (columns, types, sort order, best-hit rule). `bin/merge_antismash.py` reads it and
  fills the placeholders. `ebc-dgf`'s `build_bgc_parquet.sql` (DeltaGain, unpushed) is not
  used and can be dropped. Test: output of the new script equals the earlier output
  row for row on the 3 test genomes (standalone and through Nextflow, job 32377347).
- `nextflow/bin/build_BFD_duckDB.sh`: `mkview` lines for the four tables (optional, skipped if
  the parquet is missing) and uniqueness checks on (locustag, record, region[, gene_id | rank]).
  `bgc_mibig_pair` has no uniqueness check, because one gene can pair with several MIBiG genes.
- `sql/bfd_compat_views.sql`: `COMMENT ON VIEW` text for each table.
- `nextflow/bin/schema_doc.py`: knows the lower-case `locustag` key, so `docs/BFD_SCHEMA.md`
  will list it after the next database rebuild. That file is generated and was not edited.
- Tested on a small test database with views over the 3-genome parquet: views created,
  uniqueness queries pass, comments applied, `schema_doc.py` runs. NOT tested: a full
  `build_BFD_duckDB.sh` run, and `verify_BFD_duckDB.py` with the new views.
- `sql/table_schema.json` covers only species, asm_stats and busco_genome, so the bgc tables
  are not added to it.
