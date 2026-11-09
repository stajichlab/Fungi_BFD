# BFD database: introduction

State on 2026-11-07. This page describes the files and tables that exist now.
Items that are not finished are listed under "Known gaps".

The BFD database is a set of Parquet tables plus one DuckDB file. It holds the
gene models, assembly statistics, and protein function predictions for the
fungal genomes annotated in the BFD funannotate waves.

## 1. Where the files are

All paths are relative to the run folder:

```
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/
```

| Path | What it holds |
|---|---|
| `tables/` | **The database.** One Parquet file per table, built by the `MERGE_*` steps of the Nextflow pipeline. |
| `tables/_manifest/` | One JSON file per table: row count, build time, source, sha256. |
| `tables/_superseded_YYYYMMDD/` | Older versions of tables, kept until they are deleted. Do not query these. |
| `db/BFD.duckdb` | A DuckDB file that names the tables in `tables/` (see section 3). |
| `seqdb/` | Protein and CDS sequences, bgzip-compressed and indexed with `samtools faidx`. |
| `genome_annotation/<out>/predict_results/` | The funannotate results per genome (GFF3, GenBank, FASTA). These are the source of the gene tables. |
| `results/function/<tool>/` | Per-genome outputs of the function tools (Pfam, MEROPS, CAZy, SignalP, TMHMM, WoLF PSORT, PredGPI, TargetP, AIUPred, SwissProt, Infernal). These are the inputs to the `MERGE_*` steps. |
| `results/genome_stats/<type>/` | Per-genome outputs for assembly stats, BUSCO, telomeres, gene stats, codon and amino acid frequency. |
| `samples.csv`, `suppress.txt` | The genome list for the merges. A genome is merged when it is in `samples.csv` and not in `suppress.txt`. |

The pipeline code is in `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow`.
The merge SQL is in `nextflow/sql/merge/`, and the DuckDB build script is
`nextflow/bin/build_BFD_duckDB.sh`.

## 2. Keys

Each genome has a `LOCUSTAG`, an 8-character hex code such as `F77FFF1D`.
Every gene, transcript, and protein ID starts with the LOCUSTAG of its genome.

| Key | Example | Parent table |
|---|---|---|
| `LOCUSTAG` | `F77FFF1D` | `species` |
| `species_prefix` | `F77FFF1D` (the LOCUSTAG, in the per-gene and per-protein tables) | `species.LOCUSTAG` |
| `ASMID` | the assembly accession | `species`, `asm_stats` |
| `gene_id` | `F77FFF1D_000001` | `gene_info` |
| `transcript_id` | `F77FFF1D_000001-T1` | `gene_transcripts` |
| `protein_id` | `F77FFF1D_000001-T1` | `gene_proteins` |

The folder name `<out>` in `genome_annotation/` (for example
`Aaosphaeria_arxii_CBS_175.79`) maps to the LOCUSTAG through `samples.csv`
and `seqdb/genomes.tsv`.

## 3. How to query

Two ways give the same data.

**A. The DuckDB file.** It gives every table a short name and adds summary views.

```bash
module load duckdb
duckdb -readonly /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/db/BFD.duckdb
```

```sql
SELECT * FROM v_species_summary WHERE GENUS = 'Aspergillus';
SELECT * FROM pfam WHERE protein_id = 'F77FFF1D_000002-T1';
```

Always open it with `-readonly`. More than one reader can then use it at the same time.

**B. The Parquet files directly.** This needs no DuckDB file, and it works in
DuckDB, Python (pyarrow, polars, pandas), and R (arrow).

```sql
SELECT count(*) FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/pfam.parquet')
WHERE protein_id >= 'F77FFF1D_' AND protein_id < 'F77FFF1E';   -- one genome
```

**Use `=` or a range on the id column, not `LIKE` or `starts_with`.** The
sorted tables store the minimum and maximum id of each row group, so DuckDB
skips the row groups that cannot match. DuckDB 1.1.3 skips row groups only for
`=`, `<`, `>` and `BETWEEN`. Measured on `pfam` (368 million rows, 2026-11-08):

| Query | Time |
|---|---|
| `protein_id = 'F07D2B4E_004856-T1'` | 0.6 s |
| `protein_id >= 'F07D2B4E_' AND protein_id < 'F07D2B4F'` (17,166 rows) | 0.6 s |
| `starts_with(protein_id, 'F07D2B4E_')` (same rows) | 28.8 s |

For the range, the upper bound is the LOCUSTAG with its last character
increased by one (`...4E` becomes `...4F`; `...9` becomes `...A`).

**Fast lookups (table macros).** SwissProt hits are stored as 256 bucket
files. A `WHERE protein_id = ...` filter on the `swissprot_hits` view opens all
256 files (7-10 s per protein). These macros open only the one bucket file for
the given protein or genome:

```sql
SELECT * FROM swissprot_hits_for('F77FFF1D_000002-T1');        -- one protein
SELECT * FROM swissprot_hits_for('F77FFF1D');                  -- one genome (LOCUSTAG)
SELECT * FROM swissprot_model_fungi_for('F77FFF1D_000002-T1'); -- best hit per model fungus
```

They return in under 1 s, including DuckDB start-up. `swissprot_best` (one row
per protein) is a single sorted file and is fast through its view.

Do not filter `v_swissprot_model_fungi` by genome (`species_prefix`). The filter
is not passed through the per-protein ranking, so DuckDB ranks all 2.2 billion
hits first. On a 2026-11-08 test that ran out of memory on a login node. Use
`swissprot_model_fungi_for('<LOCUSTAG>')` instead.

**Note on the DuckDB file.** `db/BFD.duckdb` holds views that read
`tables/*.parquet` by absolute path, plus the small `swissprot_subject` table
and the macros above. It works only where `tables/` is at the same path. The
pipeline publishes it to `BFD_function/db/BFD.duckdb`;
`Fungi_BFD_runs/db/BFD.duckdb` is a link to that file. Lookups are fastest when
a table is sorted by the column you filter on (see section 8).

## 4. Tables

Row counts come from the 2026-11-06 build. The full column lists are in
`docs/BFD_SCHEMA.md`, which is generated from the database by
`nextflow/bin/schema_doc.py`.

### Genome level

| DuckDB name | Parquet file | Rows | Content |
|---|---|---|---|
| `species` (view) | `species.parquet` | 23,680 | Taxonomy (phylum to species, strain), NCBI taxon ID, BUSCO lineage, ASMID, LOCUSTAG |
| `asm_stats` (view) | `asm_stats.parquet` | 23,584 | Assembly size, contig count, N50/L50, GC %, masked %, gaps, telomere counts |
| `busco_genome` (view) | `busco_genome.parquet` | 23,582 | BUSCO completeness of the assembly |
| `telomere_summary` | `telomere_summary.parquet` | 23,582 | Telomere repeats per assembly |
| `telomere_tracts` | `telomere_tracts.parquet` | 73,242 | One row per telomere tract |
| `aa_frequency` | `aa_freq.parquet` | 477,535 | Amino acid frequency per genome |
| `codon_frequency` | `codon_freq.parquet` | 1,500,925 | Codon frequency per genome |

### Gene models

| DuckDB name | Parquet file | Rows | Content |
|---|---|---|---|
| `gene_info` | `gene_info.parquet` | 249,546,610 | Gene coordinates and type |
| `gene_transcripts` | `gene_transcripts.parquet` | 249,546,610 | Transcript coordinates, partial flags, CDS length, md5 |
| `gene_proteins` | `gene_proteins.parquet` | 244,988,294 | Protein ID, length, md5 |
| `gene_exons` | `gene_exons.parquet` | 758,181,061 | Exon coordinates and GC |
| `gene_CDS` | `gene_CDS.parquet` | 753,622,745 | CDS segment coordinates and GC |
| `gene_introns` | `gene_introns.parquet` | 508,634,451 | Intron coordinates, splice sites, phase (no sequence; see section 5) |
| `gene_trna` | `gene_trnas.parquet` | 4,558,316 | tRNA genes |
| `gene_intergenic_distances` | `gene_intergenic_distances.parquet` | 233,406,624 | Distance between neighboring genes |

The gene tables do not hold protein or CDS sequences. Get sequences from `seqdb/` (section 5).

### Protein function

| DuckDB name | Parquet file | Rows | Content |
|---|---|---|---|
| `pfam` | `pfam.parquet` | 368,209,721 | Pfam domain hits (hmmscan domtblout) |
| `merops` | `merops.parquet` | 43,462,876 | MEROPS peptidase hits |
| `cazy` | `cazy.cazymes_hmm.parquet` | 9,164,068 | dbCAN HMM hits |
| `cazy_overview` | `cazy.overview.parquet` | 9,310,657 | dbCAN consensus CAZyme calls |
| `signalp` | `signalp.signal_peptide.parquet` | 18,767,284 | SignalP 6 signal peptides (positive calls only) |
| `tmhmm` | `tmhmm.parquet` | 45,707,260 | TMHMM 2 helix counts and topology |
| `wolfpsort` | `wolfpsort.parquet` | 244,790,951 | WoLF PSORT localization |
| `predgpi` | `predgpi.parquet` | 227,238,177 | PredGPI GPI-anchor predictions |
| `targetp` | `targetP.parquet` | 178,535 | TargetP 2 (only 127 genomes; see "Known gaps") |
| `idp` | `idp.parquet` | 233,329,228 | AIUPred disordered regions |
| `idp_summary` | `idp_summary.parquet` | 144,924,594 | Disordered fraction per protein |

### Paralog Ks

| DuckDB name | Parquet file | Rows | Content |
|---|---|---|---|
| `wgd_ks` | `wgd.ks.parquet` | 25,831,001 | Paralog pair Ks from the paralogoscope workflow |
| `wgd_ks_summary` | `wgd.ks.summary.parquet` | 7,384 | Per-genome Ks summary |

### Summary views

| View | Content |
|---|---|
| `v_species_summary` | One row per genome: taxonomy, assembly stats, gene count, and counts of proteins with Pfam, CAZy, MEROPS, SignalP, TMHMM, and high disorder |
| `v_protein_annotation` | One row per protein: annotation flags from all function tables |
| `v_telomere_both_ends` | Scaffolds with a terminal telomere at both ends |

## 5. Sequences (`seqdb/`)

Protein and CDS FASTA files are copied for 23,412 genomes (217 GB). Each copy is
bgzip-compressed and indexed (`.fai` and `.gzi`). Its content was checked
against the original by md5.

```
seqdb/proteins/<bucket>/<LOCUSTAG>.proteins.fa.gz
seqdb/cds/<bucket>/<LOCUSTAG>.cds.fa.gz
```

`<bucket>` is the first two hex characters of the SHA-1 of the LOCUSTAG:

```bash
LT=F01BDC31
B=$(echo -n $LT | sha1sum | cut -c1-2)
module load samtools
samtools faidx seqdb/proteins/$B/$LT.proteins.fa.gz ${LT}_000001-T1
```

**Genome and intron sequences.** Each genome FASTA is bgzip-compressed and indexed
in place: `genome_annotation/<out>/predict_results/<out>.scaffolds.fa.gz` (+ `.fai`,
`.gzi`). An intron sequence is a region of it, from the `gene_introns` coordinates
(1-based, inclusive; use `-i` for strand -1):

```bash
samtools faidx genome_annotation/<out>/predict_results/<out>.scaffolds.fa.gz chrom:start-end      # strand 1
samtools faidx -i genome_annotation/<out>/predict_results/<out>.scaffolds.fa.gz chrom:start-end   # strand -1
```

`seqdb/genomes.tsv` lists LOCUSTAG, `<out>`, and the source protein and CDS paths.
It is built by `bin/seqdb_build.sh`.

## 6. Known gaps

- **Infernal (Rfam)**: the Infernal steps have not run for this set.
  `tables/infernal_rfam.parquet` is from the old parser, which read the wrong
  columns. Do not use it.
- **TargetP** has results for 127 genomes only. A move to DeepLoc is being considered.
- **DeepTMHMM** is not part of this release. Use `tmhmm`.
- Not converted to the new merge path yet: EggNOG, Phobius, TCDB, antiSMASH,
  tRNAscan-SE, and CAZy CGC. These have no BFD-wide tables yet.
- Genome counts differ between tables: `species` has 23,680 rows, `gene_info`
  covers 23,452 genomes, and `seqdb` covers 23,412. Genomes without finished
  predictions have no gene or function rows.

## 7. How the tables are built

1. Per-genome tools write results to `results/function/` and `results/genome_stats/`.
2. Each `MERGE_*` step reads the per-genome files for the genomes in
   `samples.csv` minus `suppress.txt`, and writes one Parquet file to `tables/`.
   Most merges use DuckDB SQL templates (`nextflow/sql/merge/*.sql`, run by
   `nextflow/bin/duckdb_merge.py`); the genome-stats merges use
   `nextflow/bin/manifest_to_parquet.py`. Both read the input files in file-name
   order (`<LOCUSTAG>.<type>...`), so rows come out grouped and roughly sorted by
   LOCUSTAG and id.
3. `BUILD_DUCKDB` runs `nextflow/bin/build_BFD_duckDB.sh` to make `db/BFD.duckdb`
   (views, the `swissprot_subject` table, the lookup macros, duplicate-id checks).
   It waits until every Parquet file of this run is fully copied to `tables/`.
   Each merge passes all of its output files to `BUILD_DUCKDB` for this check.
4. A merge whose inputs have not changed is skipped (`merge_if_stale`). A new
   table that covers fewer genomes than the old one is not published unless
   `allow_shrink` is set.

**Forcing a rebuild after a code change.** `merge_if_stale` compares file
modification times only. A change to a merge's code (for example a new sort
order) does not make its table stale. Name the merges to rebuild:

```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/BFD_function
sbatch run_BFD_merge.sh --force_merges MERGE_PFAM,MERGE_GENE_STATS
```

`skip_merges` still wins over `force_merges`.

**Cached tasks after a code change.** `run_BFD_merge.sh` uses `-resume`.
Nextflow hashes a task's script text and inputs, but not the helper files the
script calls by path (`bin/*.py`, `bin/*.sh`, `sql/merge/*.sql`). Each merge
script therefore has a comment line with an MD5 of its helper files
(`codeMd5()` in `modules/common/utils.nf`). A helper change gives a new task
hash, so a forced merge runs again and does not reuse the old output. `WORKTREE=<path> sbatch
run_BFD_merge.sh ...` runs the code in another worktree (for testing a branch).

## 8. What a pipeline rebuild applies, and what it does not

State on 2026-11-08.

**Built into the Nextflow pipeline (a rebuild applies these):**

| Item | Where |
|---|---|
| DuckDB merge templates for Pfam, MEROPS, SignalP, TMHMM, WoLF PSORT, PredGPI, TargetP, TF, CAZy, SwissProt, Infernal | `sql/merge/*.sql`, `bin/duckdb_merge.py` |
| Pfam sorted by `protein_id` | `sql/merge/pfam.sql` (`ORDER BY`) |
| SwissProt hits sorted by `protein_id`, rank; 256 bucket files; best-hit and summary tables | `sql/merge/swissprot_*.sql`, MERGE_SWISSPROT(_GATHER) |
| Input files read in file-name (LOCUSTAG) order, so row groups can be skipped on lookups | `bin/duckdb_merge.py`, `bin/manifest_to_parquet.py` |
| No sequence columns: `CDS_sequence`, `peptide`, intron `seq` | MERGE_GENE_STATS (`--exclude`) |
| View-only DuckDB file, duplicate-id checks, SwissProt lookup macros | `bin/build_BFD_duckDB.sh` |
| New predictions write `<out>.scaffolds.fa.gz` + `.fai` + `.gzi`; genome readers accept `.gz` | FUNANNOTATE_PREDICT, `resolveGenomeFile()` |
| `force_merges` for code changes | `modules/common/utils.nf` (`mergeForced`) |
| Helper code (`bin/`, `sql/merge/`) is part of the task cache key | `codeMd5()` in `modules/common/utils.nf`, one comment line per merge module |
| ZSTD compression for every merge output | `sql/merge/*.sql`, `bin/manifest_to_parquet.py` |
| `BUILD_DUCKDB` waits for every output file of every merge, not one per merge | `subworkflows/local/BFD_MERGE.nf`, `BFD_FUNCTIONAL.nf` |

**Outside the pipeline (manual steps):**

| Item | How to run it |
|---|---|
| `seqdb/` protein and CDS copies. New or re-predicted genomes are not added automatically. | `Fungi_BFD_runs/bin/seqdb_list.py`, then `sbatch --array=... bin/seqdb_build.sh` (skips copies that exist; a re-predicted genome needs its old copy removed first) |
| bgzip of `scaffolds.fa` written before 2026-11-07 (done once for 23,450 files) | `Fungi_BFD_runs/bin/bgzip_scaffolds.sh` |
| `Fungi_BFD_runs/db/BFD.duckdb` is a link to `BFD_function/db/BFD.duckdb` | made by hand; keep it |
| `docs/BFD_SCHEMA.md` | `python3 nextflow/bin/schema_doc.py --db <db> --out docs/BFD_SCHEMA.md` |
| Tables written by the old code before a template change | rerun with `--force_merges` (section 7) |

**Not converted yet:** EggNOG, Phobius, TCDB, antiSMASH, tRNAscan-SE, CAZy CGC
and DeepTMHMM merges still use the Python parsers (no BFD-wide data yet).
TELOMERES and MERGE_SAMPLES stay on Python. Infernal and TargetP are in
`skip_merges` in `param_files/params_merge_only.yaml`.
