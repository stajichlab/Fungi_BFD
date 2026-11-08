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
WHERE protein_id LIKE 'F77FFF1D_%';
```

**Note on the DuckDB file.** From this branch on, `db/BFD.duckdb` holds views
that read `tables/*.parquet` by absolute path (274 KB). It works only where
`tables/` is at the same path. The `db/BFD.duckdb` built on 2026-11-06 is the
old style, with a full copy of each table (394 GB). Views are a little slower
than the indexed copies: 8.5 s against 1.9 s to count one genome's proteins,
and 15.9 s against 13.3 s to fetch one protein's Pfam hits.

## 4. Tables

Row counts come from the 2026-11-06 build. The full column lists are in
`docs/BFD_SCHEMA.md`, which is generated from the database by
`nextflow/bin/schema_doc.py`. That file is on the `merge-csv-list` branch (PR #54).

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

- **SwissProt** is not in the database yet. `tables/swissprot.parquet` (146 GB)
  is the old format. A smaller design (hits per bucket with a subject table,
  plus best-hit and summary tables) is in PR #54 and will be built by the next
  merge run.
- **Infernal (Rfam)**: the Infernal steps have not run for this set.
  `tables/infernal_rfam.parquet` is from the old parser, which read the wrong
  columns. Do not use it.
- **TargetP** has results for 127 genomes only. A move to DeepLoc is being considered.
- **DeepTMHMM** is not part of this release. Use `tmhmm`.
- **TF inventory** (`tables/tf_inventory.parquet`) exists, but it is not yet in
  the DuckDB file. Query the Parquet file directly.
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
   `nextflow/bin/duckdb_merge.py`).
3. `BUILD_DUCKDB` runs `nextflow/bin/build_BFD_duckDB.sh` to make `db/BFD.duckdb`.
   It waits until all Parquet copies to `tables/` are complete.
4. A merge whose inputs have not changed is skipped (`merge_if_stale`). A new
   table that has fewer rows than the old one is not published unless
   `allow_shrink` is set.
