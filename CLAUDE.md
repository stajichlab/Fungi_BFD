# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Fungi5K is a comparative genomics analysis project for 5,813 re-annotated fungal genomes from NCBI (dataset frozen 2024-08-12). The project analyzes gene structure (introns, exons), functional domains (Pfam, CAZY, MEROPS), protein orthology (MMseqs2 clustering), and phylogenetics across kingdom Fungi.

## Environment

All compute jobs run on a SLURM HPC cluster using the module system. Scripts are written as SLURM batch scripts (`sbatch`) and should be submitted via `sbatch <script>` rather than run directly. Key modules used: `mmseqs2`, `duckdb`, `phyling`, `phykit`, `hmmer`, `AAFTF`, `biopython`, `fasttree`, `modeltest-ng`.

Python scripts require `biopython`, `pyfaidx`, and `pybedtools` (loaded via `module load biopython`).

## Key Data Files

- `samples.csv` — master species metadata table (ASMID, SPECIESIN, STRAIN, BIOPROJECT, NCBI_TAXONID, BUSCO_LINEAGE, PHYLUM, SUBPHYLUM, CLASS, SUBCLASS, ORDER, FAMILY, GENUS, SPECIES, LOCUSTAG)
- `sampleset.txt` — list of protein input filenames used by array jobs
- `genomes/` — genome FASTA files (`*.scaffolds.fa`)
- `gff3/` — GFF3 annotation files
- `input/` — protein FASTA files (`*.proteins.fa`), one per species
- `bigquery/` — intermediate CSV.gz files used to build DuckDB databases
- `fungi_introns/all_fungi_introns.fa.gz` — concatenated FASTA of all intron sequences with a multi-volume BLAST database (`*.nhr`, `*.nin`, `*.nsq`) and USEARCH UDB (`*.udb`)

## Identifier Conventions

- **LOCUSTAG** — 8-character hex string (e.g., `FF5840CF`) uniquely identifying each species/genome; first column of `samples.csv`
- **Gene IDs** — `LOCUSTAG_NNNNN` (e.g., `FF5840CF_00001`)
- **Transcript IDs** — `LOCUSTAG_NNNNN-TX` (e.g., `FF5840CF_00001-T1`)
- **Protein IDs** — `LOCUSTAG_NNNNN-TX.protein`
- In DuckDB tables, `locustag` is derived by splitting on `_` and taking the first element, or via `substring(id,1,8)`

## Pipeline Architecture

### 1. Genome Statistics (`pipeline/db/00_genome_stats.sh`)
Runs AAFTF assess on all genomes in `genomes/`, then calls `scripts/collect_asm_stats.py` to produce `bigquery/asm_stats.csv.gz`. Chromosome-level stats are collected by `scripts/collect_chrom_info.py` → `bigquery/chrom_info.csv.gz`.

### 2. Gene Statistics from GFF3 (`scripts/build_genestats_bigquery.py`)
Parses GFF3 + genome FASTA to produce CSV tables in `bigquery/`:
- `gene_info.csv.gz`, `gene_transcripts.csv.gz`, `gene_exons.csv.gz`, `gene_CDS.csv.gz`
- `gene_introns.csv.gz` (includes splice sites, GC content, codon position/frame, sequence)
- `gene_proteins.csv.gz`, `gene_trnas.csv.gz`

Usage: `python scripts/build_genestats_bigquery.py -g gff3/ -d genomes/ -o bigquery/`

Amino acid and codon usage frequencies are computed separately:
```bash
python scripts/calculate_AA_freq.py    → bigquery/aa_freq.csv.gz
python scripts/calculate_codon_freq.py → bigquery/codon_freq.csv.gz
```

### 3. Functional Annotation (`pipeline/function/`)
SLURM array jobs run tools on `input/*.proteins.fa`:
- `01_pfam.sh` — hmmscan against Pfam-A (produces `results/function/pfam_hmmscan/*.pfam.gz`); `10_process_pfam_TSV.sh` / `10_convert_TOR_pfam_TSV.sh` handle TSV-format Pfam output from alternative runs
- `02_cazy.sh` — run_dbcan for CAZyme annotation
- `05_merops.sh` — MEROPS protease family BLAST
- `06_signalp.sh` — signal peptide prediction
- `07_tmhmm.sh`, `07_targetp.sh`, `07_IDP_predict.sh`, `07_wolfpsort.sh`, `07_predGPI.sh` — TM helices, subcellular targeting, disorder, WoLF PSORT, GPI anchors
- `09_summarize_function_bigquery.sh` — calls `scripts/prep_for_bigquery_load.py` to consolidate results into `bigquery/` CSV.gz files

FunGuild ecological guild annotations are fetched via `scripts/download_funguild.py` and converted to `bigquery/species_funguild.csv.gz` by `scripts/build_funguild_bigquery.py`.

### 4. MMseqs2 Protein Clustering (`pipeline/`)
- `01_cluster_mmseqs.sh` — clusters all proteins at 30% identity, 70% coverage → `results/Fungi5K_cluster.tsv.gz`
- `02_process_cluster_pairwise.sh` — runs `scripts/mmseqs2pairwise.py` to generate per-cluster pairwise distances → `bigquery/gene_pairwise_distances.csv.gz`
- `03_make_mmseq_orthogroups.sh` — runs `scripts/mmseqs2bigqueryload.py` and `scripts/mmseqs2orthogroups.py` → `bigquery/mmseqs_orthogroup_clusters.csv.gz`

### 5. DuckDB Database Loading
Two databases are maintained:

**`intronDB/intron_db.duckdb`** (`pipeline/db/01_build_intronDB.sh`): gene structure tables (`gene_info`, `gene_introns`, `gene_transcripts`, `gene_exons`, `gene_proteins`) plus `gene_pairwise_distances` from inter-species intron BLAST comparisons.

**`functionalDB/function.duckdb`** (`pipeline/db/02_build_functional.sh`): full functional annotation database including all gene structure tables plus Pfam, CAZy, MEROPS, SignalP, TMHMM, TargetP, IDP, WoLF PSORT, Prosite, FunGuild, and MMseqs2 orthogroups.

Post-load processing (`pipeline/db/03_process_db.sh`) runs `scripts/Rscripts/mmseq_cluster_profile.R` and `scripts/mmseqs_genedump_to_fasta.py` to extract unannotated orthogroup proteins.

Both databases are built by loading `bigquery/*.csv.gz` directly using DuckDB's `read_csv_auto()`. See `sql/schema.sql` for all table definitions and indexes.

### 6. Phylogenetics (`Phylogeny/pipeline/`)
Run from `Phylogeny/` directory:
1. `01_phyling.sh` — align with PHYling using `fungi_odb12` markers
2. `02_phyling_filter_msa.sh` — filter MSAs requiring ≥80% taxon occupancy
3. `03_phyling_make_tree.sh` — build tree with FastTree (`-M ft`)
4. `04_make_concatpartition.sh` — concatenate alignments with phykit, run modeltest-ng for partitioned analysis

## Database Schema

See `sql/schema.sql` for the complete DuckDB table definitions. Key tables:
- `species` — from `samples.csv`, links LOCUSTAG to taxonomy
- `gene_introns` — intron features with splice sites, GC content, phase
- `pfam` — Pfam domain hits per protein (also `pfam_UoT` for an alternative Pfam run)
- `mmseqs_orthogroup_clusters` — orthogroup membership
- `mmseqs_orthogroup_cluster_count` — precomputed orthogroup sizes
- `aa_frequency`, `codon_frequency` — per-species amino acid and codon usage
- `funguild` — ecological guild annotations (trophicMode, guild, growthForm)
- `prosite` — ProSite pattern scan results
- `wolfpsort`, `targetp` — subcellular localization predictions

## Running Scripts

```bash
# Submit a pipeline step to SLURM
sbatch pipeline/01_cluster_mmseqs.sh

# Run gene stats extraction (requires biopython module)
module load biopython
python scripts/build_genestats_bigquery.py -g gff3/ -d genomes/ -o bigquery/

# Query the functional database
module load duckdb
duckdb functionalDB/function.duckdb

# Query the intron database
duckdb intronDB/intron_db.duckdb

# Extract proteins carrying a specific Pfam domain
python scripts/dump_proteins_by_domain.py -d functionalDB/function.duckdb -p PF00096
```

## Interactive Dashboard

`functional_dashboard.py` is a Dash/Plotly web app for exploring `functionalDB/function.duckdb`. It requires `duckdb`, `pandas`, `plotly`, and `dash`. Run with:
```bash
python functional_dashboard.py
```

## R Analysis

R scripts in `scripts/Rscripts/` perform exploratory analysis and visualization. Run with `Rscript`. Results and plots go to `results/` and `plots/`.

| Script | Purpose |
|---|---|
| `intron_explore.R` | Intron length distributions and splice site usage |
| `intron_size_distro.R` | Intron size distributions across taxa |
| `explore_pfam.R` | Pfam domain frequency across species |
| `summarize_pfam.R` | Pfam domain summary tables |
| `cazy_merops_ratio_counts.R` | CAZy and MEROPS family counts |
| `function_domain_IDR_explore.R` | Correlation of IDR regions with domains |
| `genome_size_explore.R` | Genome size vs gene count trends |
| `explore_stats_genomeproperties.R` | Assembly quality statistics |
| `exoncount_explore.R` | Exon count distributions |
| `gene_density_score.R` | Gene density analysis |
| `mmseq_cluster_profile.R` | Orthogroup size and conservation profiles |

### Intron length ~ Pfam domain correlation (`scripts/Rscripts/intron_pfam_correlation.R`)

Tests whether intron length is correlated with presence/absence or copy number of each Pfam domain using two complementary analyses:

- **Species-level**: Spearman correlation and Wilcoxon rank-sum test between species mean intron length and per-species domain copy number.
- **Gene-level**: Wilcoxon rank-sum test comparing per-gene mean intron length between proteins that carry vs. lack each domain.

Both analyses apply BH-FDR correction and report rank-biserial *r* as effect size.

**Stratification** — the `--stratify` flag repeats the analysis independently within each phylum or subphylum:

```bash
Rscript scripts/Rscripts/intron_pfam_correlation.R                          # all taxa
Rscript scripts/Rscripts/intron_pfam_correlation.R --stratify phylum        # per phylum
Rscript scripts/Rscripts/intron_pfam_correlation.R --stratify subphylum     # per subphylum
Rscript scripts/Rscripts/intron_pfam_correlation.R --stratify subphylum --taxon Pezizomycotina
```

Key options: `--evalue` (Pfam domain_i_evalue cutoff, default `1e-5`), `--min-species` / `--min-species-strat`, `--min-genes` / `--min-genes-strat`, `--fdr`, `--top-n`.

Outputs land in `results/pfam_intron_corr/<stratum>/` (TSV tables) and `plots/pfam_intron_corr/<stratum>/` (volcano plots, violin plots). A cross-stratum overview is written to `results/pfam_intron_corr/stratum_summary.tsv`.
