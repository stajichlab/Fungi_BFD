# Fungi5K

Comparative genomics analysis of **~22k re-annotated fungal genomes** from NCBI (dataset frozen 2025-05).
Analyzes gene structure, functional domains, protein orthology, and kingdom-wide phylogenetics across Fungi.

---

## Repository layout

```
samples.csv              ← master species table (primary input for all pipelines)
input/
  pep/                   ← protein FASTAs  ({Species_Strain}.proteins.fa)
  cds/                   ← CDS transcript FASTAs  ({Species_Strain}.cds-transcripts.fa)
  gff3/                  ← GFF3 annotation files  ({Species_Strain}.gff3)
  dna/                   ← genome FASTA files  ({Species_Strain}.scaffolds.fa)
  trna/                  ← tRNA GFF3 files  ({Species_Strain}.trna.gff3)
genome_annotation/       ← funannotate predict output  ({Species_Strain}/predict_results/)
tables/                  ← consolidated CSV.gz files loaded into DuckDB (Nextflow output)
bigquery/                ← consolidated CSV.gz files (legacy SLURM script output)
results/                 ← per-tool output files
nextflow/                ← Nextflow pipelines
scripts/                 ← Python/R helper scripts
sql/schema.sql           ← DuckDB table definitions
functionalDB/            ← function.duckdb  (full annotation database)
intronDB/                ← intron_db.duckdb (gene-structure database)
Phylogeny/               ← phylogenetic pipeline and outputs
```

The `input/` subdirectories are populated by the `SETUP_INPUT` step of `BFD.nf`,
which creates per-species symlinks into `genome_annotation/{Species_Strain}/predict_results/`.
Run `funannotate.nf` first to produce those predict results, then `BFD.nf`.

---

## Primary input file: `samples.csv`

Every pipeline derives its sample list from this file. Key columns:

| Column | Description |
|---|---|
| `ASMID` | NCBI assembly accession |
| `NCBI_TAXONID` | NCBI taxonomy ID |
| `BUSCO_LINEAGE` | BUSCO lineage for quality assessment |
| `PHYLUM` … `GENUS` | Taxonomic classification |
| `SPECIES` | Binomial species name |
| `STRAIN` | Strain identifier (`;`-delimited for multiple) |
| `LOCUSTAG` | **8-character hex string** — stable genome identifier used in all gene IDs and DB foreign keys |

---

## Nextflow pipelines

Four workflows live under `nextflow/`, all sharing one config with per-pipeline profiles.

| Pipeline | Profile | Launcher | Detail |
|---|---|---|---|
| `funannotate.nf` | `funannotate` | `run_funannotate.sh` | [README_funannotate.md](README_funannotate.md) |
| `BFD.nf` | `BFD` | `run_functional.sh` | [README_BFD.md](README_BFD.md) |
| `genome_seqstats.nf` | `BFD` | `run_seqstats.sh` | per-species sequence statistics |
| `interproscan6.nf` | `interproscan6` | — | InterProScan 6 XML for annotate_misc/ |

All workflows accept `--taxon RANK:VALUE` (e.g. `--taxon PHYLUM:Ascomycota`) and
`--n_test N` to restrict to the first N samples.

---

## Recommended run order

```bash
# 1. Syntax check (no SLURM, no tools)
bash nextflow/run_lint.sh

# 2. Stub-run: validate DAG and output structure without real tools
bash nextflow/run_test.sh

# 3. Genome annotation (funannotate) → produces genome_annotation/
sbatch nextflow/run_funannotate.sh --n_test 2   # pilot
sbatch nextflow/run_funannotate.sh               # full run

# 4. Functional annotation + input setup (BFD.nf)
sbatch nextflow/run_functional.sh --n_test 2    # pilot
sbatch nextflow/run_functional.sh               # full run

# 5. Per-species sequence statistics
sbatch nextflow/run_seqstats.sh --n_test 2      # pilot
sbatch nextflow/run_seqstats.sh                  # full run

# 6. Load results into DuckDB
sbatch pipeline/db/02_build_functional.sh
```

All Nextflow commands support `-resume` to restart from the last successful checkpoint.

---

## Loading results into DuckDB

After Nextflow pipelines complete, `tables/` contains one `.csv.gz` per table.

```bash
sbatch pipeline/db/02_build_functional.sh

module load duckdb
duckdb functionalDB/function.duckdb
```

See `sql/schema.sql` for all table definitions and indexes.

---

## Other analyses

**MMseqs2 protein clustering** (SLURM scripts, run from project root):

```bash
sbatch pipeline/01_cluster_mmseqs.sh           # cluster at 30% ID / 70% coverage
sbatch pipeline/02_process_cluster_pairwise.sh  # per-cluster pairwise distances
sbatch pipeline/03_make_mmseq_orthogroups.sh    # assign orthogroup IDs
```

**Phylogenetics** (run from `Phylogeny/`):

```bash
sbatch pipeline/01_phyling.sh
sbatch pipeline/02_phyling_filter_msa.sh
sbatch pipeline/03_phyling_make_tree.sh
sbatch pipeline/04_make_concatpartition.sh
```

---

## Contact

Jason Stajich — jasonst@ucr.edu — Stajich Lab
