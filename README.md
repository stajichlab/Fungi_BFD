# Fungi5K

Comparative genomics analysis of **22k re-annotated fungal genomes** from NCBI (dataset frozen 2025-05).  
Analyzes gene structure, functional domains, protein orthology, and kingdom-wide phylogenetics across Fungi.  
See [project_README.md](project_README.md) for a full description of all analyses and outputs.

---

## Repository layout

```
samples.csv              ← master species table (primary input for all pipelines)
input/                   ← one protein FASTA per species
genomes/                 ← genome FASTA files
gff3/                    ← GFF3 annotation files
bigquery/                ← consolidated CSV.gz files loaded into DuckDB
results/                 ← per-tool output files
pipeline/
  nextflow/              ← Nextflow functional annotation pipeline  ← YOU ARE HERE
  function/              ← original SLURM array scripts (reference)
  db/                    ← DuckDB build scripts
scripts/                 ← Python/R helper scripts
sql/schema.sql           ← DuckDB table definitions
functionalDB/            ← function.duckdb  (full annotation database)
intronDB/                ← intron_db.duckdb (gene-structure database)
Phylogeny/               ← phylogenetic pipeline and outputs
```

---

## Primary input files

### `samples.csv`

The master species table. Every pipeline — SLURM array jobs and the Nextflow workflow — derives its sample list from this file. Columns:

| Column | Description |
|---|---|
| `ASMID` | NCBI assembly accession |
| `SPECIESIN` | Full species + strain label (used in legacy scripts) |
| `STRAIN` | Strain identifier; multiple entries separated by `;` |
| `BIOPROJECT` | NCBI BioProject accession |
| `NCBI_TAXONID` | NCBI taxonomy ID |
| `BUSCO_LINEAGE` | BUSCO lineage used for quality assessment |
| `PHYLUM` through `GENUS` | Taxonomic classification |
| `SPECIES` | Binomial species name (no strain) |
| `LOCUSTAG` | **8-character hex string** uniquely identifying each genome (e.g., `FF5840CF`) |

`LOCUSTAG` is the stable genome identifier used in all gene IDs, protein IDs, and database foreign keys.

### `input/` — protein FASTA files

One file per species, named `{SPECIES}_{STRAIN}.proteins.fa`:

- `SPECIES` = the `SPECIES` column, spaces replaced with `_`
- `STRAIN` = first `;`-delimited token from the `STRAIN` column, spaces replaced with `_`, single quotes stripped
- Example: `Aaosphaeria_arxii_CBS_175.79.proteins.fa`

Protein sequence identifiers follow the pattern `LOCUSTAG_NNNNN-TX.protein` (e.g., `FF5840CF_00001-T1.protein`).

---

## Nextflow functional annotation pipeline

`pipeline/nextflow/genome_functional.nf` is the primary way to run functional annotation. It replaces the original SLURM array scripts in `pipeline/function/` with a fully tracked, resumable Nextflow workflow.

### What it runs

Nine subworkflows execute in parallel, one per tool:

| Subworkflow | Tool | Output in `results/function/` | Merged to `bigquery/` |
|---|---|---|---|
| `PFAM` | hmmscan vs Pfam-A | `pfam_hmmscan/*.pfam.gz` | `pfam.csv.gz` |
| `CAZY` | dbcanlight (cazyme + sub) | `cazy/{locustag}/overview.tsv.gz` | `cazy.overview.csv.gz`, `cazy.cazymes_hmm.csv.gz` |
| `MEROPS` | blastp vs MEROPS lib | `merops/*.blasttab.gz` | `merops.csv.gz` |
| `SIGNALP` | SignalP 6 | `signalp/*.signalp.gff3.gz` | `signalp.signal_peptide.csv.gz` |
| `TMHMM` | TMHMM | `tmhmm/*.tmhmm_short.tsv.gz` | `tmhmm.csv.gz` |
| `TARGETP` | TargetP 2 | `targetP/*_summary.targetp2.gz` | `targetP.csv.gz` |
| `IDP` | AIUPred | `aiupred/*.aiupred.txt.gz` | `idp.csv.gz`, `idp_summary.csv.gz` |
| `WOLFPSORT` | runWolfPsortSummary | `wolfpsort/*.wolfpsort.results.txt.gz` | `wolfpsort.csv.gz` |
| `PREDGPI` | predgpi.py | `predgpi/*.predgpi.gff3.gz` | `predgpi.csv.gz` |

Each `results/function/<tool>/` file is gzip-compressed. Each `bigquery/<tool>.csv.gz` is the DuckDB-ready consolidated table.

### File structure

```
pipeline/nextflow/
  genome_functional.nf   ← main workflow
  nextflow.config        ← SLURM executor, resource labels, module loading
  conf/
    test.config          ← overrides for -profile test (stub data, minimal resources)
  bin/                   ← merge scripts (auto-added to PATH by Nextflow)
    merge_cazy.py
    merge_merops.py
    merge_signalp.py
    merge_tmhmm.py
    merge_targetp.py
    merge_wolfpsort.py
    merge_predgpi.py
  tests/
    data/
      test_samples.csv         ← 2-row CSV pointing at synthetic protein FASTAs
      input/                   ← synthetic protein FASTAs for stub testing
    validate_outputs.py        ← checks bigquery CSV.gz headers after a test run
  run_functional.sh      ← SLURM launcher (submit from project root)
  run_test.sh            ← lint + stub-run + validate (no real tools required)
  run_lint.sh            ← syntax check + py_compile for all scripts
```

### Configuration files

**`nextflow.config`** — the primary config. Contains:
- SLURM executor settings (queue size, poll interval)
- Per-tool resource allocations (`label` blocks: cpus, memory, time, queue)
- `beforeScript` per label — sources `/etc/profile.d/modules.sh` and loads the required environment module for each tool
- Paths to `workDir`, trace/report/timeline outputs
- Default `params` (all overridable on the command line)

**`conf/test.config`** — included automatically by `-profile test`. Overrides:
- `params.samples` → `tests/data/test_samples.csv`
- `params.input_dir` → `tests/data/input/`
- `params.outdir` / `params.bigquery` → `tests/output/`
- All resource labels → 1 CPU / 1 GB / 10 min (suitable for stub runs)

### Parameters

All parameters can be overridden on the command line with `--param value`.

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `samples.csv` (launch dir) | Master species table |
| `--input_dir` | `input/` (launch dir) | Directory of `*.proteins.fa` files |
| `--outdir` | `results/function/` | Per-species result files |
| `--bigquery` | `bigquery/` | Merged DuckDB-loadable CSV.gz files |
| `--scripts` | `scripts/` | Location of helper Python scripts |
| `--run_pfam` … `--run_predgpi` | `true` | Toggle individual subworkflows on/off |
| `--pfam_nodes` | `1` | SLURM nodes per hmmscan job; `1` = single-node threaded, `>1` = MPI |
| `--n_test` | `0` (all samples) | Limit to first N samples (for testing) |

### Pfam MPI mode

By default each `hmmscan` job runs on a single node using 16 threads (`--cpu 16`). For large sample sets, MPI mode distributes the Pfam-A HMM database across multiple nodes so each species finishes faster.

**Enable MPI (e.g. 4 nodes per species):**
```bash
sbatch nextflow/run_functional.sh --run_pfam true --pfam_nodes 4
```

Under the hood, when `pfam_nodes > 1` the pipeline:
1. Requests `-N <n> --ntasks-per-node=1` from SLURM for each Pfam job
2. Launches hmmscan as `srun --mpi=pmi2 -N <n> -n <n> hmmscan --mpi --cut_ga --cpu 16 …`

Each MPI rank searches a stripe of Pfam-A; `--cpu 16` still gives 16 search threads per rank. Memory is allocated per-job (48 GB total), so with 4 nodes each node uses roughly 12 GB.

**Recommendations:**
- `--pfam_nodes 1` (default) — adequate for ≤500 species or when queue wait time dominates
- `--pfam_nodes 4` — suitable for large runs on the `epyc` queue where multi-node jobs are available
- Do not exceed the number of HMM partitions that HMMER would naturally create; diminishing returns set in beyond ~8 nodes

### How to run

All commands are run **from the project root directory** (where `samples.csv` lives).

**Full production run (SLURM):**
```bash
sbatch pipeline/nextflow/run_functional.sh
```

**Run only specific tools:**
```bash
sbatch pipeline/nextflow/run_functional.sh --run_pfam true --run_cazy false --run_signalp false
```

**Test first 5 species (real tools, quick check):**
```bash
sbatch pipeline/nextflow/run_functional.sh --n_test 5
```

**Resume after a failure** (`-resume` is always passed by the launcher):
```bash
sbatch pipeline/nextflow/run_functional.sh
```

**Dry-run to check channel wiring without submitting jobs:**
```bash
module load nextflow
nextflow run pipeline/nextflow/genome_functional.nf \
    -c pipeline/nextflow/nextflow.config \
    -preview
```

**Stub-run test (no real tools — exercises the full DAG with placeholder outputs):**
```bash
sbatch pipeline/nextflow/run_test.sh
# or interactively:
bash pipeline/nextflow/run_test.sh
```

**Lint only:**
```bash
bash pipeline/nextflow/run_lint.sh
```

### Recommended first-run sequence

```bash
bash pipeline/nextflow/run_lint.sh                   # 1. syntax + py_compile
bash pipeline/nextflow/run_test.sh                   # 2. stub-run + validate
sbatch pipeline/nextflow/run_functional.sh --n_test 2  # 3. real run, 2 species
sbatch pipeline/nextflow/run_functional.sh           # 4. full production run
```

### Outputs and DuckDB loading

After the pipeline completes, `bigquery/` will contain one `.csv.gz` per tool. These are loaded into `functionalDB/function.duckdb` by:

```bash
sbatch pipeline/db/02_build_functional.sh
```

See `sql/schema.sql` for table definitions. To query interactively:
```bash
module load duckdb
duckdb functionalDB/function.duckdb
```

---

## Other pipelines

### MMseqs2 protein clustering (SLURM scripts)

```bash
sbatch pipeline/01_cluster_mmseqs.sh           # cluster at 30% ID / 70% coverage
sbatch pipeline/02_process_cluster_pairwise.sh  # per-cluster pairwise distances
sbatch pipeline/03_make_mmseq_orthogroups.sh    # assign orthogroup IDs
```

Output: `bigquery/mmseqs_orthogroup_clusters.csv.gz`

### Gene structure extraction

```bash
module load biopython
python scripts/build_genestats_bigquery.py -g gff3/ -d genomes/ -o bigquery/
```

Produces: `gene_info.csv.gz`, `gene_introns.csv.gz`, `gene_exons.csv.gz`, and related tables in `bigquery/`.

### Phylogenetics (run from `Phylogeny/` directory)

```bash
sbatch pipeline/01_phyling.sh
sbatch pipeline/02_phyling_filter_msa.sh
sbatch pipeline/03_phyling_make_tree.sh
sbatch pipeline/04_make_concatpartition.sh
```

Final tree: `Phylogeny/fungi_tree/final_tree.nw`

---

## Environment and modules

The cluster uses the UCR HPCC SLURM scheduler with an environment module system. Key modules:

| Purpose | Module |
|---|---|
| Pfam HMM scan | `hmmer/3.4`, `db-pfam` |
| CAZyme annotation | `dbcanlight` |
| Protease families | `db-merops/124`, `ncbi-blast/2.16.0+` |
| Signal peptides | `signalp/6-gpu` (requires A100 GPU) |
| TM helices | `tmhmm` |
| Subcellular targeting | `targetp` |
| Disorder prediction | `aiupred` (requires A100 GPU) |
| Subcellular localization | `wolfpsort` |
| GPI anchors | `predgpi` |
| Protein clustering | `mmseqs2` |
| Database queries | `duckdb` |
| Phylogenetics | `phyling`, `phykit`, `fasttree`, `modeltest-ng` |
| Python scripts | `biopython` |
| Workflow engine | `nextflow` (≥23.04) |

The Nextflow config sources `/etc/profile.d/modules.sh` in a `beforeScript` block for each process label, so modules load correctly on worker nodes.

---

## Contact

Jason Stajich — jasonst@ucr.edu — Stajich Lab
