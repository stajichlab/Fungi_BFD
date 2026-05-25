# Fungi5K

Comparative genomics analysis of **22k re-annotated fungal genomes** from NCBI (dataset frozen 2025-05).  
Analyzes gene structure, functional domains, protein orthology, and kingdom-wide phylogenetics across Fungi.  
See [project_README.md](project_README.md) for a full description of all analyses and outputs.

---

## Repository layout

```
samples.csv              ← master species table (primary input for all pipelines)
pep/                     ← one protein FASTA per species
cds/                     ← CDS transcript FASTAs
genome/                  ← genome FASTA files
gff3/                    ← GFF3 annotation files
tables/                  ← consolidated CSV.gz files loaded into DuckDB (Nextflow output)
bigquery/                ← consolidated CSV.gz files (legacy SLURM script output)
results/                 ← per-tool output files
nextflow/                ← Nextflow pipelines  ← YOU ARE HERE
scripts/                 ← Python/R helper scripts
sql/schema.sql           ← DuckDB table definitions
functionalDB/            ← function.duckdb  (full annotation database)
intronDB/                ← intron_db.duckdb (gene-structure database)
Phylogeny/               ← phylogenetic pipeline and outputs
```

---

## Primary input files

### `samples.csv`

The master species table. Every pipeline derives its sample list from this file. Columns:

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

---

## Nextflow pipelines

Four Nextflow workflows live under `nextflow/`, all driven by a single unified config with per-pipeline profiles.

```
nextflow/
  BFD.nf    ← BFD functional annotation (Pfam, CAZy, MEROPS, SignalP, …)
  genome_seqstats.nf      ← BFD sequence statistics (AA freq, codon freq, gene stats, chrom info)
  funannotate.nf          ← genome prediction + annotation (funannotate/PASA/Augustus)
  interproscan6.nf        ← InterProScan 6 XML generation for annotate_misc/
  nextflow.config         ← unified config: shared executor + profile stanzas
  conf/
    profile_BFD.config          ← params + process labels for BFD pipelines
    profile_funannotate.config  ← params + withName blocks for funannotate pipeline
    profile_interproscan6.config← params + withName block for interproscan6
    test.config                 ← overrides for -profile test (stub data, minimal resources)
  bin/                    ← merge scripts (auto-added to PATH by Nextflow)
  run_functional.sh       ← launcher: sbatch this from project root
  run_seqstats.sh         ← launcher: sbatch this from project root
  run_funannotate.sh      ← launcher: sbatch this from project root
  run_lint.sh             ← syntax check for all scripts
  run_test.sh             ← stub-run + validate (no real tools required)
  tests/
    data/                 ← synthetic inputs for stub testing
    validate_outputs.py   ← checks CSV.gz headers after a test run
```

Each workflow is activated by its `-profile` flag; shared settings (SLURM executor, error strategy, queue defaults) are declared once at the top level of `nextflow.config`.

---

### Workflow 1: Functional annotation (`BFD.nf`)

Runs nine tools in parallel across all species in `samples.csv`. Input: `pep/*.proteins.fa`.

| Subworkflow | Tool | Output in `results/function/` | Merged to `tables/` |
|---|---|---|---|
| `PFAM` | hmmscan vs Pfam-A | `pfam_hmmscan/*.pfam.gz` | `pfam.csv.gz` |
| `CAZY` | dbcanlight | `cazy/{locustag}/overview.tsv.gz` | `cazy.overview.csv.gz`, `cazy.cazymes_hmm.csv.gz` |
| `MEROPS` | blastp vs MEROPS | `merops/*.blasttab.gz` | `merops.csv.gz` |
| `SIGNALP` | SignalP 6 | `signalp/*.signalp.gff3.gz` | `signalp.signal_peptide.csv.gz` |
| `TMHMM` | TMHMM | `tmhmm/*.tmhmm_short.tsv.gz` | `tmhmm.csv.gz` |
| `TARGETP` | TargetP 2 | `targetP/*_summary.targetp2.gz` | `targetP.csv.gz` |
| `IDP` | AIUPred | `aiupred/*.aiupred.txt.gz` | `idp.csv.gz`, `idp_summary.csv.gz` |
| `WOLFPSORT` | runWolfPsortSummary | `wolfpsort/*.wolfpsort.results.txt.gz` | `wolfpsort.csv.gz` |
| `PREDGPI` | predgpi.py | `predgpi/*.predgpi.gff3.gz` | `predgpi.csv.gz` |

**Run:**
```bash
sbatch nextflow/run_functional.sh                          # full run
sbatch nextflow/run_functional.sh --n_test 5               # first 5 species
sbatch nextflow/run_functional.sh --run_pfam true --run_cazy false  # selective
```

**Key parameters** (all overridable with `--param value`):

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `samples.csv` | Master species table |
| `--pep_dir` | `pep/` | Directory of `*.proteins.fa` files |
| `--outdir` | `results/function/` | Per-species result files |
| `--tables` | `tables/` | Merged DuckDB-loadable CSV.gz files |
| `--run_pfam` … `--run_predgpi` | `true` | Toggle individual subworkflows |
| `--pfam_nodes` | `1` | SLURM nodes per hmmscan job (`>1` = MPI mode) |
| `--n_test` | `0` (all) | Limit to first N samples |

---

### Workflow 2: Sequence statistics (`genome_seqstats.nf`)

Computes per-species amino acid frequencies, codon frequencies, gene structure tables (7 CSVs from GFF3), and chromosome/scaffold stats. Input: `pep/`, `cds/`, `gff3/`, `genome/`.

| Subworkflow | Script | Output in `results/seqstats/` | Merged to `tables/` |
|---|---|---|---|
| `AA_FREQ` | `calculate_AA_freq.py` | `aa_freq/{locustag}.aa_freq.csv` | `aa_freq.csv.gz` |
| `CODON_FREQ` | `calculate_codon_freq.py` | `codon_freq/{locustag}.codon_freq.csv` | `codon_freq.csv.gz` |
| `GENE_STATS` | `build_genestats_bigquery.py` | `gene_stats/{locustag}/` | `gene_info.csv.gz`, `gene_transcripts.csv.gz`, `gene_exons.csv.gz`, `gene_CDS.csv.gz`, `gene_introns.csv.gz`, `gene_trnas.csv.gz`, `gene_proteins.csv.gz` |
| `CHROM_INFO` | `collect_chrom_info.py` | `chrom_info/{locustag}.chrom_info.csv` | `chrom_info.csv.gz` |

**Run:**
```bash
sbatch nextflow/run_seqstats.sh                                    # full run
sbatch nextflow/run_seqstats.sh --n_test 5                         # first 5 species
sbatch nextflow/run_seqstats.sh --run_gene_stats false             # skip gene stats
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--run_aa_freq` | `true` | Toggle amino acid frequency subworkflow |
| `--run_codon_freq` | `true` | Toggle codon frequency subworkflow |
| `--run_gene_stats` | `true` | Toggle gene structure table subworkflow |
| `--run_chrom_info` | `true` | Toggle chromosome info subworkflow |

---

### Workflow 3: Genome prediction and annotation (`funannotate.nf`)

Runs the full funannotate predict/annotate pipeline: genome cleaning (NCBI FCS-GX), optional repeat masking, RNA-seq download (SRA), funannotate train → predict → annotate. Designed for 1KFG-style annotation projects.

**Run:**
```bash
sbatch nextflow/run_funannotate.sh                                        # default (predict only)
sbatch nextflow/run_funannotate.sh --run_annotate true                    # predict + annotate
sbatch nextflow/run_funannotate.sh --n_test 2 --only_clean true           # genome clean only
sbatch nextflow/run_funannotate.sh --run_repeatmasker false               # skip masking
sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Ascomycota              # Ascomycota only
sbatch nextflow/run_funannotate.sh --taxon CLASS:Sordariomycetes          # single class
sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Basidiomycota --n_test 5  # test subset
```

**Key parameters:**

| Parameter | Default | Description |
|---|---|---|
| `--target` | `annotation_22k/` | Output directory for per-species funannotate folders |
| `--source` | (1KFG NCBI ASM path) | Directory containing `{ASMID}/{ASMID}_genomic.fna.gz` |
| `--taxon` | `""` (all) | Restrict to samples matching `RANK:VALUE` (e.g. `PHYLUM:Ascomycota`, `CLASS:Sordariomycetes`). Rank must be an uppercase column name in `samples.csv` (`PHYLUM`, `SUBPHYLUM`, `CLASS`, `ORDER`, `FAMILY`, `GENUS`). |
| `--n_test` | `0` (all) | Limit to first N samples after filtering |
| `--suppress` | `""` | Path to file of ASMIDs to skip (one per line) |
| `--run_sra_fetch` | `true` | Download RNA-seq from NCBI SRA before training |
| `--run_repeatmasker` | `true` | Soft-mask genome with RepeatMasker |
| `--run_repeatmodeler` | `false` | Build de-novo repeat library with RepeatModeler |
| `--run_annotate` | `false` | Run funannotate annotate after predict |
| `--run_update` | `false` | Run funannotate update (PASA re-training) before annotate; requires `--run_sra_fetch true` |
| `--run_antismash` | `false` | Run antiSMASH before annotate |
| `--run_interpro` | `false` | Run InterProScan 5 before annotate |
| `--run_signalp` | `false` | Run SignalP 6 (GPU) before annotate |
| `--only_clean` | `false` | Stop after genome cleaning |

The pipeline uses `storeDir` to cache genome-cleaning, repeat library, SRA download, and RNA-seq preparation results — rerunning skips these expensive steps automatically.

---

### Workflow 4: InterProScan 6 (`interproscan6.nf`)

Runs the `ebi-pf-team/interproscan6` sub-pipeline (Nextflow-in-Nextflow) to produce `iprscan.xml` for each species that has predict results but no existing XML. Requires `nextflow pull ebi-pf-team/interproscan6` on the head node first.

**Run:**
```bash
nextflow run nextflow/interproscan6.nf \
    -c nextflow/nextflow.config \
    -profile interproscan6 \
    -resume
```

---

### Stub-run and testing

```bash
bash nextflow/run_lint.sh                                    # syntax + py_compile
bash nextflow/run_test.sh                                    # stub-run + validate
```

The `stub` profile limits to 2 samples and uses placeholder outputs — no real tools are invoked. This validates the full DAG and channel wiring without submitting real SLURM jobs.

### Recommended first-run sequence (BFD pipelines)

```bash
bash nextflow/run_lint.sh                                    # 1. syntax check
bash nextflow/run_test.sh                                    # 2. stub-run + validate
sbatch nextflow/run_functional.sh --n_test 2                 # 3. real run, 2 species
sbatch nextflow/run_functional.sh                            # 4. full functional run
sbatch nextflow/run_seqstats.sh                              # 5. sequence statistics
```

### Loading results into DuckDB

After the Nextflow pipelines complete, `tables/` contains one `.csv.gz` per table. Load into `functionalDB/function.duckdb`:

```bash
sbatch pipeline/db/02_build_functional.sh
```

See `sql/schema.sql` for table definitions. Query interactively:
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

### Gene structure extraction (standalone)

```bash
module load biopython
python scripts/build_genestats_bigquery.py -g gff3/ -d genome/ -o bigquery/
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
| Gene/sequence statistics | `biopython` |
| Genome prediction | `funannotate`, `augustus`, `RepeatMasker`, `RepeatModeler` |
| Genome cleaning | `AAFTF`, `taxonkit` |
| Protein clustering | `mmseqs2` |
| Database queries | `duckdb` |
| Phylogenetics | `phyling`, `phykit`, `fasttree`, `modeltest-ng` |
| Workflow engine | `nextflow` (≥25.10) |

The BFD Nextflow config sources `/etc/profile.d/modules.sh` in a `beforeScript` block for each process label, so modules load correctly on worker nodes without manual setup.

---

## Contact

Jason Stajich — jasonst@ucr.edu — Stajich Lab
