# BFD Functional Annotation Pipeline (`BFD.nf`)

Nextflow DSL2 pipeline that runs nine functional annotation tools across all species in
`samples.csv`, then consolidates each tool's per-species results into a single merged
CSV.gz in `tables/` ready for loading into DuckDB.

Run after `funannotate.nf` has produced `genome_annotation/` prediction results.

---

## Quick start

```bash
# Full run: input setup + all annotation tools + merge
sbatch nextflow/run_functional.sh

# Add new species to an existing dataset (run new ones, merge ALL into tables/)
sbatch nextflow/run_functional.sh --taxon PHYLUM:Basidiomycota

# Run specific tools only (skip others)
sbatch nextflow/run_functional.sh \
    --run_merops false --run_signalp false --run_tmhmm false \
    --run_targetp false --run_idp false --run_wolfpsort false --run_predgpi false

# Merge all existing result files without re-running any tools
sbatch nextflow/run_functional.sh \
    --run_pfam false --run_cazy false --run_merops false \
    --run_signalp false --run_tmhmm false --run_targetp false \
    --run_idp false --run_wolfpsort false --run_predgpi false \
    --run_aa_freq false --run_codon_freq false \
    --run_intergenic false --run_gene_stats false \
    --merge_all true --skip_merge false

# Run new species, defer the merge step
sbatch nextflow/run_functional.sh --skip_merge true

# Restrict to a taxonomic subset
sbatch nextflow/run_functional.sh --taxon PHYLUM:Ascomycota
sbatch nextflow/run_functional.sh --taxon CLASS:Dothideomycetes --n_test 5

# Test: first 5 species only
sbatch nextflow/run_functional.sh --n_test 5
```

---

## Pipeline stages

```
samples.csv
    │
    ▼
SETUP_INPUT              Create symlinks in input/ pointing to funannotate predict_results/
    │                    Skips species whose symlinks already exist.
    │                    (--run_setup false to skip entirely if input/ is pre-populated)
    ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Per-species annotation (all tools run in parallel per species)         │
 │                                                                          │
 │  RUN_PFAM        hmmsearch --cut_ga vs Pfam-A → <name>.pfam.gz          │
 │  RUN_CAZY        dbcanlight cazyme + substrate → <name>.overview.tsv.gz │
 │  RUN_MEROPS      blastp vs MEROPS scan.lib → <name>.blasttab.gz         │
 │  RUN_SIGNALP     SignalP 6 (GPU, fast mode) → <name>.signalp.gff3.gz    │
 │  RUN_TMHMM       TMHMM short format → <name>.tmhmm_short.tsv.gz         │
 │  RUN_TARGETP     TargetP 2 non-plant → <name>_summary.targetp2.gz       │
 │  RUN_IDP         AIUPred disorder prediction → <name>.idp.csv.gz        │
 │  RUN_WOLFPSORT   runWolfPsortSummary fungi → <name>.wolfpsort.results.gz │
 │  RUN_PREDGPI     predgpi.py GFF3 mode → <name>.predgpi.gff3.gz          │
 │                                                                          │
 │  All per-species processes use storeDir: automatically skipped if all    │
 │  output files already exist on disk (no -resume required).               │
 └──────────────────────────────────────────────────────────────────────────┘
    │
    ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  Whole-dataset bulk jobs (one SLURM job each, reads entire input/ dir)  │
 │                                                                          │
 │  CALC_AA_FREQ       calculate_AA_freq.py → tables/aa_freq.csv.gz        │
 │  CALC_CODON_FREQ    calculate_codon_freq.py → tables/codon_freq.csv.gz  │
 │  CALC_INTERGENIC    calculate_intergenic.py → tables/gene_intergenic…   │
 │  CALC_GENE_STATS    build_genestats_table.py → tables/gene_*.csv.gz     │
 └──────────────────────────────────────────────────────────────────────────┘
    │
    ▼
 ┌──────────────────────────────────────────────────────────────────────────┐
 │  MERGE steps (one per tool, controlled by --merge_all and --skip_merge) │
 │                                                                          │
 │  merge_all=true (default): glob ALL *.gz files in results/function/     │
 │    <tool>/ regardless of which species were run in this session. A sync  │
 │    barrier ensures newly-run files are on disk before the glob fires.   │
 │    MERGE is silently skipped if no files are found for a tool.          │
 │                                                                          │
 │  merge_all=false: merge only the files produced in the current run.     │
 │                                                                          │
 │  skip_merge=true: skip all MERGE steps unconditionally.                 │
 └──────────────────────────────────────────────────────────────────────────┘
```

---

## Output structure

```
input/                          ← populated by SETUP_INPUT (symlinks into genome_annotation/)
  pep/    <name>.proteins.fa
  cds/    <name>.cds-transcripts.fa
  gff3/   <name>.gff3
  dna/    <name>.scaffolds.fa
  trna/   <name>.trna.gff3

results/function/               ← per-species per-tool outputs (storeDir)
  pfam_hmmscan/
    <name>.pfam.gz              domtbl format
    <name>.tblout.gz
  cazy/<name>/
    <name>.overview.tsv.gz
    <name>.cazymes.tsv.gz
    <name>.substrates.tsv.gz
  merops/
    <name>.blasttab.gz
  signalp/
    <name>.signalp.gff3.gz
    <name>.signalp.results.txt.gz
  tmhmm/
    <name>.tmhmm_short.tsv.gz
    <name>.tmhmm_results.tsv.gz
  targetP/
    <name>_summary.targetp2.gz
  aiupred/
    <name>.aiupred.txt.gz
    <name>.idp.csv.gz
    <name>.idp_summary.csv.gz
  wolfpsort/
    <name>.wolfpsort.results.txt.gz
  predgpi/
    <name>.predgpi.gff3.gz

tables/                         ← merged DuckDB-loadable CSV.gz (MERGE output)
  pfam.csv.gz
  cazy.overview.csv.gz
  cazy.cazymes_hmm.csv.gz
  merops.csv.gz
  signalp.signal_peptide.csv.gz
  tmhmm.csv.gz
  targetP.csv.gz
  idp.csv.gz
  idp_summary.csv.gz
  wolfpsort.csv.gz
  predgpi.csv.gz
  aa_freq.csv.gz
  codon_freq.csv.gz
  gene_intergenic_distances.csv.gz
  gene_info.csv.gz
  gene_transcripts.csv.gz
  gene_exons.csv.gz
  gene_CDS.csv.gz
  gene_introns.csv.gz
  gene_trnas.csv.gz
  gene_proteins.csv.gz

logs/nextflow/
  BFD_trace.txt
  BFD_report.html
  BFD_timeline.html
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `samples.csv` | Master species/genome table |
| `--genome_annotation` | `genome_annotation/` | funannotate output root (source for SETUP_INPUT) |
| `--outdir` | `results/function/` | Per-species per-tool result files (storeDir root) |
| `--tables` | `tables/` | Merged DuckDB-loadable CSV.gz output |
| `--taxon` | `""` (all) | Restrict to `RANK:VALUE`, e.g. `PHYLUM:Ascomycota` |
| `--n_test` | `0` (all) | Limit to first N samples after taxon filter |
| `--run_setup` | `true` | Run SETUP_INPUT symlink step |
| `--run_pfam` | `true` | Run hmmsearch vs Pfam-A |
| `--run_cazy` | `true` | Run dbcanlight CAZyme annotation |
| `--run_merops` | `true` | Run blastp vs MEROPS |
| `--run_signalp` | `true` | Run SignalP 6 (GPU) |
| `--run_tmhmm` | `true` | Run TMHMM |
| `--run_targetp` | `true` | Run TargetP 2 |
| `--run_idp` | `true` | Run AIUPred intrinsic disorder (GPU) |
| `--run_wolfpsort` | `true` | Run WoLF PSORT |
| `--run_predgpi` | `true` | Run predGPI |
| `--run_aa_freq` | `true` | Compute per-species AA frequencies |
| `--run_codon_freq` | `true` | Compute per-species codon frequencies |
| `--run_intergenic` | `true` | Compute intergenic distances |
| `--run_gene_stats` | `true` | Build gene structure tables from GFF3 |
| `--merge_all` | `true` | MERGE from all files in results/ (not just current run) |
| `--skip_merge` | `false` | Skip all MERGE steps |
| `--pfam_nodes` | `1` | SLURM nodes per hmmscan job (`>1` enables MPI mode) |

---

## SLURM resource allocation

| Label | Queue | CPUs | Memory | Time | Notes |
|---|---|---|---|---|---|
| `setup` | short | 1 | 1 GB | 30 min | SETUP_INPUT symlinks |
| `pfam` | highclock | 4 | 32 GB | 8 h | MPI mode if pfam_nodes > 1 |
| `cazy` | short | 16 | 8 GB | 2 h | |
| `merops` | short | 8 | 4 GB | 2 h | |
| `signalp` | short_gpu | 8 | 64 GB | 2 h | Requires `--gres=gpu:1` |
| `tmhmm` | short | 4 | 8 GB | 2 h | |
| `targetp` | epyc | 8 | 16 GB | 4 h | |
| `idp` | short_gpu | 4 | 48 GB | 2 h | Requires `--gres=gpu:1` |
| `wolfpsort` | short | 4 | 8 GB | 2 h | |
| `predgpi` | epyc | 4 | 16 GB | 4 h | |
| `genestats` | epyc | 4 | 96 GB | 48 h | Bulk stat jobs |
| `merge` | short | 4 | 24 GB | 2 h | |

---

## Skip/cache behavior

Per-species RUN processes use `storeDir`, so Nextflow skips a species automatically
if all declared output files already exist in `results/function/<tool>/` — even
on a fresh run without `-resume`.

| Process | Skip condition |
|---|---|
| `RUN_PFAM` | `pfam_hmmscan/<name>.pfam.gz` and `.tblout.gz` both exist |
| `RUN_CAZY` | `cazy/<name>/<name>.overview.tsv.gz`, `.cazymes.tsv.gz`, `.substrates.tsv.gz` all exist |
| `RUN_MEROPS` | `merops/<name>.blasttab.gz` exists |
| `RUN_SIGNALP` | `signalp/<name>.signalp.gff3.gz` and `.results.txt.gz` both exist |
| `RUN_TMHMM` | `tmhmm/<name>.tmhmm_short.tsv.gz` and `.tmhmm_results.tsv.gz` both exist |
| `RUN_TARGETP` | `targetP/<name>_summary.targetp2.gz` exists |
| `RUN_IDP` | `aiupred/<name>.aiupred.txt.gz`, `.idp.csv.gz`, `.idp_summary.csv.gz` all exist |
| `RUN_WOLFPSORT` | `wolfpsort/<name>.wolfpsort.results.txt.gz` exists |
| `RUN_PREDGPI` | `predgpi/<name>.predgpi.gff3.gz` exists |

To force re-running a species for a specific tool, delete its output file(s) from
`results/function/<tool>/` before re-submitting.

---

## Adding new species to an existing dataset

The typical workflow when adding new genomes to an already-annotated dataset:

```bash
# 1. Annotate new genomes first (funannotate.nf)
sbatch nextflow/run_funannotate.sh --taxon PHYLUM:Mucoromycota

# 2. Run functional annotation on new species only.
#    Existing species are automatically skipped (storeDir).
#    merge_all=true (default) means MERGE will include all existing + new files.
sbatch nextflow/run_functional.sh --taxon PHYLUM:Mucoromycota

# 3. Load updated tables into DuckDB
sbatch pipeline/db/02_build_functional.sh
```

To run new species and defer the merge until later (e.g., when multiple batches
are running in parallel):

```bash
# Run species without merging
sbatch nextflow/run_functional.sh --taxon PHYLUM:Mucoromycota --skip_merge true

# Merge all results once all batches are done
sbatch nextflow/run_functional.sh \
    --run_pfam false --run_cazy false --run_merops false \
    --run_signalp false --run_tmhmm false --run_targetp false \
    --run_idp false --run_wolfpsort false --run_predgpi false \
    --run_aa_freq false --run_codon_freq false \
    --run_intergenic false --run_gene_stats false
    # merge_all=true and skip_merge=false are the defaults
```

---

## Merge behavior: `merge_all` vs `merge_all=false`

| Mode | Which files go into `tables/*.csv.gz` |
|---|---|
| `--merge_all true` (default) | ALL `*.gz` files in `results/function/<tool>/` regardless of which species were run this session |
| `--merge_all false` | Only files produced by RUN processes in the current Nextflow invocation |

`merge_all=true` is the right default for incremental runs: it ensures `tables/`
always reflects the full dataset, not just the subset processed in a given run.
A sync barrier in the workflow ensures newly-run species are flushed to disk before
the directory glob fires.

---

## Required modules

| Tool | Module |
|---|---|
| Pfam HMM scan | `hmmer/3.4`, `db-pfam` |
| CAZyme annotation | `dbcanlight` |
| Protease families | `db-merops/124`, `ncbi-blast/2.16.0+` |
| Signal peptides | `signalp/6-gpu` (A100 GPU) |
| TM helices | `tmhmm` |
| Subcellular targeting | `targetp` |
| Disorder prediction | `aiupred` (A100 GPU) |
| Subcellular localization | `wolfpsort` |
| GPI anchors | `predgpi` |
| Gene/sequence statistics | `biopython` |

---

## Loading results into DuckDB

```bash
# After BFD.nf completes and tables/ is populated:
sbatch pipeline/db/02_build_functional.sh

# Interactive query
module load duckdb
duckdb functionalDB/function.duckdb
```

See `sql/schema.sql` for all table definitions and indexes. Key tables:
`pfam`, `cazy_overview`, `cazy_cazymes_hmm`, `merops`, `signalp`, `tmhmm`,
`targetp`, `idp`, `idp_summary`, `wolfpsort`, `predgpi`,
`aa_frequency`, `codon_frequency`, `gene_intergenic_distances`,
`gene_info`, `gene_introns`, `gene_exons`, `gene_CDS`,
`gene_transcripts`, `gene_trnas`, `gene_proteins`.

---

## Relationship to funannotate.nf

`funannotate.nf` produces `genome_annotation/<Species_Strain>/predict_results/`.
`BFD.nf` reads from that directory via its `SETUP_INPUT` step (symlinks into `input/`).
Always run `funannotate.nf` first, then `BFD.nf`.
