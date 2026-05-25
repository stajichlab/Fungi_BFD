# PHYling Phylogenomics Pipeline

Multi-locus phylogenomic tree construction for subsets of the BFD fungal genomes.
Uses [PHYling](https://github.com/stajichlab/PHYling) to identify and align BUSCO
single-copy markers, filter by informativeness, then build a partitioned ML tree
with IQ-TREE, RAxML, or FastTree.

Pipeline: `nextflow/phyling.nf` — launcher: `nextflow/run_phyling.sh`

---

## Prerequisites

Input FASTAs must already exist in `input/pep/` or `input/cds/` — these are
populated by the `SETUP_INPUT` step of `BFD.nf`. Run `BFD.nf` first (or at least
`sbatch nextflow/run_functional.sh` with all annotation tools disabled and
`--run_setup true`).

```bash
# If input/ is not yet populated, run setup only:
sbatch nextflow/run_functional.sh \
    --run_pfam false --run_cazy false --run_merops false \
    --run_signalp false --run_tmhmm false --run_targetp false \
    --run_idp false --run_wolfpsort false --run_predgpi false \
    --run_aa_freq false --run_codon_freq false \
    --run_intergenic false --run_gene_stats false
```

---

## Sample sheet: `phylo.csv`

By default the pipeline reads `phylo.csv` in the project root. If `phylo.csv` is
absent it falls back to `samples.csv`. `phylo.csv` should be a curated subset of
`samples.csv` — for example, a particular clade, an outgroup-inclusive set, or a
set of high-quality assemblies. It must have the same column format as `samples.csv`.

```bash
# Create phylo.csv for Dothideomycetes + outgroup Eurotiomycetes
head -1 samples.csv > phylo.csv
awk -F',' '$7=="Ascomycota" && ($9=="Dothideomycetes" || $9=="Eurotiomycetes")' \
    samples.csv >> phylo.csv
```

Override on the command line with `--samples my_subset.csv`.

---

## Output structure

Results land under `phylogeny/` organised by taxon group, sequence type, and markerset:

```
phylogeny/
  {taxon_slug}/            ← value from --taxon (e.g. Ascomycota), or 'all'
    {seq_type}/            ← 'protein' or 'cds'
      {markerset}/
        align/             ← phyling align output (one .fa per marker)
        filter/            ← top-N marker MSAs selected by treeness/RCV
        tree/
          concat.treefile  ← ML tree in Newick format  ← PRIMARY RESULT
          concat.partition ← partition file
  markersets/              ← cached extracted markerset HMM directories
```

The **primary result file** is `tree/concat.treefile` — open in FigTree, iTOL, or
ggtree for visualisation.

---

## Running the pipeline

All examples are submitted from the **project root** directory.

### Full run (all taxa in phylo.csv)

```bash
sbatch nextflow/run_phyling.sh
```

### Restrict to a taxonomic group

The `--taxon RANK:VALUE` flag filters `phylo.csv` to the specified group and names
the output directory after the value (e.g. `phylogeny/Ascomycota/`).

```bash
sbatch nextflow/run_phyling.sh --taxon PHYLUM:Ascomycota
sbatch nextflow/run_phyling.sh --taxon PHYLUM:Basidiomycota
sbatch nextflow/run_phyling.sh --taxon PHYLUM:Mucoromycota
sbatch nextflow/run_phyling.sh --taxon CLASS:Dothideomycetes
sbatch nextflow/run_phyling.sh --taxon CLASS:Sordariomycetes
sbatch nextflow/run_phyling.sh --taxon CLASS:Agaricomycetes
sbatch nextflow/run_phyling.sh --taxon ORDER:Pleosporales
sbatch nextflow/run_phyling.sh --taxon GENUS:Aspergillus
```

### Use CDS instead of protein

```bash
sbatch nextflow/run_phyling.sh --seq_type cds
sbatch nextflow/run_phyling.sh --taxon CLASS:Dothideomycetes --seq_type cds
```

### Multiple markersets (run in parallel)

Each markerset runs as an independent parallel branch producing its own subtree
under `phylogeny/{taxon}/{seq_type}/{markerset}/`. Quote multi-value strings.

```bash
# Two broad markersets for cross-validation
sbatch nextflow/run_phyling.sh \
    --markerset "fungi_odb12,ascomycota_odb12"

# Class-level precision + broad context
sbatch nextflow/run_phyling.sh --taxon CLASS:Dothideomycetes \
    --markerset "fungi_odb12,ascomycota_odb12,dothideomycetes_odb12"

# Basidiomycota with class resolution
sbatch nextflow/run_phyling.sh --taxon PHYLUM:Basidiomycota \
    --markerset "fungi_odb12,basidiomycota_odb12,agaricomycetes_odb12"
```

### Change tree-building method

```bash
# IQ-TREE (default, recommended for publication)
sbatch nextflow/run_phyling.sh --tree_method iqtree

# RAxML-NG
sbatch nextflow/run_phyling.sh --tree_method raxml

# FastTree (fastest; good for pilots and large datasets)
sbatch nextflow/run_phyling.sh --tree_method ft
```

### Pilot runs

```bash
# First 20 taxa, FastTree (quick sanity check)
sbatch nextflow/run_phyling.sh --n_test 20 --tree_method ft

# First 50 Ascomycota with IQ-TREE
sbatch nextflow/run_phyling.sh --taxon PHYLUM:Ascomycota --n_test 50

# Custom sample sheet
sbatch nextflow/run_phyling.sh --samples my_subset.csv
```

### Resume an interrupted run

```bash
sbatch nextflow/run_phyling.sh -resume
sbatch nextflow/run_phyling.sh --taxon PHYLUM:Ascomycota -resume
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--samples` | `phylo.csv` (falls back to `samples.csv`) | Sample sheet |
| `--seq_type` | `protein` | `protein` (reads `input/pep/`) or `cds` (reads `input/cds/`) |
| `--markerset` | `fungi_odb12` | Comma-separated BUSCO lineage name(s) |
| `--markerset_db` | `/srv/projects/db/BUSCO/v12/lineages` | Path to BUSCO lineage database |
| `--tree_method` | `iqtree` | `ft`, `iqtree`, or `raxml` |
| `--top_n` | `50` | Top-N markers to retain after treeness/RCV filtering |
| `--taxon` | `""` (all) | `RANK:VALUE` filter on sample sheet |
| `--n_test` | `0` (all) | Limit to first N taxa (applied after `--taxon`) |
| `--phylo_outdir` | `phylogeny/` | Root output directory |

---

## Available BUSCO fungal lineages

Lineages at `/srv/projects/db/BUSCO/v12/lineages/`. Pre-extracted directories
(✓) are used directly; tarballs-only (tar) are extracted into `phylogeny/markersets/`
on first use and cached there for subsequent runs.

| Lineage | Scope | Status |
|---|---|---|
| `fungi_odb12` | All Fungi (broadest) | tar — extracted on first run |
| `ascomycota_odb12` | Ascomycota | ✓ ready |
| `basidiomycota_odb12` | Basidiomycota | ✓ ready |
| `mucoromycota_odb12` | Mucoromycota | ✓ ready |
| `agaricomycetes_odb12` | Agaricomycetes | ✓ ready |
| `dothideomycetes_odb12` | Dothideomycetes | ✓ ready |
| `sordariomycetes_odb12` | Sordariomycetes | ✓ ready |
| `tremellomycetes_odb12` | Tremellomycetes | ✓ ready |
| `pucciniomycetes_odb12` | Pucciniomycetes | ✓ ready |
| `xylariales_odb12` | Xylariales | ✓ ready |
| `onygenales_odb12` | Onygenales | ✓ ready |
| `ajellomycetaceae_odb12` | Ajellomycetaceae | ✓ ready |
| `aspergillus_odb12` | Aspergillus | ✓ ready |
| `agaricales_odb12` | Agaricales | tar — extracted on first run |
| `eurotiomycetes_odb12` | Eurotiomycetes | tar — extracted on first run |
| `hypocreales_odb12` | Hypocreales | tar — extracted on first run |
| `leotiomycetes_odb12` | Leotiomycetes | tar — extracted on first run |
| `saccharomycetes_odb12` | Saccharomycetes | tar — extracted on first run |
| `saccharomycetaceae_odb12` | Saccharomycetaceae | tar — extracted on first run |

**Recommended markerset combinations:**

| Study scope | Recommended `--markerset` |
|---|---|
| All fungi (broad) | `fungi_odb12` |
| Ascomycota | `fungi_odb12,ascomycota_odb12` |
| Basidiomycota | `fungi_odb12,basidiomycota_odb12` |
| Mucoromycota | `fungi_odb12,mucoromycota_odb12` |
| Dothideomycetes | `fungi_odb12,ascomycota_odb12,dothideomycetes_odb12` |
| Sordariomycetes | `fungi_odb12,ascomycota_odb12,sordariomycetes_odb12` |
| Agaricomycetes | `fungi_odb12,basidiomycota_odb12,agaricomycetes_odb12` |

---

## Pipeline steps (what runs inside each job)

### 1. `MARKERSET_PREPARE` (per markerset, cached)

Resolves the markerset HMM directory. Runs once per markerset across all pipeline
invocations; subsequent runs find the cached result in `phylogeny/markersets/` and
skip this step entirely.

- Pre-extracted directory → creates a symlink in `phylogeny/markersets/`
- Tarball only → extracts to `phylogeny/markersets/{name}/`

### 2. `PHYLING_ALIGN` (per markerset)

Runs `phyling align`. Before calling phyling, the process renames staged FASTAs to
clean taxon labels — stripping `.proteins` and `.cds-transcripts` suffixes so the
tree leaf labels read as `{Species_Strain}` rather than `{Species_Strain}.proteins`.

```
phyling align -I staged/ -m markerset_dir/ -o align/ --seqtype {dna|pep} -t 32
```

Output: one `.fa` MSA file per marker found in ≥ 4 taxa.

**Resources:** 32 CPUs, 256 GB RAM, `highmem` queue, 24 h walltime.
This is the heaviest step — hmmsearch against all BUSCO HMMs for every taxon.

### 3. `PHYLING_FILTER` (per markerset)

Selects the top-N most informative markers using the treeness/RCV score computed
by PhyKIT. Lower-scoring markers (compositionally biased or poorly resolved) are
discarded. Default `--top_n 50`; increase for denser datasets, decrease for sparse ones.

```
phyling filter -I align/ -n 50 -o filter/ --seqtype {dna|pep} -t 16
```

**Resources:** 16 CPUs, 32 GB RAM, `short` queue, 4 h.

### 4. `PHYLING_TREE` (per markerset)

Concatenates the filtered marker MSAs, writes a partition file, then builds a
partitioned ML tree with the chosen method.

```
phyling tree -I filter/ -M iqtree --concat --partition -o tree/ --seqtype {dna|pep} -t 32
```

**Resources:** 32 CPUs, 128 GB RAM, `epyc` queue, 48 h.

---

## SLURM resource summary

| Process | CPUs | Memory | Queue | Time |
|---|---|---|---|---|
| `MARKERSET_PREPARE` | 1 | 2 GB | short | 30 min |
| `PHYLING_ALIGN` | 32 | 256 GB | highmem | 24 h |
| `PHYLING_FILTER` | 16 | 32 GB | short | 4 h |
| `PHYLING_TREE` | 32 | 128 GB | epyc | 48 h |

---

## Monitoring

```bash
# Watch the Nextflow log
tail -f .nextflow.log

# Check SLURM jobs
squeue -u $USER

# Check a failed process
ls work/??/*/
cat work/<hash>/.command.log   # stdout + stderr for that process
cat work/<hash>/.command.sh    # exact shell command that was run
```

Nextflow trace, report, and timeline HTML files are written to `logs/nextflow/`
after each run.

---

## Troubleshooting

**`phyling align` finds 0 markers**
- Confirm `input/pep/` or `input/cds/` is populated (run BFD.nf SETUP_INPUT first)
- Confirm `--seq_type` matches the files present (`protein` → `.proteins.fa`, `cds` → `.cds-transcripts.fa`)
- Check the markerset name is spelled exactly as listed above

**`phyling filter` produces empty output**
- Too few taxa share any single marker — reduce `--top_n` (e.g. `--top_n 20`)
- Very small pilot datasets (`--n_test 10`) may not have enough taxa for reliable filtering

**`phyling tree` is slow**
- Normal for IQ-TREE and RAxML on large concatenated alignments (hundreds of taxa, 50 markers)
- Switch to `--tree_method ft` for a rapid exploratory tree; run IQ-TREE separately once satisfied

**`fungi_odb12` extraction takes a long time on first run**
- Expected — the tarball unpacks ~4,000 HMM files
- `storeDir` caches it; all subsequent runs skip this step entirely

**Resume after failure**
- Always add `-resume` — Nextflow caches completed processes in `work/` and skips them
- Do not delete `work/` between runs

---

## Visualising the tree

The output `tree/concat.treefile` is in Newick format and can be opened directly in:

- **FigTree** — `module load figtree; figtree tree/concat.treefile`
- **iTOL** — upload at [itol.embl.de](https://itol.embl.de) for annotated interactive display
- **ggtree (R)** — `library(ggtree); tr <- read.tree("concat.treefile"); ggtree(tr)`
- **ETE3 (Python)** — `from ete3 import Tree; t = Tree("concat.treefile")`
