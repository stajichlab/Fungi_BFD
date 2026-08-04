# Clean Mask Quick Start

## For Runs Directory Split Setup

Since you've split code from runs, here's the quick reference:

```
/bigdata/stajichlab/shared/projects/BFD/
├── Fungi_BFD/              ← Code & scripts (this repo)
└── Fungi_BFD_runs/         ← Data & execution
    ├── input_clean_genomes/
    ├── samples.csv
    ├── tables/
    ├── results/
    └── work/
```

## Run It

```bash
# Method 1: From runs directory (recommended)
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh

# Method 2: From code directory (specify run dir)
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
nextflow run nextflow/main.nf -profile earlgrey \
  -params-file nextflow/conf/params_clean_mask.yaml \
  --genome_dir ../Fungi_BFD_runs/input_clean_genomes \
  --outdir ../Fungi_BFD_runs/results/repeatlibrary
```

## Suppress Genomes

**File:** `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt` (CSV format)

```bash
# Add to suppress list (CSV: ASMID,REASON)
echo "GCA_000976515.2_Sc_YJM981_v1,pending_qa" >> \
  /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/data/curation/suppress.txt

# Run (will skip that ASMID)
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh
```

## Common Commands

```bash
# Test on 5 species
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh --n-test 5

# Process single species
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh \
  --asmid GCA_000976515.2_Sc_YJM981_v1

# Resume failed run
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh --resume

# Use different size cutoff
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh --cutoff_mb 250
```

## Configuration

**Params file:** `Fungi_BFD/nextflow/conf/params_clean_mask.yaml`
**Suppress list:** `Fungi_BFD/data/curation/suppress.txt`

Edit these to change defaults.

## Output

- **Repeat libraries:** `Fungi_BFD_runs/results/repeatlibrary/`
- **Masked genomes:** `Fungi_BFD_runs/input_clean_genomes/*.masked.fasta.gz`
- **Logs:** `Fungi_BFD_runs/logs/nextflow/`

## Full Documentation

See: `Fungi_BFD/CLEAN_MASK_PARAMS_GUIDE.md`
