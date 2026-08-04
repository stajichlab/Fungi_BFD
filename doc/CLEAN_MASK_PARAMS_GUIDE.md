# Clean Genome + RepeatMask Params Guide

## Overview

This guide explains how to run the `earlgrey_mask` workflow using the new params-file system with support for genome suppression.

## Quick Start

### 1. Basic Run (All genomes)
```bash
./nextflow/run_clean_mask.sh
```

### 2. Test Mode (First 10 species)
```bash
./nextflow/run_clean_mask.sh --n-test 10
```

### 3. Single ASMID
```bash
./nextflow/run_clean_mask.sh --asmid GCA_000976515.2_Sc_YJM981_v1
```

### 4. With Resume
```bash
./nextflow/run_clean_mask.sh --resume
```

### 5. Custom Suppress List
```bash
./nextflow/run_clean_mask.sh --suppress-list /path/to/custom_suppress.txt
```

## Files Overview

### `nextflow/conf/params_clean_mask.yaml`
Main parameters file for the earlgrey_mask workflow.

**Key parameters:**

| Parameter | Default | Purpose |
|-----------|---------|---------|
| `pipeline` | `earlgrey_mask` | Pipeline selection (do not change) |
| `samples` | `./samples.csv` | Sample metadata file |
| `asm_stats` | `./tables/asm_stats.tsv.gz` | Assembly statistics |
| `genome_dir` | `./input_clean_genomes` | Input genomes directory |
| `genome_suffix` | `.fa` | Genome file suffix |
| `suppress_list` | `./data/curation/suppress.txt` | Suppression file |
| `cutoff_mb` | `200` | Minimum representative genome size (Mb) |
| `n_test` | `0` | Test mode: run first N species (0 = all) |
| `asmid` | `~` | Filter to single ASMID (empty = all) |
| `outdir` | `./results/repeatlibrary` | Output library directory |
| `masked_dir` | `./input_clean_genomes` | Output masked genomes directory |

### `data/curation/suppress.txt`
CSV file listing ASMIDs to exclude from processing.

**Format:** CSV with header row
```
ASMID,REASON
GCA_000149305.1_RO3,quality_issue
GCA_020081605.1_ASM2008160v1,in_progress
GCA_000987654.1_TEST,test_data
```

**Example:**
```csv
ASMID,REASON
GCA_000149305.1_RO3,quality_issue
GCA_020081605.1_ASM2008160v1,in_progress
GCA_000987654.1_TEST,test_data
```

**Reasons:**
- `quality_issue` — Poor assembly quality
- `in_progress` — Currently being reprocessed
- `manual_hold` — Held pending decision
- `test_data` — Test/scratch genome
- `redundant` — Duplicate of better assembly
- `corrupted` — Corrupted or incomplete file

### `nextflow/run_clean_mask.sh`
Wrapper script for running the earlgrey_mask workflow.

**Features:**
- Uses params file for configuration
- Supports command-line argument overrides
- Automatic trace/timeline/report logging
- Resume capability
- Pretty-printed output

## Workflow Steps

The earlgrey_mask workflow runs these processes:

1. **SELECT_REPS** — Select per-species representatives and members
   - Filters by minimum genome size (cutoff_mb)
   - Respects suppress.txt exclusions
   - Requires clean genomes to exist in genome_dir
   - Outputs: `misc/repeat_representatives.csv`

2. **EARLGREY_BUILD_LIB** — Build repeat library from representatives
   - Runs EarlGrey on each representative
   - Creates consensus repeat library
   - Masks representative genome
   - Output: `results/repeatlibrary/<SPECIES>/`

3. **REPEATMASK_STRAIN** — Apply library to member strains
   - Takes library from step 2
   - Masks each member strain genome
   - Output: `results/repeatlibrary/<SPECIES>/strains/`

4. **DELIVER_MASK** — Stage final masked genomes
   - Copies masked genomes to output location
   - Output: `input_clean_genomes/<ASMID>.masked.fasta.gz`

## Advanced Usage

### Modify Params Inline

Override any parameter from command line:

```bash
# Use different genome directory
./nextflow/run_clean_mask.sh --genome_dir ./data/custom_genomes

# Change size cutoff
./nextflow/run_clean_mask.sh --cutoff_mb 150

# Use alternate params file
nextflow run nextflow/main.nf -profile earlgrey \
  -params-file nextflow/conf/params_custom.yaml
```

### Workflow Control

**Process only species with reps over 250 Mb:**
```bash
./nextflow/run_clean_mask.sh --cutoff_mb 250
```

**Resume failed workflow:**
```bash
./nextflow/run_clean_mask.sh --resume
```

**Test on small dataset:**
```bash
./nextflow/run_clean_mask.sh --n-test 5 --genome_dir ./data/test_genomes
```

### Suppress Multiple Genomes

Edit `data/curation/suppress.txt`:

```bash
# Add to suppress list (CSV format)
echo "GCA_000976515.2_Sc_YJM981_v1,pending_qa" >> data/curation/suppress.txt

# Run (will skip that ASMID)
./nextflow/run_clean_mask.sh
```

### Create Custom Suppress List

```bash
# Create a new suppression file for specific project (CSV format)
cat > suppress_batch_2026_08.txt << 'EOF'
ASMID,REASON
GCA_000149305.1_RO3,in_progress
GCA_020081605.1_ASM2008160v1,in_progress
EOF

# Use it
./nextflow/run_clean_mask.sh --suppress-list suppress_batch_2026_08.txt
```

## Monitoring & Outputs

### During Execution

```bash
# Watch trace file (in another terminal)
tail -f logs/nextflow/earlgrey_trace.*.txt

# Monitor work directory
ls -lh work/earlgrey/

# Check Nextflow logs
tail -50 .nextflow.log
```

### After Completion

**Check results:**
```bash
# Repeat libraries built
ls results/repeatlibrary/*/

# Masked genomes
ls input_clean_genomes/*.masked.fasta.gz

# Reports
open logs/nextflow/earlgrey_report.*.html
open logs/nextflow/earlgrey_timeline.*.html
```

**View suppression log:**
```bash
# ASMIDs that were suppressed
grep "\[suppress\]" .nextflow.log | tail -20

# Species with libraries built
grep "\[info\].*Species" .nextflow.log
```

## Troubleshooting

### Issue: "representative genome not found"

**Cause:** Genome doesn't exist in `genome_dir`

**Solution:**
```bash
# Check if file exists
ls input_clean_genomes/GCA_000976515.2* 

# If missing, make sure clean genomes were generated first
# Or check if filename suffix is correct
grep "genome_suffix" nextflow/conf/params_clean_mask.yaml
```

### Issue: All species suppressed, nothing ran

**Cause:** suppress.txt excludes too many genomes

**Solution:**
```bash
# Check suppress list
grep -c "^[^#]" data/curation/suppress.txt

# Comment out suppressed entries temporarily
sed -i.bak 's/^GCA_/# GCA_/' data/curation/suppress.txt

# Run again
./nextflow/run_clean_mask.sh --n-test 5
```

### Issue: "no representative genome in clean genome directory"

**Cause:** Representative exists in samples.csv but clean genome file missing

**Solution:**
```bash
# Check what genomes are available
ls input_clean_genomes/ | head -20

# Check what the script is looking for
grep "REP_ASMID" misc/repeat_representatives.csv | head -5

# Either:
# 1. Generate clean genome first
# 2. Add that species to suppress.txt
# 3. Check if filename suffix is correct
```

## Integration with Splits

Since you've split code from runs (`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs`), adjust paths:

**For runs directory:**
```bash
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/run_clean_mask.sh

# Or with explicit paths
nextflow run \
  /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/main.nf \
  -profile earlgrey \
  -params-file /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/nextflow/conf/params_clean_mask.yaml \
  --genome_dir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/input_clean_genomes \
  --outdir /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/results/repeatlibrary
```

## Example Workflows

### Workflow 1: QA New Genomes

```bash
# 1. Add new genomes to input_clean_genomes/
# 2. Update samples.csv with new ASMIDs
# 3. Run on small test set first
./nextflow/run_clean_mask.sh --n-test 5

# 4. If successful, suppress old versions
echo "GCA_OLD_VERSION_1.0  redundant" >> data/curation/suppress.txt

# 5. Run full pipeline
./nextflow/run_clean_mask.sh --resume
```

### Workflow 2: Reprocess Failed Genomes

```bash
# 1. Create suppress list for successful genomes
grep "COMPLETED" logs/nextflow/earlgrey_trace.*.txt | \
  cut -d',' -f3 > suppress_completed.txt

# 2. Run only failed genomes
./nextflow/run_clean_mask.sh \
  --suppress-list suppress_completed.txt \
  --resume
```

### Workflow 3: Batch Processing

```bash
# 1. Process by species group
for cutoff in 150 200 250 300; do
  echo "Processing genomes > ${cutoff} Mb"
  ./nextflow/run_clean_mask.sh \
    --cutoff_mb $cutoff \
    -resume
done
```

## Key Script Changes

### SELECT_REPS (`nextflow/bin/select_repeat_representatives.py`)

**Added support for:**
- `--suppress-list` parameter to exclude ASMIDs
- `load_suppress_list()` function to read suppression file
- Automatic filtering of suppressed ASMIDs before processing
- Logging of suppressed count

### SELECT_REPS Process (`nextflow/modules/earlgrey/SELECT_REPS/main.nf`)

**Updated to:**
- Pass `params.suppress_list` to script
- Conditionally include `--suppress-list` arg only if defined

### Params File (`nextflow/conf/params_clean_mask.yaml`)

**New parameter:**
- `suppress_list: ./data/curation/suppress.txt`

## Quick Reference

| Task | Command |
|------|---------|
| Run all | `./nextflow/run_clean_mask.sh` |
| Test run | `./nextflow/run_clean_mask.sh --n-test 10` |
| Single ASMID | `./nextflow/run_clean_mask.sh --asmid GCA_xxx` |
| Resume | `./nextflow/run_clean_mask.sh --resume` |
| Custom suppress | `./nextflow/run_clean_mask.sh --suppress-list file.txt` |
| View help | `./nextflow/run_clean_mask.sh --help` |
| Custom size cutoff | `./nextflow/run_clean_mask.sh --cutoff_mb 250` |

## Questions?

Check these files for more details:
- Params format: `nextflow/conf/params_clean_mask.yaml`
- Suppress format: `data/curation/suppress.txt`
- Script usage: `./nextflow/run_clean_mask.sh --help`
- Workflow logic: `nextflow/workflows/earlgrey_mask.nf`
- Script logic: `nextflow/bin/select_repeat_representatives.py`
