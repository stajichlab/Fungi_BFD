# Update SRA Query Cache with Platform Information

`nextflow/bin/update_platform_sra_query_cache.py` backfills the `platform` column
into pre-cached SRA query files that were written before platform tracking was added.

## What it does

Iterates all `*.sra_query.csv` files in `rnaseq_reads/sra_query/` and for each file that:
- has at least one data record (not header-only), and
- lacks a `platform` column

it calls `efetch -db sra -id <accessions> -format runinfo` to retrieve platform info
and rewrites the file in-place with `platform` appended as a new column.

BGI platform detection matches any value containing `BGI` or `BGISEQ` in the Platform field.

### Output files

| File | Contents |
|---|---|
| `UPDATED_SRA_CACHE.txt` | Filenames of cache files that were updated |
| `UPDATED_WITH_BGI.txt` | Species names (derived from filename) where any record used a BGI platform |

## Usage

```bash
module load ncbi_edirect
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD

# Dry run first to see what would change
python nextflow/bin/update_platform_sra_query_cache.py --dry-run

# Run for real (outputs go to current directory)
python nextflow/bin/update_platform_sra_query_cache.py

# Custom paths
python nextflow/bin/update_platform_sra_query_cache.py \
    --sra-dir rnaseq_reads/sra_query \
    --outdir results/
```

## Options

| Flag | Default | Description |
|---|---|---|
| `--sra-dir` | `rnaseq_reads/sra_query` | Directory containing `*.sra_query.csv` files |
| `--outdir` | `.` | Directory for summary output files |
| `--dry-run` | off | Report what would change without modifying files |

## Notes

- The script adds a 0.4 s pause between files and retries up to 3 times with exponential
  backoff to stay within NCBI rate limits.
- There are ~650 files that need updating, so run this in a `screen`/`tmux` session
  or wrap it in a short SLURM job to avoid it being killed on a login node.
- Files already containing a `platform` column are skipped unconditionally.
- Files with only a header line (no data rows) are also skipped.
