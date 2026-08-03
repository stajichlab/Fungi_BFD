# QUICKSTART — BFD Trace Profile Runtime Stats

Summarize runtime and memory statistics from Nextflow trace files across BFD runs.

## Quick Start

```bash
# Run on a directory of trace files (default glob: *_trace.*.txt)
python3 scripts-profile/bfd_trace_stats.py do_ANI/logs/nextflow/pre-scratch/

# Filter to specific trace files (e.g., only BFD runs, only ANI runs)
python3 scripts-profile/bfd_trace_stats.py do_ANI/logs/nextflow/pre-scratch/ --glob 'BFD_trace.*.txt'
python3 scripts-profile/bfd_trace_stats.py do_ANI/logs/nextflow/ --glob 'ANI_trace.*.txt'
```

## What It Produces

A single terminal report with:

| Section | Details |
|---|---|
| **File inventory** | List of matched trace files and row counts |
| **Status counts** | Total COMPLETED, FAILED, ABORTED tasks |
| **Runtime table** | Per-process: N, Tasks, Min, P10, Median, Mean, P90, P99, Max |
| **Memory table** | Per-process: N, Min, P10, Median, Mean, P90, P99, Max (RSS) |
| **Failure summary** | Failed process names with counts and percentages |
| **Overall summary** | Aggregate stats across all tasks |

## Example Output

```
========================================================================================================================
BFD TRACE STATISTICS — /bigdata/.../pre-scratch
Pattern: '*_trace.*.txt'
========================================================================================================================

Files analyzed: 4
  BFD_trace.2026-07-28_15_36_30.txt: 16091 rows
  BFD_trace.2026-08-01_16_00_14.txt: 465 rows
  ...

Overall status counts:
  COMPLETED: 27037
  FAILED: 191

STATISTICS BY BASE PROCESS NAME
Process                                                      N   Tasks     Min Time          P10       Median         Mean          P90          P99     Max Time
BUSCO_GENOME                                            14371   14286          2ms        9m 7s      15m 50s       16m 9s      23m 18s      32m 13s    1h 4m 50s
CALC_ASM_STATS                                          12991   12987          2ms         9.0s        20.2s        39.8s       1m 21s       5m 35s      12m 43s
MERGE_SAMPLES                                               1       1         4.7s         4.7s         4.7s         4.7s         4.7s         4.7s         4.7s
...
```

## Arguments

| Argument | Description | Default |
|---|---|---|
| `trace_dir` | Directory containing Nextflow `.txt` trace files | *(required)* |
| `--glob PATTERN` | Glob pattern to match trace files | `*_trace.*.txt` |

## Tips

- **Compare runs:** Point at different directories to compare BFD vs ANI, or pre-scratch vs main log folders
- **Filter by pipeline:** Use `--glob 'BFD_trace.*.txt'` to exclude ANI traces, or vice versa
- **Single run:** Point at a directory with just one trace file to profile a single Nextflow execution