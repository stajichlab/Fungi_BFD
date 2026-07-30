#!/usr/bin/env python3
"""
profile_nextflow_memory.py — profile Nextflow trace memory usage across runs
to identify under- and over-provisioned process memory allocations.

Scans trace files (funannotate_trace*.txt, *_trace*.txt) from multiple run
roots (default: ., do_annotation, do_annotation_asco, do_annotation_Rhod,
do_annotation_zoo, do_annotation_zzF), parses per-task RSS, dedupes on
task hash, groups by process name, and compares peak RSS against the attempt-1
memory allocation from nextflow/conf/profile_funannotate.config.

Outputs a CSV summary (process, n_completed, n_failed, n_oom, min/median/p90/max
rss, config memory, risk flag, suggested memory) and a printed table sorted by
risk level.
"""

import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path


# ─────────────────────────────────────────────────────────────────────────────
# Memory parsing (Nextflow units: KB, MB, GB, etc.)
# ─────────────────────────────────────────────────────────────────────────────

def parse_realtime_minutes(s):
    """Parse Nextflow realtime strings ('1h 20m 30s', '45m 10s', '30s', '500ms')
    to minutes. Returns None if unparseable/empty.
    """
    if not s or s == '-':
        return None
    s = s.strip()
    if not s:
        return None

    units = {'ms': 1 / 1000, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}
    matches = re.findall(r'([\d.]+)\s*(ms|s|m|h|d)', s)
    if not matches:
        return None

    total_seconds = sum(float(val) * units[unit] for val, unit in matches)
    return total_seconds / 60.0


def parse_memory_string(s):
    """Parse Nextflow memory strings ('3.5 MB', '64 GB', etc.) to bytes.

    Returns None if s is empty, '-', or unparseable.
    """
    if not s or s == '-':
        return None
    s = s.strip()
    if not s:
        return None

    units = {
        'KB': 1024,
        'MB': 1024 ** 2,
        'GB': 1024 ** 3,
        'TB': 1024 ** 4,
    }

    m = re.match(r'^([\d.]+)\s*([A-Z]+)$', s)
    if not m:
        return None

    val = float(m.group(1))
    unit = m.group(2)
    if unit not in units:
        return None

    return int(val * units[unit])


# ─────────────────────────────────────────────────────────────────────────────
# Trace file parsing
# ─────────────────────────────────────────────────────────────────────────────

# Trace files live directly under <root>/logs/nextflow/ (and, for legacy runs,
# <root>/logs_archive/**/logs/nextflow/). Do NOT recurse from the run root itself
# with '**' — these run directories also contain work/, rnaseq_reads/, and other
# multi-hundred-GB trees, so an unbounded recursive glob is extremely slow.
TRACE_GLOB_PATTERNS = [
    'logs/nextflow/*trace*.txt',
    'logs_archive/*/logs/nextflow/*trace*.txt',
    'logs_archive/logs/nextflow/*trace*.txt',
]


def parse_traces(roots):
    """Scan roots for trace files, parse rows, dedupe on hash.

    Returns dict: process_name -> [row_dict, ...] (across all roots, deduped)
    """
    expected_header = {'task_id', 'hash', 'name', 'status', 'exit', 'realtime',
                       '%cpu', 'rss', 'tag'}

    seen_hashes = set()
    rows_by_process = defaultdict(list)
    files_parsed = 0

    trace_files = []
    for root in roots:
        root_path = Path(root)
        if not root_path.is_dir():
            continue
        for pattern in TRACE_GLOB_PATTERNS:
            trace_files.extend(root_path.glob(pattern))

    for trace_file in sorted(set(trace_files)):
        try:
            with open(trace_file, 'r') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                if reader.fieldnames is None:
                    continue

                header_set = set(reader.fieldnames)
                if not expected_header.issubset(header_set):
                    # Not a trace file, skip
                    continue

                files_parsed += 1
                for row in reader:
                    hash_val = row.get('hash', '')
                    if not hash_val or hash_val in seen_hashes:
                        continue

                    seen_hashes.add(hash_val)

                    # Extract process name: everything before '('
                    name_full = row.get('name', '')
                    m = re.match(r'^([A-Za-z0-9_]+)', name_full)
                    process = m.group(1) if m else name_full

                    rows_by_process[process].append(row)
        except Exception as e:
            print(f"Warning: failed to parse {trace_file}: {e}", file=sys.stderr)

    print(f"Parsed {files_parsed} trace files, {sum(len(v) for v in rows_by_process.values())} rows, "
          f"{len(seen_hashes)} unique task hashes", file=sys.stderr)

    return rows_by_process


# ─────────────────────────────────────────────────────────────────────────────
# Config parsing (extract attempt-1 memory per process)
# ─────────────────────────────────────────────────────────────────────────────

def parse_config_memory(config_path):
    """Parse profile_funannotate.config; extract attempt-1 memory per process.

    Returns dict: process_name -> memory_bytes (or None if unparseable)
    """
    result = {}

    if not Path(config_path).is_file():
        return result

    try:
        text = Path(config_path).read_text()
    except Exception as e:
        print(f"Warning: failed to read {config_path}: {e}", file=sys.stderr)
        return result

    # withName: 'PROCESS' { ... } blocks
    # Match pattern: withName: 'PROCNAME' { ... memory = ... }

    # Simple approach: find each withName: 'PROC' block, extract process name,
    # then look for the next memory = ... line before the closing }

    with_name_re = re.compile(r"withName:\s*['\"]([A-Za-z0-9_]+)['\"]")

    lines = text.split('\n')
    for i, line in enumerate(lines):
        m = with_name_re.search(line)
        if not m:
            continue

        process = m.group(1)

        # Scan forward until we find memory = or the closing }
        for j in range(i, min(i + 20, len(lines))):
            block_line = lines[j]

            if re.search(r'^\s*}', block_line):
                # End of block
                break

            # Try to extract memory value
            # Pattern 1: memory = '64 GB' or memory = "64 GB"
            mem_match = re.search(r"memory\s*=\s*['\"]([^'\"]+)['\"]", block_line)
            if mem_match:
                mem_str = mem_match.group(1)
                mem_bytes = parse_memory_string(mem_str)
                if mem_bytes is not None:
                    result[process] = mem_bytes
                break

            # Pattern 2: memory = { ... ? 'X GB' : ... } (closure)
            # Extract the first quoted size (attempt-1 branch)
            if 'memory' in block_line and '{' in block_line:
                # Complex case: closure. Try to find the first quoted memory in this block
                # Scan this line and next few lines for quoted strings in the closure
                closure_text = block_line
                for k in range(j + 1, min(j + 10, len(lines))):
                    closure_text += ' ' + lines[k]
                    if '}' in lines[k]:
                        break

                # Find first quoted size in the closure
                sizes = re.findall(r"['\"](\d+\s*(?:KB|MB|GB|TB))['\"]", closure_text)
                if sizes:
                    mem_bytes = parse_memory_string(sizes[0])
                    if mem_bytes is not None:
                        result[process] = mem_bytes
                break

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Analysis and flagging
# ─────────────────────────────────────────────────────────────────────────────

def analyze_processes(rows_by_process, config_memory, min_realtime_min=0.0):
    """Analyze RSS per process, compute flags.

    min_realtime_min: exclude tasks with realtime below this threshold from the
    RSS/CPU distribution. Short tasks are often a process's fast-path/no-op branch
    (e.g. FUNANNOTATE_TRAIN's "already trained, just extract shared files" branch)
    rather than a genuine full run, and including them drags the percentile stats
    down without reflecting what a real run needs. n_completed/n_failed/n_oom are
    always computed on the full, unfiltered row set — those are correctness/failure
    signals, not memory-sizing signals.

    Returns list of dicts with keys:
      process, n_completed, n_failed, n_oom, n_rss_samples, n_short_excluded,
      min_rss_mb, median_rss_mb, p90_rss_mb, max_rss_mb, config_memory_gb, flag,
      suggested_memory_gb
    """

    results = []

    for process, rows in sorted(rows_by_process.items()):
        n_completed = sum(1 for r in rows if r.get('status') == 'COMPLETED')
        n_failed = sum(1 for r in rows if r.get('status') == 'FAILED')
        n_oom = sum(1 for r in rows if r.get('exit') in ('137', '140'))

        # Collect parseable RSS values, excluding short (likely fast-path) tasks
        rss_bytes = []
        n_short_excluded = 0
        for r in rows:
            rt_min = parse_realtime_minutes(r.get('realtime', ''))
            if min_realtime_min > 0 and rt_min is not None and rt_min < min_realtime_min:
                n_short_excluded += 1
                continue

            rss_str = r.get('rss', '')
            rss_b = parse_memory_string(rss_str)
            if rss_b is not None:
                rss_bytes.append(rss_b)

        if not rss_bytes:
            # No usable RSS data for this process
            continue

        rss_bytes.sort()
        n_samples = len(rss_bytes)
        min_rss = rss_bytes[0]
        median_rss = statistics.median(rss_bytes)
        max_rss = rss_bytes[-1]
        p90_idx = int(n_samples * 0.9)
        p90_rss = rss_bytes[p90_idx] if p90_idx < n_samples else rss_bytes[-1]

        config_mem_bytes = config_memory.get(process)
        config_mem_gb = config_mem_bytes / (1024 ** 3) if config_mem_bytes else None

        # Flag
        flag = "OK"
        if n_oom > 0 or (config_mem_bytes and max_rss > 0.85 * config_mem_bytes):
            flag = "OOM_RISK"
        elif config_mem_bytes and config_mem_bytes > (256 * 1024 * 1024) and \
             p90_rss < 0.4 * config_mem_bytes:
            flag = "OVERPROVISIONED"

        # Suggested memory: max_rss * 1.3, min 256 MB, round up to next GB
        suggested_bytes = max(max_rss * 1.3, 256 * 1024 * 1024)
        suggested_gb = int((suggested_bytes + (1024 ** 3) - 1) / (1024 ** 3))

        results.append({
            'process': process,
            'n_completed': n_completed,
            'n_failed': n_failed,
            'n_oom': n_oom,
            'n_rss_samples': n_samples,
            'n_short_excluded': n_short_excluded,
            'min_rss_mb': min_rss / (1024 ** 2),
            'median_rss_mb': median_rss / (1024 ** 2),
            'p90_rss_mb': p90_rss / (1024 ** 2),
            'max_rss_mb': max_rss / (1024 ** 2),
            'config_memory_gb': config_mem_gb if config_mem_gb else 'N/A',
            'flag': flag,
            'suggested_memory_gb': suggested_gb,
        })

    # Sort by risk: OOM_RISK first, then OVERPROVISIONED, then OK
    flag_order = {'OOM_RISK': 0, 'OVERPROVISIONED': 1, 'OK': 2}
    results.sort(key=lambda r: (flag_order.get(r['flag'], 99), r['process']))

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Output
# ─────────────────────────────────────────────────────────────────────────────

def write_csv(results, out_path):
    """Write results to CSV."""
    fieldnames = [
        'process', 'n_completed', 'n_failed', 'n_oom',
        'n_rss_samples', 'n_short_excluded',
        'min_rss_mb', 'median_rss_mb', 'p90_rss_mb', 'max_rss_mb',
        'config_memory_gb', 'flag', 'suggested_memory_gb'
    ]

    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            writer.writerow({k: r.get(k) for k in fieldnames})


def print_summary(results):
    """Print a summary table to stdout."""
    print("\n" + "=" * 100)
    print("Nextflow Memory Profile Summary")
    print("=" * 100)
    print(f"{'Process':<30} {'N':<8} {'Excl':<6} {'Flag':<20} {'p90 RSS (MB)':<14} {'Max RSS (MB)':<15} {'Config (GB)':<12} {'Suggested (GB)':<15}")
    print("-" * 100)

    for r in results:
        n = r['n_completed'] + r['n_failed']
        max_rss = r['max_rss_mb']
        p90_rss = r['p90_rss_mb']
        excl = r['n_short_excluded']
        config_gb = f"{r['config_memory_gb']:.1f}" if isinstance(r['config_memory_gb'], float) else r['config_memory_gb']
        print(f"{r['process']:<30} {n:<8} {excl:<6} {r['flag']:<20} {p90_rss:<14.1f} {max_rss:<15.1f} {config_gb:<12} {r['suggested_memory_gb']:<15}")

    # Summary by flag
    print("\n" + "-" * 100)
    flag_counts = {}
    for r in results:
        flag = r['flag']
        flag_counts[flag] = flag_counts.get(flag, 0) + 1

    print("Summary by risk level:")
    for flag in ['OOM_RISK', 'OVERPROVISIONED', 'OK']:
        count = flag_counts.get(flag, 0)
        if count > 0:
            print(f"  {flag}: {count}")

    print("=" * 100 + "\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument(
        '--roots',
        nargs='+',
        default=['.', 'do_annotation', 'do_annotation_asco', 'do_annotation_Rhod',
                 'do_annotation_zoo', 'do_annotation_zzF'],
        help='Root directories to scan for trace files'
    )
    ap.add_argument(
        '--config',
        default='nextflow/conf/profile_funannotate.config',
        help='Path to Nextflow config file'
    )
    ap.add_argument(
        '--out',
        default='analysis/nextflow_memory_profile/memory_profile_summary.csv',
        help='Output CSV path'
    )
    ap.add_argument(
        '--min-realtime-min',
        type=float,
        default=0.0,
        help=('Exclude tasks with realtime below this many minutes from the '
              'RSS/CPU stats (default: 0, no filtering). Use this to drop '
              'fast-path/no-op runs (e.g. a process that short-circuits when '
              'output already exists) that would otherwise pull percentile '
              'stats down without reflecting a genuine full run.')
    )
    args = ap.parse_args()

    # Parse traces
    rows_by_process = parse_traces(args.roots)
    if not rows_by_process:
        print("No trace data found.", file=sys.stderr)
        return 1

    # Parse config
    config_memory = parse_config_memory(args.config)

    # Analyze
    results = analyze_processes(rows_by_process, config_memory, args.min_realtime_min)

    # Output
    out_path = Path(args.out)
    write_csv(results, out_path)
    print(f"Wrote {len(results)} processes to {out_path}", file=sys.stderr)

    # Print summary
    print_summary(results)

    return 0


if __name__ == '__main__':
    sys.exit(main())
