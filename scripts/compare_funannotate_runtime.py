#!/usr/bin/env python3
"""
compare_funannotate_runtime.py — compare per-process wall-clock runtime between
two generations of the funannotate pipeline.

The funannotate pipeline has been migrating from the host module (`module load
funannotate/dev-1.9`, perl-based PASA/EVM/Trinity) to a Singularity container
(`funannotate-1.9.0b6.sif`, rust-optimized PASA/EVM/Trinity + AVX2 bowtie2).
This script scans Nextflow trace files (`logs/nextflow/funannotate_trace*.txt`),
splits tasks into two groups by the trace file's timestamp (a container go-live
cutoff date), and reports per-process runtime stats (n, median, p90, max realtime)
for each group side by side, plus matplotlib graphs.

Default container go-live cutoff = 2026-08-15 (FUNANNOTATE_PREDICT/ANNOTATE were
actually wired into the rust container in commit 33fb4b0 "containerize funannotate
and substeps" / 7bd6d9f "roll out containerized funannotate"; FUNANNOTATE_TRAIN
landed a day earlier on 08-14 in 695f47e). Trace runs >= this date are container.

Usage:
  compare_funannotate_runtime.py \
      --roots do_annotation do_annotation_asco do_annotation_Rhod \
      --cutoff 2026-08-13 \
      --outdir analysis/funannotate_runtime_comparison

Outputs (per --outdir):
  1. runtime_comparison.csv   — per-process (container vs module) median/p90/max
     realtime in minutes, plus ratio and % faster.
  2. runtime_by_process.csv   — long-form: one row per (process, group, task-observ)
  3. runtime_barplot.png      — grouped bar: median realtime per process, per group
  4. runtime_p90_barplot.png  — grouped bar: p90 realtime per process, per group
  5. runtime_boxplot.png      — boxplots per process, overlaid module vs container
"""

import argparse
import csv
import re
import statistics
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Trace globs relative to a run root. Do NOT recurse the whole run root (work/,
# rnaseq_reads/ etc. are huge). Same rationale as profile_nextflow_memory.py.
TRACE_GLOB_PATTERNS = [
    'logs/nextflow/funannotate_trace*.txt',
    'logs/nextflow/funannotate_trace.txt',
]
EXPECTED_HEADER = {'task_id', 'hash', 'name', 'status', 'exit', 'realtime',
                   '%cpu', 'rss', 'tag'}

TS_RE = re.compile(r'(\d{4})-(\d{2})-(\d{2})_')
REALTIME_RE = re.compile(r'([\d.]+)\s*(ms|s|m|h|d)')


def parse_realtime_seconds(s):
    """Parse Nextflow realtime strings ('1h 20m 30s', '45m 10s', '500ms') to
    seconds. Returns None if unparseable/empty."""
    if not s or s == '-':
        return None
    units = {'ms': 1 / 1000, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}
    matches = REALTIME_RE.findall(s)
    if not matches:
        return None
    return sum(float(val) * units[unit] for val, unit in matches)


def trace_date_from_name(path):
    """Extract the run date (2026-08-13) from a trace filename like
    funannotate_trace.2026-08-13_09_05_01.txt. Returns date object or None."""
    m = TS_RE.search(str(path))
    if not m:
        return None
    try:
        return datetime(int(m.group(1)), int(m.group(2)), int(m.group(3))).date()
    except ValueError:
        return None


def process_basename(name_full):
    """FUNANNOTATE:FUNANNOTATE_RNASEQ:SRA_QUERY_BATCH (...) -> SRA_QUERY_BATCH"""
    name_no_tag = re.sub(r'\s*\(.*$', '', name_full).strip()
    return name_no_tag.split(':')[-1].strip()


def collect_trace_files(roots):
    files = []
    for root in roots:
        root_path = Path(root)
        if not root_path.is_dir():
            continue
        for pattern in TRACE_GLOB_PATTERNS:
            for p in root_path.glob(pattern):
                if p.is_file():
                    files.append(p)
    return sorted(set(files))


def load_tasks(trace_paths):
    """Parse trace files into grouped row dicts keyed by (date, process_name).

    Returns dict: (date, process) -> list of (realtime_seconds, status, exit)
    and the set of completed process names. Dedupes on task hash across files.
    """
    seen_hashes = set()
    data = defaultdict(list)
    files_parsed = 0
    for trace_path in trace_paths:
        date = trace_date_from_name(trace_path)
        if date is None:
            continue
        try:
            with open(trace_path) as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                if reader.fieldnames is None:
                    continue
                if not EXPECTED_HEADER.issubset(set(reader.fieldnames)):
                    continue
                for row in reader:
                    h = row.get('hash', '')
                    if not h or h in seen_hashes:
                        continue
                    seen_hashes.add(h)
                    proc = process_basename(row.get('name', ''))
                    if not proc:
                        continue
                    rt = parse_realtime_seconds(row.get('realtime', ''))
                    if rt is None:
                        continue
                    data[(date, proc)].append((rt, row.get('status', ''),
                                                row.get('exit', '')))
            files_parsed += 1
        except Exception as e:
            print(f"Warning: failed to parse {trace_path}: {e}", file=sys.stderr)
    print(f"Parsed {files_parsed} trace files, "
          f"{len(seen_hashes)} unique task hashes", file=sys.stderr)
    return data


def group_by_generation(data, cutoff):
    """Split tasks into two generations by trace date vs cutoff.

    Returns (container_tasks, module_tasks) dicts:
      dict: process -> list of realtime_seconds (only COMPLETED attempts)
    plus per-group per-process counts of completed rows.
    """
    container = defaultdict(list)
    module = defaultdict(list)
    container_n = defaultdict(int)
    module_n = defaultdict(int)
    for (date, proc), rows in data.items():
        dst = container if date >= cutoff else module
        nctr = container_n if date >= cutoff else module_n
        for rt, status, exit_code in rows:
            if status != 'COMPLETED':
                continue
            dst[proc].append(rt)
            nctr[proc] += 1
    return container, module, container_n, module_n


def summarize_series(values):
    if not values:
        return None
    values = sorted(values)
    n = len(values)
    p90 = values[min(int(n * 0.9), n - 1)]
    return {
        'n': n,
        'total': sum(values),
        'median': statistics.median(values),
        'p90': p90,
        'max': max(values),
    }


def build_comparison(container, module):
    """Build a comparison table across all processes seen in either group."""
    processes = sorted(set(container) | set(module))
    rows = []
    for proc in processes:
        c = summarize_series(container.get(proc, []))
        m = summarize_series(module.get(proc, []))
        if c is None:
            ratio = pct_faster = None
        elif m is None:
            ratio = pct_faster = None
        elif m['median'] > 0:
            ratio = c['median'] / m['median']
            pct_faster = (1 - c['median'] / m['median']) * 100
        else:
            ratio = pct_faster = None
        rows.append({
            'process': proc,
            'container_n': c['n'] if c else 0,
            'container_median_min': (c['median'] / 60) if c else None,
            'container_p90_min': (c['p90'] / 60) if c else None,
            'container_max_min': (c['max'] / 60) if c else None,
            'module_n': m['n'] if m else 0,
            'module_median_min': (m['median'] / 60) if m else None,
            'module_p90_min': (m['p90'] / 60) if m else None,
            'module_max_min': (m['max'] / 60) if m else None,
            'container_vs_module_ratio': ratio,
            'container_percent_faster': pct_faster,
        })
    rows.sort(key=lambda r: (r['process']))
    return rows


def write_csv(rows, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ['process', 'container_n', 'container_median_min', 'container_p90_min',
              'container_max_min', 'module_n', 'module_median_min', 'module_p90_min',
              'module_max_min', 'container_vs_module_ratio', 'container_percent_faster']
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})


def print_summary(rows):
    def _fmt(v, suffix=''):
        if v is None:
            return '  -  '
        return f"{v:>6.1f}{suffix}"

    print("\n" + "=" * 112)
    print("Funannotate runtime comparison: CONTAINER (new) vs MODULE (old)")
    print("Median realtime per process, in minutes")
    print("=" * 112)
    hdr = (f"{'Process':<30}{'C:n':>5}{'C-med':>8}{'C-p90':>8}"
           f"{'M:n':>5}{'M-med':>8}{'M-p90':>8}{'C/M':>7}{'%faster':>9}")
    print(hdr)
    print("-" * 112)
    for r in rows:
        if r['process'] != 'FUNANNOTATE_PREDICT' and r['process'] != 'FUNANNOTATE_TRAIN':
            continue
        print(f"{r['process']:<30}"
              f"{r['container_n']:>5}{_fmt(r['container_median_min']):>8}"
              f"{_fmt(r['container_p90_min']):>8}"
              f"{r['module_n']:>5}{_fmt(r['module_median_min']):>8}"
              f"{_fmt(r['module_p90_min']):>8}"
              f"{_fmt(r['container_vs_module_ratio']):>7}"
              f"{_fmt(r['container_percent_faster']):>9}")
    print("-" * 112)
    print(f"{'ALL PROCESSES':<30}{'':>5}{'':>8}{'':>8}{'':>5}{'':>8}{'':>8}{'':>7}{'':>9}")
    for r in rows:
        if r['process'] in ('FUNANNOTATE_PREDICT', 'FUNANNOTATE_TRAIN'):
            continue
        print(f"{r['process']:<30}"
              f"{r['container_n']:>5}{_fmt(r['container_median_min']):>8}"
              f"{_fmt(r['container_p90_min']):>8}"
              f"{r['module_n']:>5}{_fmt(r['module_median_min']):>8}"
              f"{_fmt(r['module_p90_min']):>8}"
              f"{_fmt(r['container_vs_module_ratio']):>7}"
              f"{_fmt(r['container_percent_faster']):>9}")
    print("=" * 112)


def make_barplot(rows, processes, value_key, title, out_path, ylabel='minutes'):
    procs = [p for p in processes if any(r['process'] == p for r in rows)]
    container_vals = []
    module_vals = []
    for r in rows:
        if r['process'] in procs:
            container_vals.append(r[value_key] if r[value_key] is not None else 0)
            module_vals.append(r[value_key] if r[value_key] is not None else 0)
    if not procs:
        return
    x = np.arange(len(procs))
    width = 0.38
    fig, ax = plt.subplots(figsize=(0.45 * max(len(procs), 2) + 4, 6))
    ax.bar(x - width / 2, module_vals, width, label='MODULE (old, perl)', color='#c44e52')
    ax.bar(x + width / 2, container_vals, width, label='CONTAINER (new, rust)', color='#4c72b0')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xticks(x)
    ax.set_xticklabels(procs, rotation=45, ha='right')
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def make_boxplot(container, module, processes, out_path):
    fig, ax = plt.subplots(figsize=(0.5 * len(processes) + 4, 6))
    positions = []
    labels = []
    offset = 0
    for i, proc in enumerate(processes):
        for grp, style in (('C', dict(facecolor='#4c72b0', position=offset + i * 3 - 0.4)),
                           ('M', dict(facecolor='#c44e52', position=offset + i * 3 + 0.4))):
            data = (container.get(proc, []) if grp == 'C' else module.get(proc, []))
            data = [d / 60 for d in data]
            if not data:
                continue
            bp = ax.boxplot(data, positions=[style['position']], widths=0.6,
                            patch_artist=True, showfliers=False)
            bp['boxes'][0].set_facecolor(style['facecolor'])
            positions.append(style['position'])
        labels.append((offset + i * 3, proc))
        offset += 3
    ax.set_xticks([p for p, _ in labels])
    ax.set_xticklabels([l for _, l in labels], rotation=45, ha='right')
    ax.set_ylabel('realtime (minutes)')
    ax.set_title('Per-process realtime: CONTAINER (blue) vs MODULE (red)')
    import matplotlib.patches as mpatches
    ax.legend(handles=[mpatches.Patch(color='#4c72b0', label='CONTAINER'),
                       mpatches.Patch(color='#c44e52', label='MODULE')])
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--roots', nargs='+', required=True,
                    help='Run root dirs to scan for logs/nextflow/funannotate_trace*.txt')
    ap.add_argument('--cutoff', default='2026-08-15',
                    help='Container go-live date; trace runs >= this date are CONTAINER '
                         '(default 2026-08-15, when FUNANNOTATE_PREDICT/ANNOTATE were '
                         'actually wired into the rust container)')
    ap.add_argument('--outdir', default='analysis/funannotate_runtime_comparison',
                    help='Output directory for CSV + PNG')
    args = ap.parse_args()

    cutoff = datetime.fromisoformat(args.cutoff).date()

    trace_paths = collect_trace_files(args.roots)
    if not trace_paths:
        print("No trace files found under roots.", file=sys.stderr)
        return 1
    print(f"Found {len(trace_paths)} trace files", file=sys.stderr)

    data = load_tasks(trace_paths)
    container, module, c_n, m_n = group_by_generation(data, cutoff)
    print(f"Container tasks (COMPLETED): {sum(c_n.values())} across {len(container)} processes",
          file=sys.stderr)
    print(f"Module tasks (COMPLETED):    {sum(m_n.values())} across {len(module)} processes",
          file=sys.stderr)

    rows = build_comparison(container, module)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    write_csv(rows, outdir / 'runtime_comparison.csv')

    print_summary(rows)

    # Bar plots: only for processes observed in BOTH generations (meaningful A/B),
    # and only those with >= 3 completed samples in each group.
    both = [p for p in sorted(set(container) & set(module))
            if c_n[p] >= 3 and m_n[p] >= 3]
    for label, key, fname in (('median', 'container_median_min', 'runtime_barplot.png'),
                              ('p90', 'container_p90_min', 'runtime_p90_barplot.png')):
        make_barplot(rows, both, key,
                     f'Per-process {label.upper()} realtime (minutes)',
                     outdir / fname, ylabel='minutes')

    make_boxplot(container, module, both, outdir / 'runtime_boxplot.png')

    print(f"\nWrote comparison CSV and graphs to {outdir}", file=sys.stderr)
    print(f"Processes with >=3 samples in both generations, compared in plots: "
          f"{', '.join(both) if both else 'none'}", file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())