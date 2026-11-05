#!/usr/bin/env python3
"""
compare_ani_methods.py — post-hoc join of two --ani_method result trees.

Run the compare_ANI workflow once per method (results are segregated under
results/ANI/<method>/<RANK>/), then join the two trees on the same groupings so
you can compare how two ANI tools agree on the %ID of every pair, and which is
faster, with no changes to the pipeline itself.

Reads, per (method, rank, group):
    <ani-dir>/<method>/<rank>/<group>/**/*.ani.tsv   (skani/mash/sourmash
    emit <group>.full.ani.tsv; fastani emits <group>.ani.tsv)
and Nextflow trace files (logs/nextflow/ANI_trace.*.txt by default) for the
wall-clock timing of each method's compare task per group.

Usage:
  python3 scripts/compare_ani_methods.py \
      --ani-dir results/ANI --compare GENUS \
      --methods skani,fastani \
      --thresholds 95.0,98.0

Outputs (under --out-dir, default <ani-dir>/method_comparison/<RANK>/):
  1. comparison_pairs.tsv    — every pair in either method, one row each:
                               group, query, reference, <A>_%ID, <B>_%ID,
                               delta_%ID (= B - A; blank where a method had no pair)
  2. comparison_summary.tsv  — per group: n common/only-A/only-B pairs,
                               mean/median/max |delta|, Pearson r over common
                               pairs, and per-threshold agreement counts
  3. method_timing.tsv       — per group: wall-clock realtime (sec) of each
                               method's compare task + speedup (A / B)
  4. <group>/<group>.pairs.tsv — per-group pairs tables (same cols as 1)
"""

import argparse
import csv
import glob
import math
import os
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

# Trace: process name is like 'SKANI_COMPARE (Fusarium [100 genomes])' (or the
# fully-qualified nlh variant 'COMPARE_ANI:ANI_COMPARE_METHOD:SKANI_COMPARE').
TRACE_GLOB    = 'logs/nextflow/ANI_trace.*.txt'
REALTIME_RE   = re.compile(r'([\d.]+)\s*(ms|s|m|h|d)')
TAG_SUFFIX_TS = re.compile(r'\s*n=\d+$')          # old tag 'Dothideomycetes n=3'
TAG_SUFFIX_N  = re.compile(r'\s+\[\d+\s+genomes\]$')  # current tag 'Fusarium [100 genomes]'


def parse_realtime_seconds(s):
    """Parse Nextflow realtime to seconds. Accepts '1h 20m 30s', '851ms' and
    bare integers (ms). Returns None if unparseable/empty."""
    if not s or s == '-':
        return None
    if s.strip().lstrip('+-').isdigit():
        return int(s) / 1000.0
    units = {'ms': 1 / 1000, 's': 1, 'm': 60, 'h': 3600, 'd': 86400}
    matches = REALTIME_RE.findall(s)
    if not matches:
        return None
    return sum(float(val) * units[unit] for val, unit in matches)


def load_pairs(tsv_path):
    """Return dict[(q, r)] = max_ani for a normalized ANI TSV (rows may be in
    either order; both directions already collapsed upstream but keep max)."""
    pairs = {}
    with open(tsv_path, newline='') as fh:
        for row in csv.reader(fh, delimiter='\t'):
            if len(row) < 3:
                continue
            q, r = row[0], row[1]
            try:
                ani = float(row[2])
            except ValueError:
                continue  # header line or malformed row
            if q == r:
                continue
            key = tuple(sorted([q, r]))
            if key not in pairs or pairs[key] < ani:
                pairs[key] = ani
    return pairs


def find_group_tsv(animethod_dir, group):
    """Find the per-group ANI TSV under <dir>/<group>/ (batches/ is only one
    possible location). Returns None if absent. Prefers *.full.ani.tsv over
    *.ani.tsv."""
    group_dir = Path(animethod_dir) / group
    candidates = sorted(group_dir.rglob('*.full.ani.tsv'))
    candidates += sorted(group_dir.rglob('*.ani.tsv'))
    for c in candidates:
        if group in c.name:
            return c
    return candidates[0] if candidates else None


def process_basename(name):
    """'COMPARE_ANI:ANI_COMPARE_METHOD:SKANI_COMPARE (tag)' → 'SKANI_COMPARE'."""
    base = name.split('(')[0].strip()
    return base.split(':')[-1].strip()


def extract_group(tag):
    """Nextflow task tag → group name ('Fusarium [100 genomes]' → 'Fusarium',
    legacy 'Dothideomycetes n=3' → 'Dothideomycetes')."""
    t = str(tag).strip()
    t = TAG_SUFFIX_N.sub('', t)
    t = TAG_SUFFIX_TS.sub('', t)
    return t


def load_timing(trace_glob, process_names):
    """Scan trace file(s); return dict group -> dict(process -> max realtime_s)
    over COMPLETED tasks of the given process names."""
    per_group = defaultdict(dict)
    files = sorted(glob.glob(os.path.expanduser(str(trace_glob))))
    for tf in files:
        try:
            with open(tf, newline='') as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                if not reader.fieldnames or 'name' not in reader.fieldnames:
                    continue
                for row in reader:
                    proc = process_basename(row.get('name', ''))
                    if proc not in process_names:
                        continue
                    if row.get('status') != 'COMPLETED':
                        continue
                    rt = parse_realtime_seconds(row.get('realtime', ''))
                    if rt is None:
                        continue
                    grp = extract_group(row.get('tag', ''))
                    if not grp:
                        continue
                    # max: prefer the longer (final) attempt for a retried task
                    if rt > per_group[grp].get(proc, 0):
                        per_group[grp][proc] = rt
        except OSError as e:
            print(f"Warning: cannot read trace {tf}: {e}", file=sys.stderr)
    if not files:
        print("Warning: no trace files matched; timing table will be empty "
              f"(glob: {trace_glob})", file=sys.stderr)
    return per_group


def pearson_r(xs, ys):
    """Pearson correlation coefficient of paired lists; None if <2 points or
    zero variance."""
    if len(xs) < 2:
        return None
    nx, ny = len(xs), len(ys)
    if nx != ny:
        return None
    mx, my = statistics.mean(xs), statistics.mean(ys)
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    dx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    dy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if dx == 0 or dy == 0:
        return None
    r = num / (dx * dy)
    return max(-1.0, min(1.0, r))


def fmt(x, nd=4):
    return 'n/a' if x is None else f'{x:.{nd}f}'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--ani-dir',    default='results/ANI',
                    help='Root of the method-segregated results trees (default: results/ANI)')
    ap.add_argument('--compare',    default='GENUS',
                    help='Rank subdir(s) to join, comma-separated or "all" (default: GENUS)')
    ap.add_argument('--methods',    default='skani,fastani',
                    help='Two methods to join, comma-separated, method A first (default: skani,fastani)')
    ap.add_argument('--thresholds', default='95.0,98.0',
                    help='Comma-separated ANI%% thresholds for agreement counts (default: 95.0,98.0)')
    ap.add_argument('--trace-glob', default=TRACE_GLOB,
                    help='Glob of Nextflow trace files for timing (default: logs/nextflow/ANI_trace.*.txt)')
    ap.add_argument('--out-dir',    default=None,
                    help='Output root (default: <ani-dir>/method_comparison)')
    args = ap.parse_args()

    methods = [m.strip() for m in args.methods.split(',') if m.strip()]
    if len(methods) != 2:
        ap.error('--methods must name exactly two methods (e.g. skani,fastani)')
    mA, mB = methods
    thresholds = [float(t) for t in args.thresholds.split(',') if t.strip()]

    ani_dir = Path(args.ani_dir)
    out_root = Path(args.out_dir) if args.out_dir else Path(ani_dir) / 'method_comparison'

    # discover common rank dirs
    ranks = None
    if args.compare.strip().lower() != 'all':
        ranks = [r.strip() for r in args.compare.split(',') if r.strip()]
    else:
        dA = ani_dir / mA
        ranks = sorted(dA.iterdir()) if dA.is_dir() else []
        ranks = [r.name for r in ranks if (dA / r.name).is_dir() and (ani_dir / mB / r.name).is_dir()]
    if not ranks:
        print(f"No rank dirs found under {ani_dir}/{mA}/ (or none shared with {mB})",
              file=sys.stderr)
        sys.exit(1)

    timing = load_timing(args.trace_glob, {f'{mA.upper()}_COMPARE', f'{mB.upper()}_COMPARE'})

    pair_header = ['group', 'query', 'reference', f'{mA}_%ID', f'{mB}_%ID', 'delta_%ID']

    for rank in ranks:
        dA, dB = Path(ani_dir) / mA / rank, Path(ani_dir) / mB / rank
        rank_out = out_root / rank

        groups_a = {g for g in os.listdir(dA) if (dA / g).is_dir()} if dA.is_dir() else set()
        groups_b = {g for g in os.listdir(dB) if (dB / g).is_dir()} if dB.is_dir() else set()
        shared = sorted(groups_a & groups_b)

        all_pairs_rows, summary_rows = [], []
        n_common = n_only_a = n_only_b = 0
        deltas_all = []

        for group in shared:
            tsvA = find_group_tsv(dA, group)
            tsvB = find_group_tsv(dB, group)
            if tsvA is None or tsvB is None:
                print(f"  skip {group}: missing {mA} tsv={tsvA} / {mB} tsv={tsvB}",
                      file=sys.stderr)
                continue
            pairsA = load_pairs(tsvA)
            pairsB = load_pairs(tsvB)

            grp_pairs, grp_deltas = [], []
            grp_only_a = grp_only_b = 0
            for key in sorted(set(pairsA) | set(pairsB)):
                a = pairsA.get(key)
                b = pairsB.get(key)
                if a is None:
                    grp_only_b += 1
                elif b is None:
                    grp_only_a += 1
                delta = (b - a) if (a is not None and b is not None) else None
                if a is not None and b is not None:
                    grp_deltas.append(delta)
                    deltas_all.append(delta)
                row = [group, key[0], key[1],
                       fmt(a), fmt(b), fmt(delta)]
                grp_pairs.append(row)
            n_common += len(grp_deltas)
            n_only_a += grp_only_a
            n_only_b += grp_only_b

            if grp_pairs:
                (rank_out / group).mkdir(parents=True, exist_ok=True)
                with open(rank_out / group / f'{group}.pairs.tsv', 'w', newline='') as fh:
                    w = csv.writer(fh, delimiter='\t')
                    w.writerow(pair_header)
                    w.writerows(grp_pairs)
            all_pairs_rows.extend(grp_pairs)

            common_a = [pairsA[k] for k in sorted(set(pairsA) & set(pairsB))]
            common_b = [pairsB[k] for k in sorted(set(pairsA) & set(pairsB))]
            r = pearson_r(common_a, common_b)
            mean_d = statistics.mean(grp_deltas) if grp_deltas else None
            med_d  = statistics.median(grp_deltas) if grp_deltas else None
            max_a  = max(abs(d) for d in grp_deltas) if grp_deltas else None

            tA = timing.get(group, {}).get(f'{mA.upper()}_COMPARE')
            tB = timing.get(group, {}).get(f'{mB.upper()}_COMPARE')
            speedup = (tA / tB) if (tA and tB and tB > 0) else None

            row = {'group': group,
                   'n_pairs_common': len(grp_deltas),
                   'n_pairs_only_A': grp_only_a,
                   'n_pairs_only_B': grp_only_b,
                   'mean_delta': mean_d, 'median_delta': med_d, 'max_abs_delta': max_a,
                   'pearson_r': r,
                   f'{mA}_sec': tA, f'{mB}_sec': tB,
                   'speedup_A_over_B': speedup}
            for th in thresholds:
                row[f'n_both>= {th}'] = sum(1 for k in set(pairsA) & set(pairsB)
                                            if pairsA[k] >= th and pairsB[k] >= th)
                row[f'n_only_{mA}>= {th}'] = sum(1 for k in set(pairsA) & set(pairsB)
                                                 if pairsA[k] >= th > pairsB[k])
                row[f'n_only_{mB}>= {th}'] = sum(1 for k in set(pairsA) & set(pairsB)
                                                 if pairsB[k] >= th > pairsA[k])
                row[f'n_both< {th}'] = sum(1 for k in set(pairsA) & set(pairsB)
                                           if pairsA[k] < th and pairsB[k] < th)
            summary_rows.append(row)

        rank_out.mkdir(parents=True, exist_ok=True)
        with open(rank_out / 'comparison_pairs.tsv', 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(pair_header)
            w.writerows(all_pairs_rows)

        sum_keys = ['group', 'n_pairs_common', 'n_pairs_only_A', 'n_pairs_only_B',
                    'mean_delta', 'median_delta', 'max_abs_delta', 'pearson_r',
                    f'{mA}_sec', f'{mB}_sec', 'speedup_A_over_B']
        with open(rank_out / 'comparison_summary.tsv', 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            header = list(sum_keys)
            for th in thresholds:
                header += [f'n_both>= {th}', f'n_only_{mA}>= {th}',
                           f'n_only_{mB}>= {th}', f'n_both< {th}']
            w.writerow(header)
            for r in summary_rows:
                cells = []
                for k in header:
                    v = r[k]
                    if v is None:
                        cells.append('')
                    elif k in ('mean_delta', 'median_delta', 'max_abs_delta',
                               'pearson_r', 'speedup_A_over_B'):
                        cells.append(fmt(v, 4 if k != 'speedup_A_over_B' else 3))
                    elif k.endswith('_sec'):
                        cells.append(fmt(v, 1))
                    else:
                        cells.append(v)
                w.writerow(cells)

        with open(rank_out / 'method_timing.tsv', 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['group', f'{mA}_sec', f'{mB}_sec', 'speedup_A_over_B'])
            for r in summary_rows:
                w.writerow([r['group'],
                            fmt(r[f'{mA}_sec'], 1), fmt(r[f'{mB}_sec'], 1),
                            fmt(r['speedup_A_over_B'], 3)])

        mean_d = statistics.mean(deltas_all) if deltas_all else None
        print(f"[{rank}] {len(shared)} shared groups | {n_common} common pairs, "
              f"{n_only_a} only-{mA}, {n_only_b} only-{mB} | "
              f"mean delta_%ID = {fmt(mean_d, 3)}")
        print(f"  wrote {rank_out}/comparison_pairs.tsv, comparison_summary.tsv, "
              f"method_timing.tsv")


if __name__ == '__main__':
    main()
