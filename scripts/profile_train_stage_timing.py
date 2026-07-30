#!/usr/bin/env python3
"""
profile_train_stage_timing.py — parse funannotate-train.log stage timestamps
to estimate wall-clock time spent per pipeline stage (Trinity-GG, seqclean,
transcript alignment, PASA, TransDecoder, Kallisto, final selection).

Reads a random sample of logfiles/funannotate-train.log under given root dirs,
pairs consecutive '[MM/DD/YY HH:MM:SS]:' timestamp lines, and attributes the
elapsed gap to a canonical stage by regex-matching the line that STARTS the
gap. Aggregates per-stage total/median/p90/max seconds and % of total run
time across the sample.
"""

import argparse
import random
import re
import statistics
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path

TS_RE = re.compile(r'^\[(\d{2}/\d{2}/\d{2} \d{2}:\d{2}:\d{2})\]:\s*(.*)$')
TS_FMT = '%m/%d/%y %H:%M:%S'

# Ordered (pattern, canonical stage name) — first match wins.
STAGE_PATTERNS = [
    (r'^Building Hisat2 genome index', 'hisat2_build'),
    (r'^Aligning reads to genome using Hisat2', 'hisat2_align'),
    (r'^Running genome-guided Trinity', 'trinity_gg'),
    (r'transcripts derived from Trinity', 'trinity_gg_finalize'),
    (r'^Removing poly-A sequences from trinity transcripts using seqclean', 'seqclean_trinity'),
    (r'^/.*bin/seqclean trinity\.fasta', 'seqclean_trinity_cmd'),
    (r'^minimap2 -ax splice', 'minimap2_transcript_align'),
    (r'^Converting transcript alignments to GFF3', 'gff3_convert'),
    (r'^Converting Trinity transcript alignments to GFF3', 'gff3_convert_trinity'),
    (r'^Running PASA alignment step', 'pasa_launch'),
    (r'PASA assigned .* transcripts to .* loci', 'pasa_launch_done'),
    (r'^Getting PASA models for training with TransDecoder', 'transdecoder_training_set'),
    (r'^PASA finished\.', 'pasa_finished'),
    (r'^Building Kallisto index', 'kallisto_index'),
    (r'^Mapping reads using pseudoalignment in Kallisto', 'kallisto_quant'),
    (r'^Parsing expression value results', 'select_best_model'),
    (r'^Wrote .* PASA gene models', 'select_best_model_done'),
    (r'^Read trimming', 'read_trim_skip'),
    (r'^Read normalization', 'read_norm_skip'),
]


def classify(msg):
    for pat, name in STAGE_PATTERNS:
        if re.search(pat, msg):
            return name
    return None


def parse_log(path):
    """Return list of (canonical_stage_or_None, duration_seconds) for gaps
    between consecutive timestamped lines, plus total wall time (first to
    last timestamp)."""
    events = []
    try:
        with open(path, errors='replace') as fh:
            for line in fh:
                m = TS_RE.match(line)
                if not m:
                    continue
                ts = datetime.strptime(m.group(1), TS_FMT)
                events.append((ts, m.group(2).strip()))
    except Exception:
        return [], None

    if len(events) < 2:
        return [], None

    # Stateful attribution: many timestamped lines are raw subprocess command
    # echoes (e.g. the full Launch_PASA_pipeline.pl invocation, or a piped
    # hisat2|samtools command) that don't match any STAGE_PATTERNS themselves,
    # but occur *inside* a stage that was announced by an earlier line (e.g.
    # "Running PASA alignment step..."). Attributing each gap only by
    # matching the immediately-preceding line (stateless) misclassifies the
    # bulk of these multi-hour gaps as 'unclassified'. Instead, track a
    # sticky "current stage" that updates whenever a line matches a pattern,
    # and carries forward across unmatched lines until the next match.
    gaps = []
    current_stage = None
    for i in range(len(events) - 1):
        ts0, msg0 = events[i]
        ts1, _msg1 = events[i + 1]
        matched = classify(msg0)
        if matched:
            current_stage = matched
        dur = (ts1 - ts0).total_seconds()
        if dur < 0:
            continue  # clock oddity, skip
        gaps.append((current_stage, dur))

    total_wall = (events[-1][0] - events[0][0]).total_seconds()
    return gaps, total_wall


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--roots', nargs='+', required=True,
                     help='Root dirs to search for */logfiles/funannotate-train.log')
    ap.add_argument('--sample-size', type=int, default=400)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--min-lines', type=int, default=0,
                     help='Skip logs with fewer timestamped-or-total lines than this (0=no filter)')
    ap.add_argument('--out', required=True, help='Output CSV path (per-stage summary)')
    ap.add_argument('--out-per-run', default=None,
                     help='Optional CSV path for per-run total wall time + top stage')
    args = ap.parse_args()

    log_paths = []
    for root in args.roots:
        root_path = Path(root)
        if not root_path.is_dir():
            continue
        log_paths.extend(root_path.glob('*/logfiles/funannotate-train.log'))

    log_paths = sorted(set(log_paths))
    print(f"Found {len(log_paths)} funannotate-train.log files under {args.roots}", file=sys.stderr)

    rng = random.Random(args.seed)
    sample = log_paths if len(log_paths) <= args.sample_size else rng.sample(log_paths, args.sample_size)
    print(f"Sampling {len(sample)} logs (seed={args.seed})", file=sys.stderr)

    stage_durations = defaultdict(list)
    per_run_rows = []
    n_parsed = 0
    n_skipped = 0

    for p in sample:
        gaps, total_wall = parse_log(p)
        if total_wall is None:
            n_skipped += 1
            continue
        n_parsed += 1
        run_stage_totals = defaultdict(float)
        for stage, dur in gaps:
            key = stage if stage else 'unclassified'
            stage_durations[key].append(dur)
            run_stage_totals[key] += dur

        top_stage = max(run_stage_totals.items(), key=lambda kv: kv[1])[0] if run_stage_totals else 'NA'
        per_run_rows.append((str(p), total_wall, top_stage, run_stage_totals.get('pasa_launch', 0.0),
                              run_stage_totals.get('trinity_gg', 0.0) + run_stage_totals.get('hisat2_align', 0.0)
                              + run_stage_totals.get('hisat2_build', 0.0)))

    print(f"Parsed {n_parsed} logs, skipped {n_skipped} (unparseable/too few timestamps)", file=sys.stderr)

    # Per-stage summary
    rows = []
    grand_total = sum(sum(v) for v in stage_durations.values())
    for stage, durs in sorted(stage_durations.items(), key=lambda kv: -sum(kv[1])):
        durs_sorted = sorted(durs)
        n = len(durs_sorted)
        p90_idx = min(int(n * 0.9), n - 1)
        rows.append({
            'stage': stage,
            'n_occurrences': n,
            'total_seconds': sum(durs_sorted),
            'median_seconds': statistics.median(durs_sorted),
            'p90_seconds': durs_sorted[p90_idx],
            'max_seconds': durs_sorted[-1],
            'pct_of_sampled_wallclock': 100.0 * sum(durs_sorted) / grand_total if grand_total else 0.0,
        })

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    import csv
    with open(out_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['stage', 'n_occurrences', 'total_seconds', 'median_seconds',
                                            'p90_seconds', 'max_seconds', 'pct_of_sampled_wallclock'])
        w.writeheader()
        for r in rows:
            w.writerow(r)
    print(f"Wrote {len(rows)} stage rows to {out_path}", file=sys.stderr)

    if args.out_per_run:
        with open(args.out_per_run, 'w', newline='') as fh:
            w = csv.writer(fh)
            w.writerow(['path', 'total_wall_seconds', 'dominant_stage', 'pasa_launch_seconds', 'trinity_related_seconds'])
            for row in per_run_rows:
                w.writerow(row)
        print(f"Wrote {len(per_run_rows)} per-run rows to {args.out_per_run}", file=sys.stderr)

    # Print human summary
    print("\n" + "=" * 90)
    print(f"{'Stage':<28} {'N':<6} {'Total (h)':<10} {'Median (s)':<12} {'p90 (s)':<10} {'% of sampled wall':<10}")
    print("-" * 90)
    for r in rows:
        print(f"{r['stage']:<28} {r['n_occurrences']:<6} {r['total_seconds']/3600:<10.2f} "
              f"{r['median_seconds']:<12.1f} {r['p90_seconds']:<10.1f} {r['pct_of_sampled_wallclock']:<10.1f}")
    print("=" * 90)

    if per_run_rows:
        totals = [r[1] for r in per_run_rows]
        pasa = [r[3] for r in per_run_rows]
        trinity = [r[4] for r in per_run_rows]
        n = len(totals)
        print(f"\nPer-run total wall time: n={n} median={statistics.median(totals)/60:.1f} min "
              f"p90={sorted(totals)[int(n*0.9)]/60:.1f} min max={max(totals)/60:.1f} min")
        pasa_frac = [p_/t if t > 0 else 0 for p_, t in zip(pasa, totals) if t > 0]
        trinity_frac = [tg/t if t > 0 else 0 for tg, t in zip(trinity, totals) if t > 0]
        if pasa_frac:
            print(f"PASA (Launch_PASA_pipeline.pl) share of run wall time: median={statistics.median(pasa_frac)*100:.1f}% "
                  f"mean={100*sum(pasa_frac)/len(pasa_frac):.1f}%")
        if trinity_frac:
            print(f"Trinity-GG (hisat2+trinity) share of run wall time: median={statistics.median(trinity_frac)*100:.1f}% "
                  f"mean={100*sum(trinity_frac)/len(trinity_frac):.1f}%")


if __name__ == '__main__':
    main()
