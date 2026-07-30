#!/usr/bin/env python3
"""
profile_predict_stage_timing.py — parse funannotate-predict.log stage timestamps
to estimate wall-clock time spent per pipeline stage (RNA-seq hints prep,
GeneMark-ES self-training, Augustus training, Augustus gene prediction, SNAP
train+predict, EVM, tRNA/tbl2asn finishing), mirroring
scripts/profile_train_stage_timing.py's approach for funannotate train.

Purpose: determine which stages inside `funannotate predict` are actually
expensive, to evaluate the cost premise behind sharing ab-initio (AUGUSTUS/
SNAP/GeneMark-ES) parameters across strains (todo/species_level_abinitio_reuse.md).
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

# Ordered (pattern, canonical stage name) — first match wins. Derived from a
# manual read of a real funannotate-predict.log
# (genome_annotation/Aaosphaeria_arxii_CBS_175.79/logfiles/funannotate-predict.log).
STAGE_PATTERNS = [
    (r'^Loading genome assembly', 'genome_load'),
    (r'transcript alignments from:', 'transcript_alignment_parse'),
    (r'^Creating transcript EVM alignments', 'transcript_evm_hints'),
    (r'^Extracting hints from RNA-seq BAM file using bam2hints', 'rnaseq_hints_extract'),
    (r'filterIntronsFindStrand\.pl', 'rnaseq_hints_filter'),
    (r'^join_mult_hints\.pl', 'rnaseq_hints_join'),
    (r'^Running GeneMark-ES on assembly', 'genemark_es_train'),
    (r'predictions from GeneMark', 'genemark_es_done'),
    (r'^Filtering PASA data for suitable training set', 'pasa_filter_training_set'),
    (r'funannotate-BUSCO2\.py', 'augustus_busco_seed_run'),
    (r'^Training Augustus using PASA gene models', 'augustus_train_pasa'),
    (r'^Training Augustus using BUSCO gene models', 'augustus_train_busco'),
    (r'Augustus initial training results', 'augustus_train_done'),
    (r'^Running Augustus gene prediction using .* parameters', 'augustus_predict_run'),
    (r'high quality predictions from Augustus', 'augustus_predict_done'),
    (r'^Running SNAP gene prediction', 'snap_train_and_predict'),
    (r'predictions from SNAP', 'snap_done'),
    (r'^Launching EVM via funannotate-runEVM\.py', 'evm_launch'),
    (r'total gene models from EVM', 'evm_done'),
    (r'^Predicting tRNAs', 'trna_predict'),
    (r'^Generating GenBank tbl annotation file', 'tbl_generate'),
    (r'^Converting to final Genbank format', 'tbl2asn_convert'),
    (r'^Funannotate predict is finished', 'finished_marker'),
]

# Stages that ab-initio parameter reuse (funannotate predict -p params.json)
# could plausibly skip or shrink: the self-training portions only, not the
# per-strain prediction passes that use whatever model (fresh or reused).
ABINITIO_TRAINING_STAGES = {'genemark_es_train', 'augustus_train_pasa', 'augustus_train_busco',
                             'augustus_busco_seed_run', 'snap_train_and_predict'}
AUGUSTUS_PASA_STAGES = {'augustus_train_pasa'}
AUGUSTUS_BUSCO_STAGES = {'augustus_train_busco', 'augustus_busco_seed_run'}


def classify(msg):
    for pat, name in STAGE_PATTERNS:
        if re.search(pat, msg):
            return name
    return None


def parse_log(path):
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
            continue
        gaps.append((current_stage, dur))

    total_wall = (events[-1][0] - events[0][0]).total_seconds()
    return gaps, total_wall


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--roots', nargs='+', required=True,
                     help='Root dirs to search for */logfiles/funannotate-predict.log')
    ap.add_argument('--sample-size', type=int, default=400)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--out', required=True, help='Output CSV path (per-stage summary)')
    ap.add_argument('--out-per-run', default=None,
                     help='Optional CSV path for per-run total wall time + ab-initio-training fraction')
    args = ap.parse_args()

    log_paths = []
    for root in args.roots:
        root_path = Path(root)
        if not root_path.is_dir():
            continue
        log_paths.extend(root_path.glob('*/logfiles/funannotate-predict.log'))

    log_paths = sorted(set(log_paths))
    print(f"Found {len(log_paths)} funannotate-predict.log files under {args.roots}", file=sys.stderr)

    rng = random.Random(args.seed)
    sample = log_paths if len(log_paths) <= args.sample_size else rng.sample(log_paths, args.sample_size)
    print(f"Sampling {len(sample)} logs (seed={args.seed})", file=sys.stderr)

    stage_durations = defaultdict(list)
    per_run_rows = []
    n_parsed = 0
    n_skipped = 0

    for p in sample:
        gaps, total_wall = parse_log(p)
        if total_wall is None or total_wall <= 0:
            n_skipped += 1
            continue
        n_parsed += 1
        run_stage_totals = defaultdict(float)
        for stage, dur in gaps:
            key = stage if stage else 'unclassified'
            stage_durations[key].append(dur)
            run_stage_totals[key] += dur

        abinitio_total = sum(v for k, v in run_stage_totals.items() if k in ABINITIO_TRAINING_STAGES)
        hints_total = (run_stage_totals.get('rnaseq_hints_extract', 0.0)
                       + run_stage_totals.get('rnaseq_hints_filter', 0.0)
                       + run_stage_totals.get('rnaseq_hints_join', 0.0))
        augustus_pasa_total = sum(v for k, v in run_stage_totals.items() if k in AUGUSTUS_PASA_STAGES)
        augustus_busco_total = sum(v for k, v in run_stage_totals.items() if k in AUGUSTUS_BUSCO_STAGES)
        evidence_mode = 'pasa' if augustus_pasa_total > 0 else ('busco' if augustus_busco_total > 0 else 'unknown')
        top_stage = max(run_stage_totals.items(), key=lambda kv: kv[1])[0] if run_stage_totals else 'NA'
        per_run_rows.append((str(p), total_wall, top_stage, abinitio_total, hints_total,
                              run_stage_totals.get('genemark_es_train', 0.0),
                              augustus_pasa_total, augustus_busco_total, evidence_mode,
                              run_stage_totals.get('snap_train_and_predict', 0.0)))

    print(f"Parsed {n_parsed} logs, skipped {n_skipped} (unparseable/zero-duration)", file=sys.stderr)

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
            w.writerow(['path', 'total_wall_seconds', 'dominant_stage', 'abinitio_training_seconds',
                        'rnaseq_hints_seconds', 'genemark_es_train_seconds', 'augustus_pasa_seconds',
                        'augustus_busco_seconds', 'augustus_evidence_mode', 'snap_train_and_predict_seconds'])
            for row in per_run_rows:
                w.writerow(row)
        print(f"Wrote {len(per_run_rows)} per-run rows to {args.out_per_run}", file=sys.stderr)

    print("\n" + "=" * 90)
    print(f"{'Stage':<28} {'N':<6} {'Total (h)':<10} {'Median (s)':<12} {'p90 (s)':<10} {'% of sampled wall':<10}")
    print("-" * 90)
    for r in rows:
        print(f"{r['stage']:<28} {r['n_occurrences']:<6} {r['total_seconds']/3600:<10.2f} "
              f"{r['median_seconds']:<12.1f} {r['p90_seconds']:<10.1f} {r['pct_of_sampled_wallclock']:<10.1f}")
    print("=" * 90)

    if per_run_rows:
        totals = [r[1] for r in per_run_rows]
        abinitio = [r[3] for r in per_run_rows]
        hints = [r[4] for r in per_run_rows]
        n = len(totals)
        print(f"\nPer-run total wall time: n={n} median={statistics.median(totals)/60:.1f} min "
              f"p90={sorted(totals)[int(n*0.9)]/60:.1f} min max={max(totals)/60:.1f} min")
        abinitio_frac = [a/t if t > 0 else 0 for a, t in zip(abinitio, totals) if t > 0]
        hints_frac = [h/t if t > 0 else 0 for h, t in zip(hints, totals) if t > 0]
        if abinitio_frac:
            print(f"Ab-initio TRAINING (GeneMark-ES self-train + Augustus train [PASA or BUSCO path] "
                  f"+ SNAP fathom/forge) share of PREDICT wall time: "
                  f"median={statistics.median(abinitio_frac)*100:.1f}% "
                  f"mean={100*sum(abinitio_frac)/len(abinitio_frac):.1f}%")
        if hints_frac:
            print(f"RNA-seq hints prep (bam2hints+filterIntronsFindStrand+join_mult_hints) "
                  f"share of PREDICT wall time: median={statistics.median(hints_frac)*100:.1f}% "
                  f"mean={100*sum(hints_frac)/len(hints_frac):.1f}%")

        # Stratify by Augustus evidence mode (PASA vs BUSCO-seeded training path).
        pasa_rows = [r for r in per_run_rows if r[8] == 'pasa']
        busco_rows = [r for r in per_run_rows if r[8] == 'busco']
        unknown_rows = [r for r in per_run_rows if r[8] == 'unknown']
        print(f"\nAugustus evidence-mode split: pasa={len(pasa_rows)} busco={len(busco_rows)} "
              f"unknown/unmatched={len(unknown_rows)} (of n={n})")
        for label, rows_ in (('pasa', pasa_rows), ('busco', busco_rows)):
            if not rows_:
                continue
            aug_secs = [r[6] if label == 'pasa' else r[7] for r in rows_]
            aug_secs_all = [(r[6] + r[7]) for r in rows_]
            totals_ = [r[1] for r in rows_]
            frac = [a / t if t > 0 else 0 for a, t in zip(aug_secs_all, totals_) if t > 0]
            print(f"  [{label}] n={len(rows_)} median Augustus-train time={statistics.median(aug_secs)/60:.1f} min "
                  f"p90={sorted(aug_secs)[int(len(aug_secs)*0.9)]/60:.1f} min "
                  f"median % of that run's PREDICT wall time={statistics.median(frac)*100:.1f}%")


if __name__ == '__main__':
    main()
