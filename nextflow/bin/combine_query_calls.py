#!/usr/bin/env python3
"""
combine_query_calls.py — Concatenate every per-group *_query_calls.tsv written by
query_ANI.nf's REPORT_QUERY_ANI into one master orphan-classification table.

Row order: same_genus_high_confidence, then closely_related_review, then
no_close_match, then no_alignment; best_ani descending within each tier. This
puts the strongest, most actionable calls (adopt this genus) at the top and the
genomes still needing manual/taxonomic follow-up at the bottom.
"""

import argparse
import csv
import glob
import os


TIER_ORDER = {
    'same_genus_high_confidence': 0,
    'closely_related_review': 1,
    'no_close_match': 2,
    'no_alignment': 3,
}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--calls-dir', required=True)
    ap.add_argument('--output', required=True)
    args = ap.parse_args()

    rows = []
    cols = None
    for path in sorted(glob.glob(os.path.join(args.calls_dir, '*_query_calls.tsv'))):
        with open(path) as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            if cols is None:
                cols = reader.fieldnames
            rows.extend(reader)

    def sort_key(r):
        ani = float(r['best_ani']) if r['best_ani'] else -1.0
        return (TIER_ORDER.get(r['tier'], 9), -ani)

    rows.sort(key=sort_key)

    with open(args.output, 'w', newline='') as out:
        writer = csv.DictWriter(out, fieldnames=cols or [])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} orphan classification rows to {args.output}")


if __name__ == '__main__':
    main()
