#!/usr/bin/env python3
"""Merge per-species CAZy results into tables-loadable CSV files."""

import csv
import gzip
import sys

def merge_overview(files, outfile):
    with open(outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'EC', 'cazyme_fam',
                    'sub_fam', 'diamond_fam', 'substrate', 'toolcount'])
        for f in sorted(files):
            with gzip.open(f, "rt") as fh:
                reader = csv.reader(fh, delimiter="\t")
                next(reader)  # skip header
                for row in reader:
                    w.writerow([row[0].split("_")[0]] + row)

def merge_hmm(files, outfile):
    with open(outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'HMM_id', 'profile_length', 'protein_id',
                    'protein_length', 'evalue', 'q_start', 'q_end',
                    's_start', 's_end', 'coverage'])
        for f in sorted(files):
            with gzip.open(f, "rt") as fh:
                reader = csv.reader(fh, delimiter="\t")
                next(reader)  # skip header
                for row in reader:
                    w.writerow([row[2].split("_")[0]] + row)

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Merge CAZy per-species results")
    parser.add_argument("--overviews", nargs="+", required=True,
                        help="overview.tsv.gz files")
    parser.add_argument("--cazymes", nargs="+", required=True,
                        help="cazymes.tsv.gz files")
    parser.add_argument("--out-overview", default="cazy.overview.csv")
    parser.add_argument("--out-hmm", default="cazy.cazymes_hmm.csv")
    args = parser.parse_args()

    merge_overview(args.overviews, args.out_overview)
    merge_hmm(args.cazymes, args.out_hmm)
    print(f"Written: {args.out_overview}, {args.out_hmm}", file=sys.stderr)
