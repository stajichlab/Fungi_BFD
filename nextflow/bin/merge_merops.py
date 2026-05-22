#!/usr/bin/env python3
"""Merge per-species MEROPS blasttab results into a tables-loadable CSV."""

import csv
import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge MEROPS blasttab results")
    parser.add_argument("blasttabs", nargs="+", help="*.blasttab.gz files")
    parser.add_argument("-o", "--outfile", default="merops.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'merops_id',
                    'percent_identity', 'aln_length', 'mismatches',
                    'gap_openings', 'q_start', 'q_end',
                    's_start', 's_end', 'evalue', 'bitscore'])
        for f in sorted(args.blasttabs):
            with gzip.open(f, "rt") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if not row:
                        continue
                    w.writerow([row[0].split("_")[0]] + row)

    print(f"Written: {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
