#!/usr/bin/env python3
"""Merge per-species predGPI GFF3 results into a tables-loadable CSV."""

import csv
import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge predGPI GFF3 results")
    parser.add_argument("gff3s", nargs="+", help="*.predgpi.gff3.gz files")
    parser.add_argument("-o", "--outfile", default="predgpi.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'source', 'feature',
                    'start', 'end', 'score', 'strand', 'phase', 'attributes'])
        for f in sorted(args.gff3s):
            with gzip.open(f, "rt") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if not row or row[0].startswith("#"):
                        continue
                    pid = row[0]
                    w.writerow([pid.split("_")[0]] + row)

    print(f"Written: {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
