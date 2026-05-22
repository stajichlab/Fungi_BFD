#!/usr/bin/env python3
"""Merge per-species TMHMM short results into a tables-loadable CSV."""

import csv
import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge TMHMM short results")
    parser.add_argument("tsvs", nargs="+", help="*.tmhmm_short.tsv.gz files")
    parser.add_argument("-o", "--outfile", default="tmhmm.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id',
                    'len', 'ExpAA', 'First60', 'PredHel', 'Topology'])
        for f in sorted(args.tsvs):
            with gzip.open(f, "rt") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if not row or row[0].startswith("#"):
                        continue
                    pid = row[0]
                    parsed = [pid.split("_")[0], pid]
                    for field in row[1:]:
                        parsed.append(field.split("=")[1])
                    if parsed[-2] == "0":   # skip proteins with 0 TM helices
                        continue
                    w.writerow(parsed)

    print(f"Written: {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
