#!/usr/bin/env python3
"""Merge per-species phobius.pl -short output into a tables-loadable CSV.

-short format columns: SEQENCE_ID  TM  SP  Prediction
"""

import csv
import gzip
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Merge Phobius -short results")
    parser.add_argument("short_txts", nargs="+", help="*.phobius.short.txt.gz files")
    parser.add_argument("-o", "--outfile", default="phobius.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'tm_count', 'sp_predicted', 'prediction'])
        for f in sorted(args.short_txts):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    fields = line.split()
                    if len(fields) < 4 or fields[0] == "SEQENCE_ID":
                        continue
                    protein_id, tm_count, sp, prediction = fields[0], fields[1], fields[2], fields[3]
                    if tm_count == "0" and sp == "0":
                        continue
                    w.writerow([protein_id.split("_")[0], protein_id, tm_count,
                                "Y" if sp != "0" else "N", prediction])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
