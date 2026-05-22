#!/usr/bin/env python3
"""Merge per-species SignalP GFF3 results into a tables-loadable CSV."""

import csv
import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge SignalP GFF3 results")
    parser.add_argument("gff3s", nargs="+", help="*.signalp.gff3.gz files")
    parser.add_argument("-o", "--outfile", default="signalp.signal_peptide.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id',
                    'peptide_start', 'peptide_end', 'probability'])
        for f in sorted(args.gff3s):
            with gzip.open(f, "rt") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if not row or row[0].startswith("#"):
                        continue
                    pid = row[0].split(" ")[0]
                    w.writerow([pid.split("_")[0], pid, row[3], row[4], row[5]])

    print(f"Written: {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
