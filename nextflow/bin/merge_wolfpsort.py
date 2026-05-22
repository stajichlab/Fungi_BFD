#!/usr/bin/env python3
"""Merge per-species WoLF PSORT results into a tables-loadable CSV (best hit only)."""

import csv
import gzip
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(description="Merge WoLF PSORT results")
    parser.add_argument("results", nargs="+",
                        help="*.wolfpsort.results.txt.gz files")
    parser.add_argument("-o", "--outfile", default="wolfpsort.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'localization', 'score'])
        for f in sorted(args.results):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    line = line.strip()
                    if not line:
                        continue
                    pid, resultstr = line.split(" ", 1)
                    first = resultstr.split(", ")[0]
                    code, score = first.split(" ")
                    w.writerow([pid.split("_")[0], pid, code, score])

    print(f"Written: {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
