#!/usr/bin/env python3
"""Merge per-species DeepTMHMM TMRs.gff3 output into a tables-loadable CSV.

DeepTMHMM's TMRs.gff3 is a simplified 4-column format (not real GFF3):
  protein_id  region_type(TMhelix|inside|outside|signal)  start  end
plus '#' comment lines (e.g. summary "# proteinID Length TM Count").
"""

import csv
import gzip
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Merge DeepTMHMM TMRs.gff3 results")
    parser.add_argument("gff3s", nargs="+", help="*.deeptmhmm.gff3.gz files")
    parser.add_argument("-o", "--outfile", default="deeptmhmm.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'feature', 'start', 'end'])
        for f in sorted(args.gff3s):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 4:
                        continue
                    protein_id, feature, start, end = fields[0], fields[1], fields[2], fields[3]
                    if feature != "TMhelix":
                        continue
                    w.writerow([protein_id.split("_")[0], protein_id, feature, start, end])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
