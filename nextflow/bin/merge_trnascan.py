#!/usr/bin/env python3
"""Merge per-genome tRNAscan-SE (overlap-filtered) GFF3 into a tables-loadable CSV."""

import csv
import gzip
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Merge tRNAscan-SE no-overlaps GFF3 results")
    parser.add_argument("gff3s", nargs="+", help="*.trnascan.no-overlaps.gff3.gz files")
    parser.add_argument("-o", "--outfile", default="trnascan.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'contig', 'start', 'end', 'strand', 'product', 'note'])
        for f in sorted(args.gff3s):
            species_prefix = f.split("/")[-1].split(".")[0]
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9 or cols[2] != "tRNA":
                        continue
                    attrs = dict(kv.split("=", 1) for kv in cols[8].split(";") if "=" in kv)
                    w.writerow([species_prefix, cols[0], cols[3], cols[4], cols[6],
                                attrs.get("product", ""), attrs.get("note", "")])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
