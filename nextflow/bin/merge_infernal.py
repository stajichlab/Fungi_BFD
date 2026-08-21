#!/usr/bin/env python3
"""Merge per-genome Infernal cmscan --tblout Rfam hits into a tables-loadable CSV.

cmscan --tblout columns (space-delimited, fixed width):
  target_name target_acc query_name query_acc mdl mdl_from mdl_to
  seq_from seq_to strand trunc pass gc bias score evalue inc description...
"""

import csv
import gzip
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Merge Infernal cmscan --tblout Rfam results")
    parser.add_argument("tblouts", nargs="+", help="*.rfam.tblout.gz files")
    parser.add_argument("-o", "--outfile", default="infernal_rfam.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'target_name', 'target_acc', 'query_name',
                    'contig', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to',
                    'strand', 'evalue', 'score'])
        for f in sorted(args.tblouts):
            species_prefix = f.split("/")[-1].split(".")[0]
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if line.startswith("#") or not line.strip():
                        continue
                    fields = line.split(None, 17)
                    if len(fields) < 16:
                        continue
                    target_name, target_acc, query_name = fields[0], fields[1], fields[2]
                    mdl_from, mdl_to = fields[5], fields[6]
                    seq_from, seq_to, strand = fields[7], fields[8], fields[9]
                    score, evalue = fields[14], fields[15]
                    w.writerow([species_prefix, target_name, target_acc, query_name,
                                query_name, mdl_from, mdl_to, seq_from, seq_to,
                                strand, evalue, score])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
