#!/usr/bin/env python3
"""Merge per-species TCDB blastp results into a tables-loadable CSV.

TCDB subject headers are '<TC number>|<accession>...' (e.g.
'gnl|TC-DB|Q9X2H0|2.A.1.2.28'); the TC number is the last '|'-delimited
field.
"""

import csv
import gzip
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Merge TCDB blastp results")
    parser.add_argument("blasttabs", nargs="+", help="*.blasttab.gz files")
    parser.add_argument("-o", "--outfile", default="tcdb.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'tc_number',
                    'percent_identity', 'aln_length', 'mismatches',
                    'gap_openings', 'q_start', 'q_end',
                    's_start', 's_end', 'evalue', 'bitscore'])
        for f in sorted(args.blasttabs):
            with gzip.open(f, "rt") as fh:
                for row in csv.reader(fh, delimiter="\t"):
                    if not row:
                        continue
                    protein_id = row[0]
                    tc_number = row[1].split("|")[-1]
                    w.writerow([protein_id.split("_")[0], protein_id, tc_number] + row[2:])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
