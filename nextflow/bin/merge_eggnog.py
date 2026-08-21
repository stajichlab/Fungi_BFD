#!/usr/bin/env python3
"""Merge per-species eggnog-mapper .emapper.annotations into a tables-loadable CSV.

Column order matches emapper.py's standard `--output` annotations file:
 0 query          1 seed_ortholog  2 evalue        3 score
 4 eggNOG_OGs     5 max_annot_lvl  6 COG_category  7 Description
 8 Preferred_name 9 GOs            10 EC           11 KEGG_ko
 12 KEGG_Pathway  13 KEGG_Module   14 KEGG_Reaction 15 KEGG_rclass
 16 BRITE         17 KEGG_TC       18 CAZy          19 BiGG_Reaction
 20 PFAMs
"""

import csv
import gzip
import sys
import argparse


def main():
    parser = argparse.ArgumentParser(description="Merge eggnog-mapper annotation results")
    parser.add_argument("annotations", nargs="+", help="*.emapper.annotations.gz files")
    parser.add_argument("-o", "--outfile", default="eggnog.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'seed_ortholog', 'evalue', 'score',
                    'cog_category', 'description', 'preferred_name', 'go_terms', 'ec',
                    'kegg_ko', 'kegg_pathway', 'cazy', 'pfams'])
        for f in sorted(args.annotations):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if not line.strip() or line.startswith("#"):
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 21:
                        continue
                    protein_id = fields[0]
                    w.writerow([
                        protein_id.split("_")[0], protein_id,
                        fields[1], fields[2], fields[3], fields[6], fields[7],
                        fields[8], fields[9], fields[10], fields[11], fields[12],
                        fields[18], fields[20],
                    ])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
