#!/usr/bin/env python3
"""Merge per-genome run_dbcan CGC (CAZyme Gene Cluster) output into a
tables-loadable CSV. Column layout follows run_dbcan's cgc_standard_out.tsv
(CGC_id, contig, gene coordinates/annotations per cluster member) -- NOT YET
validated against real output (see modules/BFD/CAZY_CGC/main.nf); adjust the
column indices here once a real cgc_standard_out.tsv is available.
"""

import csv
import gzip
import sys
import argparse
from collections import defaultdict


def main():
    parser = argparse.ArgumentParser(description="Merge run_dbcan CGC results")
    parser.add_argument("cgc_tsvs", nargs="+", help="*.cgc.tsv.gz files")
    parser.add_argument("-o", "--outfile", default="cazy_cgc.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'cgc_id', 'contig', 'gene_ids'])
        for f in sorted(args.cgc_tsvs):
            species_prefix = f.split("/")[-1].split(".")[0]
            genes_by_cgc = defaultdict(list)
            contig_by_cgc = {}
            with gzip.open(f, "rt") as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                for row in reader:
                    cgc_id = row.get("CGC_id") or row.get("CGC#")
                    if not cgc_id:
                        continue
                    gene_id = row.get("Gene_ID") or row.get("Protein_ID", "")
                    contig = row.get("Contig_ID") or row.get("Contig", "")
                    genes_by_cgc[cgc_id].append(gene_id)
                    contig_by_cgc[cgc_id] = contig
            for cgc_id, genes in genes_by_cgc.items():
                w.writerow([species_prefix, cgc_id, contig_by_cgc[cgc_id], ";".join(genes)])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
