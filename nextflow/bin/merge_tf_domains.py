#!/usr/bin/env python3
"""Classify transcription factors by filtering RUN_PFAM domtblout hits
against a curated list of fungal TF-associated Pfam accessions
(assets/fungal_tf_pfam_domains.csv). No FTFD download needed -- ftfd.snu.ac.kr
does not resolve (confirmed 2026-08-20), so this derives TF calls from
Pfam output already computed by RUN_PFAM instead of a second external tool.
"""

import csv
import gzip
import sys
import argparse


def load_tf_domains(path):
    domains = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            domains[row["pfam_acc"]] = (row["pfam_name"], row["tf_family"])
    return domains


def main():
    parser = argparse.ArgumentParser(description="Classify TF domains from Pfam domtblout hits")
    parser.add_argument("domtbls", nargs="+", help="*.domtblout.gz files (hmmsearch --domtbl output)")
    parser.add_argument("--tf-domains", required=True, help="assets/fungal_tf_pfam_domains.csv")
    parser.add_argument("-o", "--outfile", default="tf_inventory.csv")
    args = parser.parse_args()

    tf_domains = load_tf_domains(args.tf_domains)

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'protein_id', 'pfam_acc', 'pfam_name',
                    'tf_family', 'domain_i_evalue', 'domain_score'])
        for f in sorted(args.domtbls):
            with gzip.open(f, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        continue
                    fields = line.strip().split()
                    if len(fields) < 21:
                        continue
                    protein_id = fields[0]
                    hmm_acc = fields[4].split(".")[0]
                    if hmm_acc not in tf_domains:
                        continue
                    pfam_name, tf_family = tf_domains[hmm_acc]
                    domain_i_evalue = fields[12]
                    domain_score = fields[13]
                    w.writerow([protein_id.split("_")[0], protein_id, hmm_acc,
                                pfam_name, tf_family, domain_i_evalue, domain_score])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
