#!/usr/bin/env python3
"""Audit how SPECIES/STRAIN are transformed by the sanitizer.

Reads the raw NCBI accession + taxonomy tables and writes a record with columns:
    ASMID, ORIG SPECIES, MODIFIED SPECIES, ORIG STRAIN, MODIFIED STRAIN
so the original NCBI assignment and the cleaned value used in samples.csv can be
compared row by row.

Uses the shared sanitizer (scripts/sample_sanitize.py) so the audit always
reflects exactly what create_samples_file.py will produce.
"""

import argparse
import csv
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sample_sanitize import clean_species, clean_strain, clean_asmid, backfill_strain

DEF_NCBI = "../../1KFG/2026/NCBI_fungi"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--accessions", default=f"{DEF_NCBI}/ncbi_accessions.csv",
                    help="ncbi_accessions.csv (ACCESSION,SPECIES,STRAIN,...,ASM_NAME)")
    ap.add_argument("--taxonomy", default=f"{DEF_NCBI}/ncbi_accessions_taxonomy.csv",
                    help="ncbi_accessions_taxonomy.csv (ASM_ACCESSION,...,SPECIES)")
    ap.add_argument("--outfile", default="data/curation/strain_species_audit.csv")
    ap.add_argument("--changed-only", action="store_true",
                    help="only emit rows where species or strain changed")
    args = ap.parse_args()

    # binomial SPECIES per assembly, deduped (the taxonomy file may carry
    # duplicate lines; first occurrence wins).
    binomial = {}
    with open(args.taxonomy, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            asm = row["ASM_ACCESSION"]
            binomial.setdefault(asm, (row.get("SPECIES", "") or "").strip())

    rows = []
    n_sp = n_st = 0
    with open(args.accessions, newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            asm_base = f"{row['ACCESSION']}_{row['ASM_NAME']}"
            asmid = clean_asmid(asm_base)
            species_in = (row.get("SPECIES", "") or "").strip()
            orig_species = binomial.get(asm_base, species_in)
            orig_strain = (row.get("STRAIN", "") or "").strip()

            mod_species = clean_species(orig_species)
            mod_strain = backfill_strain(
                clean_strain(orig_strain),
                clean_species(species_in),
                mod_species,
            )

            changed_sp = mod_species != orig_species
            changed_st = mod_strain != orig_strain
            if changed_sp:
                n_sp += 1
            if changed_st:
                n_st += 1
            if args.changed_only and not (changed_sp or changed_st):
                continue
            rows.append([asmid, orig_species, mod_species, orig_strain, mod_strain])

    os.makedirs(os.path.dirname(args.outfile) or ".", exist_ok=True)
    with open(args.outfile, "w", newline="") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(["ASMID", "ORIG SPECIES", "MODIFIED SPECIES",
                    "ORIG STRAIN", "MODIFIED STRAIN"])
        w.writerows(rows)

    print(f"wrote {len(rows)} rows -> {args.outfile}")
    print(f"  species changed: {n_sp}")
    print(f"  strain  changed: {n_st}")


if __name__ == "__main__":
    main()
