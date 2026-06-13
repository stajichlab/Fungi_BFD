#!/usr/bin/env python3
"""Build species.csv.gz — the taxonomy table consumed by the DuckDB build/MCP server.

Projects samples.csv down to one row per genome carrying the identifier and
taxonomy fields the BFD DuckDB build expects: LOCUSTAG, ASMID, the taxonomy ranks
(PHYLUM → SPECIES), and BUSCO_LINEAGE. Restricted to the genomes matched/processed
in a run (via --matched), so a --taxon run yields a clade-restricted species table.

This is intentionally distinct from samples.csv.gz (the run manifest, which keeps
every original column): species.csv.gz is the normalized table the database loads
and joins to asm_stats on ASMID / to the functional tables on LOCUSTAG.
"""

import argparse
import csv
import gzip
import os
import sys

# Columns the DuckDB build and MCP server reference on the `species` table.
REQUIRED_COLUMNS = [
    "ASMID", "LOCUSTAG",
    "PHYLUM", "SUBPHYLUM", "CLASS", "ORDER", "FAMILY", "GENUS", "SPECIES",
    "BUSCO_LINEAGE",
]
# Included when present in samples.csv, otherwise silently skipped.
OPTIONAL_COLUMNS = ["SUBCLASS", "NCBI_TAXONID", "STRAIN"]


def main():
    parser = argparse.ArgumentParser(
        description="Build species.csv.gz (taxonomy table) from samples.csv.")
    parser.add_argument("--samples", required=True,
                        help="full samples CSV (with header)")
    parser.add_argument("--matched",
                        help="file with one matched key per line; if omitted, keep all rows")
    parser.add_argument("--key", default="LOCUSTAG",
                        help="samples column to match against --matched values [LOCUSTAG]")
    parser.add_argument("-o", "--outfile", default="tables/All_Taxa/species.csv.gz",
                        help="output CSV; gzipped when name ends in .gz [tables/All_Taxa/species.csv.gz]")
    args = parser.parse_args()

    keep = None
    if args.matched:
        with open(args.matched) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    opener = gzip.open if args.outfile.endswith(".gz") else open
    n = 0
    with open(args.samples, newline="") as ifh, opener(args.outfile, "wt", newline="") as ofh:
        reader = csv.DictReader(ifh)
        fields = reader.fieldnames or []

        missing = [c for c in REQUIRED_COLUMNS if c not in fields]
        if missing:
            print(f"ERROR: {args.samples} missing required columns {missing} "
                  f"(have: {fields})", file=sys.stderr)
            sys.exit(1)
        if args.key not in fields:
            print(f"ERROR: key column '{args.key}' not in {args.samples}", file=sys.stderr)
            sys.exit(1)

        out_columns = REQUIRED_COLUMNS + [c for c in OPTIONAL_COLUMNS if c in fields]
        writer = csv.DictWriter(ofh, fieldnames=out_columns, extrasaction="ignore")
        writer.writeheader()
        for row in reader:
            if keep is not None and (row.get(args.key) or "").strip() not in keep:
                continue
            writer.writerow({c: (row.get(c) or "") for c in out_columns})
            n += 1

    print(f"Wrote {n} rows to {args.outfile}")


if __name__ == "__main__":
    main()
