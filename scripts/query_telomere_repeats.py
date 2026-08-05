#!/usr/bin/env python3
"""Query telomere repeat data from db/BFD.duckdb, filtered by taxon or ID.

Joins `species` to `telomere_summary` (one row per genome: monomer, scaffold
counts, both-ends count) and, with --tracts, also to `telomere_tracts` (one
row per tract: scaffold, end, strand, sequence) for sequence/length-level
mining.

Examples:
    # Per-genome telomere monomer summary for a genus
    query_telomere_repeats.py --genus Coccidioides

    # Per-tract detail (actual repeat sequence, per scaffold end) for a genus
    query_telomere_repeats.py --genus Coccidioides --tracts

    # Any other taxonomic rank, a single species, or a single assembly
    query_telomere_repeats.py --family Onygenaceae
    query_telomere_repeats.py --species "Coccidioides immitis"
    query_telomere_repeats.py --asmid GCA_000149895.1_ASM14989v1

    # Write to CSV instead of printing
    query_telomere_repeats.py --genus Coccidioides --csv out.csv
"""

import argparse
import csv
import sys

import duckdb

RANK_COLUMNS = {
    "phylum": "PHYLUM",
    "subphylum": "SUBPHYLUM",
    "class": "CLASS",
    "subclass": "SUBCLASS",
    "order": '"ORDER"',
    "family": "FAMILY",
    "genus": "GENUS",
    "species": "SPECIES",
}

SUMMARY_QUERY = """
SELECT
    sp.GENUS, sp.SPECIES, sp.STRAIN, sp.LOCUSTAG, sp.ASMID,
    ts.telomere_scaffolds,
    ts.telomere_scaffolds_both_ends,
    ts.telomere_tracts,
    ts.telomere_top_monomer,
    ts.telomere_monomers
FROM species sp
JOIN telomere_summary ts USING(ASMID)
WHERE {where}
ORDER BY sp.SPECIES, sp.STRAIN
"""

TRACTS_QUERY = """
SELECT
    sp.GENUS, sp.SPECIES, sp.STRAIN, sp.LOCUSTAG, sp.ASMID,
    tt.scaffold, tt.end_type, tt.strand, tt.terminal,
    tt.monomer, tt.repeat_count, tt.tract_length, tt.tract_seq
FROM species sp
JOIN telomere_tracts tt USING(ASMID)
WHERE {where}
ORDER BY sp.SPECIES, sp.STRAIN, tt.scaffold, tt.end_type
"""


def build_filter(args):
    """Return (where_clause, params) from whichever filter flag was given."""
    if args.asmid:
        return "sp.ASMID = ?", [args.asmid]
    if args.locustag:
        return "sp.LOCUSTAG = ?", [args.locustag]
    for rank, column in RANK_COLUMNS.items():
        value = getattr(args, rank)
        if value:
            return f"sp.{column} = ?", [value]
    raise ValueError("no filter given")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--db", default="db/BFD.duckdb", help="path to BFD.duckdb (default: db/BFD.duckdb)")
    filt = p.add_mutually_exclusive_group(required=True)
    filt.add_argument("--asmid")
    filt.add_argument("--locustag")
    for rank in RANK_COLUMNS:
        filt.add_argument(f"--{rank}")
    p.add_argument("--tracts", action="store_true",
                    help="show per-tract detail (scaffold/end/sequence) instead of per-genome summary")
    p.add_argument("--csv", metavar="PATH", help="write results to PATH instead of printing a table")
    args = p.parse_args()

    where, params = build_filter(args)
    query = (TRACTS_QUERY if args.tracts else SUMMARY_QUERY).format(where=where)

    con = duckdb.connect(args.db, read_only=True)
    result = con.execute(query, params)
    columns = [d[0] for d in result.description]
    rows = result.fetchall()

    if not rows:
        print(f"No telomere data found for the given filter.", file=sys.stderr)
        sys.exit(1)

    if args.csv:
        with open(args.csv, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(columns)
            writer.writerows(rows)
        print(f"Wrote {len(rows)} rows to {args.csv}")
    else:
        widths = [max(len(str(c)), *(len(str(r[i])) for r in rows)) for i, c in enumerate(columns)]
        fmt = "  ".join(f"{{:<{w}}}" for w in widths)
        print(fmt.format(*columns))
        print(fmt.format(*("-" * w for w in widths)))
        for r in rows:
            print(fmt.format(*(str(v) if v is not None else "" for v in r)))
        print(f"\n{len(rows)} rows")


if __name__ == "__main__":
    main()
