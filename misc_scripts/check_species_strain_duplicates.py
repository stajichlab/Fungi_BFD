#!/usr/bin/env python3
"""
Check samples.csv for duplicate SPECIES+STRAIN combinations.

Outputs each collision group with line numbers, ASMID, BIOPROJECT, and
assembly stats (N50, total length) pulled from results/genome_stats/asm_reports/.

Also writes samples.dupremoved.csv alongside the input file, keeping only the
highest-N50 assembly for each duplicate group.

Usage:
    python check_species_strain_duplicates.py [samples.csv] [asm_reports_dir]
"""

import csv
import os
import re
import sys
from collections import defaultdict


def load_asm_stats(asm_reports_dir: str, asmid: str) -> dict:
    """Return N50 and total_length from {asmid}.stats.txt, or None values if missing."""
    stats_file = os.path.join(asm_reports_dir, f"{asmid}.stats.txt")
    result = {"n50": None, "total_length": None}
    if not os.path.isfile(stats_file):
        return result
    with open(stats_file) as fh:
        for line in fh:
            m = re.match(r'\s+N50\s*=\s*(\d+)', line)
            if m:
                result["n50"] = int(m.group(1))
            m = re.match(r'\s+TOTAL LENGTH\s*=\s*(\d+)', line)
            if m:
                result["total_length"] = int(m.group(1))
    return result


def fmt_bp(n):
    if n is None:
        return "N/A"
    if n >= 1_000_000:
        return f"{n/1_000_000:.2f} Mb"
    if n >= 1_000:
        return f"{n/1_000:.1f} kb"
    return str(n)


def main():
    infile      = sys.argv[1] if len(sys.argv) > 1 else "samples.csv"
    asm_dir_arg = sys.argv[2] if len(sys.argv) > 2 else None

    if asm_dir_arg:
        asm_reports_dir = asm_dir_arg
    else:
        base = os.path.dirname(os.path.abspath(infile))
        asm_reports_dir = os.path.join(base, "results", "genome_stats", "asm_reports")

    asm_dir_exists = os.path.isdir(asm_reports_dir)
    if not asm_dir_exists:
        print(f"Warning: asm_reports dir not found: {asm_reports_dir}", file=sys.stderr)

    # Read all rows, preserving order
    all_rows = []
    fieldnames = None
    with open(infile, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames
        for lineno, row in enumerate(reader, start=2):
            all_rows.append((lineno, row))

    # Group by SPECIES+STRAIN
    groups: dict[str, list] = defaultdict(list)
    for lineno, row in all_rows:
        species = (row.get("SPECIES") or "").strip()
        strain  = (row.get("STRAIN")  or "").strip()
        key = f"{species}\t{strain}"
        asmid = (row.get("ASMID") or "").strip()
        stats = load_asm_stats(asm_reports_dir, asmid) if asm_dir_exists else {"n50": None, "total_length": None}
        groups[key].append({
            "lineno":     lineno,
            "asmid":      asmid,
            "locustag":   (row.get("LOCUSTAG")   or "").strip(),
            "bioproject": (row.get("BIOPROJECT") or "").strip(),
            "species":    species,
            "strain":     strain,
            "n50":        stats["n50"],
            "total_length": stats["total_length"],
            "row":        row,
        })

    duplicates = {k: v for k, v in groups.items() if len(v) > 1}

    if not duplicates:
        print("No duplicate SPECIES+STRAIN combinations found.")
        return

    # Determine which LOCUSTAGs to drop (all but the best N50 in each group)
    drop_locustags: set[str] = set()
    for key, entries in duplicates.items():
        # Sort by N50 descending; treat None as -1 so entries with stats rank above missing
        ranked = sorted(entries, key=lambda e: e["n50"] if e["n50"] is not None else -1, reverse=True)
        for loser in ranked[1:]:
            drop_locustags.add(loser["locustag"])

    # Report
    print(f"Found {len(duplicates)} duplicate group(s):\n")
    for key, entries in sorted(duplicates.items(), key=lambda x: x[1][0]["lineno"]):
        species, strain = key.split("\t", 1)
        label = f"{species} / {strain}" if strain else species
        ranked = sorted(entries, key=lambda e: e["n50"] if e["n50"] is not None else -1, reverse=True)
        print(f"  SPECIES+STRAIN: {label!r}  ({len(entries)} entries)")
        for e in ranked:
            kept = e["locustag"] not in drop_locustags
            flag = "  [KEEP]" if kept else "  [DROP]"
            print(f"    line {e['lineno']:5d}  ASMID={e['asmid']:<45s}  "
                  f"LOCUSTAG={e['locustag']}  BIOPROJECT={e['bioproject']:<15s}  "
                  f"N50={fmt_bp(e['n50']):<12s}  total={fmt_bp(e['total_length'])}{flag}")
        print()

    print(f"Total duplicate groups : {len(duplicates)}")
    print(f"Total affected entries : {sum(len(v) for v in duplicates.values())}")
    print(f"Entries to drop        : {len(drop_locustags)}")

    # Write filtered output
    outfile = os.path.join(
        os.path.dirname(os.path.abspath(infile)),
        "samples.dupremoved.csv"
    )
    kept = 0
    with open(outfile, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for _lineno, row in all_rows:
            locustag = (row.get("LOCUSTAG") or "").strip()
            if locustag not in drop_locustags:
                writer.writerow(row)
                kept += 1

    print(f"\nWrote {kept} rows to {outfile}  (dropped {len(drop_locustags)})")


if __name__ == "__main__":
    main()
