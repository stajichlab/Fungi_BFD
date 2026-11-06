#!/usr/bin/env python3
"""check_table_coverage.py -- refuse a merged table that covers far fewer
genomes than the table it replaces.

Usage:
    check_table_coverage.py --new pfam.parquet --old /abs/tables/pfam.parquet \\
        [--min-frac 0.9] [--key COLUMN]

The genome key is the first column found among species_prefix, LOCUSTAG,
ASMID, genome (or --key). Exit 0 when the old table is missing, has no key
column, or the new table has at least min_frac x the old table's distinct
genomes. Exit 1 otherwise, with a message naming both counts.

Why: on 2026-11-05 a merge-only run replaced targetP.parquet (728 genomes)
with a 127-genome table, because TargetP had been set aside and only part of
its result files were on disk. Row-count checks miss this for tables with many
rows per genome; distinct genomes catch it.

Uses the duckdb CLI, like publish_table.py.
"""
import argparse
import os
import subprocess
import sys

KEYS = ["species_prefix", "LOCUSTAG", "ASMID", "genome"]


def duck(sql):
    out = subprocess.run(["duckdb", "-csv", "-noheader", "-c", sql],
                         capture_output=True, text=True, check=True)
    return out.stdout.strip()


def columns(path):
    return duck(f"SELECT column_name FROM (DESCRIBE SELECT * FROM read_parquet('{path}'));").split()


def distinct(path, key):
    return int(duck(f'SELECT COUNT(DISTINCT "{key}") FROM read_parquet(\'{path}\');'))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--new", required=True)
    ap.add_argument("--old", required=True)
    ap.add_argument("--min-frac", type=float, default=0.9)
    ap.add_argument("--key", default=None)
    a = ap.parse_args()
    if not os.path.isfile(a.old):
        print(f"coverage: no previous {os.path.basename(a.old)}; nothing to compare")
        return 0
    old_cols, new_cols = columns(a.old), columns(a.new)
    key = a.key or next((k for k in KEYS if k in old_cols and k in new_cols), None)
    if key is None:
        print(f"coverage: {os.path.basename(a.new)} has no genome key column; check skipped")
        return 0
    n_old, n_new = distinct(a.old, key), distinct(a.new, key)
    print(f"coverage: {os.path.basename(a.new)} {key}: new={n_new} old={n_old} (min {a.min_frac:.2f} x old)")
    if n_old and n_new < a.min_frac * n_old:
        print(f"ERROR: {os.path.basename(a.new)} covers {n_new} genomes, the table it replaces covers "
              f"{n_old}. Refusing to publish. If the shrink is deliberate, rerun with "
              f"--allow_shrink true (params file).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
