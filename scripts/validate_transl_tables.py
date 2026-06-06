#!/usr/bin/env python3
"""
Validate that the translation table used by funannotate predict matches
the TRANSL_TABLE in samples.csv for every completed genome annotation.

Reads: samples.csv (LOCUSTAG → TRANSL_TABLE)
Scans: genome_annotation/*/predict_results/*.stats.json
       - extracts --name <LOCUSTAG> and --table <N> from the stored command
       - absent --table means funannotate defaulted to 1

Outputs a TSV report to stdout with one row per annotation checked, and
prints a summary of mismatches/missing at the end.

Usage:
    python scripts/validate_transl_tables.py [--annot-dir genome_annotation] [--samples samples.csv]
"""

import argparse
import csv
import json
import re
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annot-dir", default="genome_annotation",
                   help="Root genome_annotation directory (default: genome_annotation)")
    p.add_argument("--samples", default="samples.csv",
                   help="samples.csv path (default: samples.csv)")
    return p.parse_args()


def load_samples(csv_path):
    """Return dict locustag -> expected_table (int)."""
    expected = {}
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            lt = row.get("LOCUSTAG", "").strip()
            tt = row.get("TRANSL_TABLE", "1").strip()
            if not lt:
                continue
            try:
                expected[lt] = int(tt)
            except ValueError:
                # malformed — record as-is so we can flag it
                expected[lt] = tt
    return expected


def parse_stats_json(stats_path):
    """
    Return (locustag, table_used) extracted from stats.json command field.
    table_used defaults to 1 if --table is absent from the command.
    Returns (None, None) if the file is unparseable or has no command.
    """
    try:
        with open(stats_path) as fh:
            data = json.load(fh)
    except (json.JSONDecodeError, OSError):
        return None, None

    command = data.get("command", "")
    if not command:
        return None, None

    # --name <LOCUSTAG>  (8-char hex string used as the locus tag prefix)
    name_m = re.search(r"--name\s+(\S+)", command)
    locustag = name_m.group(1) if name_m else None

    # --table <N>  (absent → funannotate default = 1)
    table_m = re.search(r"--table\s+(\d+)", command)
    table_used = int(table_m.group(1)) if table_m else 1

    return locustag, table_used


def main():
    args = parse_args()
    annot_dir = Path(args.annot_dir)
    samples_path = Path(args.samples)

    if not samples_path.exists():
        sys.exit(f"ERROR: samples file not found: {samples_path}")
    if not annot_dir.is_dir():
        sys.exit(f"ERROR: annotation directory not found: {annot_dir}")

    expected = load_samples(samples_path)

    print("\t".join(["folder", "locustag", "expected_table", "used_table", "status"]))

    n_ok = n_mismatch = n_no_locustag = n_not_in_csv = n_no_stats = 0

    for stats_json in sorted(annot_dir.glob("*/predict_results/*.stats.json")):
        folder = stats_json.parts[-3]  # genome_annotation/<folder>/predict_results/...
        locustag, table_used = parse_stats_json(stats_json)

        if locustag is None:
            status = "NO_COMMAND"
            n_no_locustag += 1
            print("\t".join([folder, "?", "?", "?", status]))
            continue

        if locustag not in expected:
            status = "NOT_IN_CSV"
            n_not_in_csv += 1
            print("\t".join([folder, locustag, "?", str(table_used), status]))
            continue

        exp = expected[locustag]
        if table_used == exp:
            status = "OK"
            n_ok += 1
        else:
            status = "MISMATCH"
            n_mismatch += 1

        print("\t".join([folder, locustag, str(exp), str(table_used), status]))

    # Locustags in CSV that have no annotation folder at all (not necessarily an error)
    annotated_locustags = set()
    for stats_json in annot_dir.glob("*/predict_results/*.stats.json"):
        lt, _ = parse_stats_json(stats_json)
        if lt:
            annotated_locustags.add(lt)
    missing_annot = {lt for lt in expected if lt not in annotated_locustags}

    print(f"\n# Summary", file=sys.stderr)
    print(f"# Checked:        {n_ok + n_mismatch + n_no_locustag + n_not_in_csv}", file=sys.stderr)
    print(f"# OK:             {n_ok}", file=sys.stderr)
    print(f"# MISMATCH:       {n_mismatch}", file=sys.stderr)
    print(f"# NOT_IN_CSV:     {n_not_in_csv}  (locustag in annotation but not in samples.csv)", file=sys.stderr)
    print(f"# NO_COMMAND:     {n_no_locustag}  (stats.json present but no command field)", file=sys.stderr)
    print(f"# No annotation:  {len(missing_annot)}  (in samples.csv but no predict_results yet)", file=sys.stderr)

    if n_mismatch:
        print(f"\n# MISMATCHES (filter with: grep MISMATCH):", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
