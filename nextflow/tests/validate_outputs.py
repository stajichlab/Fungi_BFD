#!/usr/bin/env python3
"""
Validate that genome_functional.nf produced the expected tables CSV.gz outputs.
Used by run_test.sh after a -stub-run.
"""

import argparse
import csv
import gzip
import sys
from pathlib import Path

EXPECTED_BIGQUERY = {
    "pfam.csv.gz": [
        "protein_id", "hmm_id", "hmm_acc",
    ],
    "cazy.overview.csv.gz": [
        "species_prefix", "protein_id", "EC", "cazyme_fam",
    ],
    "cazy.cazymes_hmm.csv.gz": [
        "species_prefix", "HMM_id", "protein_id",
    ],
    "merops.csv.gz": [
        "species_prefix", "protein_id", "merops_id",
    ],
    "signalp.signal_peptide.csv.gz": [
        "species_prefix", "protein_id", "peptide_start",
    ],
    "tmhmm.csv.gz": [
        "species_prefix", "protein_id", "PredHel",
    ],
    "targetP.csv.gz": [
        "species_prefix", "protein_id", "prediction",
    ],
    "idp.csv.gz": [
        "protein_id",
    ],
    "idp_summary.csv.gz": [
        "protein_id",
    ],
    "wolfpsort.csv.gz": [
        "species_prefix", "protein_id", "localization",
    ],
    "predgpi.csv.gz": [
        "species_prefix", "protein_id", "feature",
    ],
}


def check_csv_gz(path: Path, required_cols: list[str]) -> list[str]:
    """Return list of error strings (empty = pass)."""
    errors = []
    if not path.exists():
        return [f"MISSING: {path}"]
    try:
        with gzip.open(path, "rt") as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames or []
            for col in required_cols:
                if col not in header:
                    errors.append(f"{path.name}: missing column '{col}' (header={header})")
    except Exception as exc:
        errors.append(f"{path.name}: could not read — {exc}")
    return errors


def main():
    parser = argparse.ArgumentParser(
        description="Validate genome_functional.nf stub-run outputs"
    )
    parser.add_argument("--tables", required=True,
                        help="Path to tables output directory")
    parser.add_argument("--outdir", required=True,
                        help="Path to results/function output directory (spot-checked)")
    args = parser.parse_args()

    bq = Path(args.tables)
    errors = []

    print(f"Checking tables outputs in: {bq}")
    for fname, cols in EXPECTED_BIGQUERY.items():
        errs = check_csv_gz(bq / fname, cols)
        if errs:
            errors.extend(errs)
        else:
            print(f"  OK  {fname}")

    if errors:
        print("\nFAILURES:")
        for e in errors:
            print(f"  {e}")
        sys.exit(1)
    else:
        print(f"\nAll {len(EXPECTED_BIGQUERY)} tables files present and valid.")


if __name__ == "__main__":
    main()
