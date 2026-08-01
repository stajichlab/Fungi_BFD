#!/usr/bin/env python3
"""
Validate that genome_functional.nf produced the expected tables Parquet outputs.
Used by run_test.sh after a -stub-run. Per T-014 §D.1/#27, tables/ outputs are
Parquet (one unpartitioned file per type), not CSV.gz.
"""

import argparse
import sys
from pathlib import Path

import duckdb

EXPECTED_BIGQUERY = {
    "pfam.parquet": [
        "protein_id", "pfam_id", "pfam_acc",
    ],
    "cazy.overview.parquet": [
        "species_prefix", "protein_id", "EC", "cazyme_fam",
    ],
    "cazy.cazymes_hmm.parquet": [
        "species_prefix", "HMM_id", "protein_id",
    ],
    "merops.parquet": [
        "species_prefix", "protein_id", "merops_id",
    ],
    "signalp.signal_peptide.parquet": [
        "species_prefix", "protein_id", "peptide_start",
    ],
    "tmhmm.parquet": [
        "species_prefix", "protein_id", "PredHel",
    ],
    "targetP.parquet": [
        "species_prefix", "protein_id", "prediction",
    ],
    "idp.parquet": [
        "protein_id",
    ],
    "idp_summary.parquet": [
        "protein_id",
    ],
    "wolfpsort.parquet": [
        "species_prefix", "protein_id", "localization",
    ],
    "predgpi.parquet": [
        "species_prefix", "protein_id", "feature",
    ],
}


def check_parquet(path: Path, required_cols: list[str]) -> list[str]:
    """Return list of error strings (empty = pass)."""
    errors = []
    if not path.exists():
        return [f"MISSING: {path}"]
    try:
        header = duckdb.sql(f"SELECT * FROM read_parquet('{path}') LIMIT 0").columns
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
        errs = check_parquet(bq / fname, cols)
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
