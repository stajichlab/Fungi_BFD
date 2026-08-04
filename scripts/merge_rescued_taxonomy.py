#!/usr/bin/env python3
"""
Merge rescued taxonomy data back into samples.csv to produce updated file.

This script takes samples.csv and fills in missing taxonomy fields using data
from the rescue files (taxonomy_rescue.csv and basidiobolus_taxonomy_rescue.csv).
"""

import csv
from pathlib import Path
from typing import Dict, List

def load_rescue_data(rescue_csv: str) -> Dict[str, Dict[str, str]]:
    """Load rescue CSV indexed by ASMID."""
    rescue_data = {}
    try:
        with open(rescue_csv, 'r', newline='', encoding='utf-8') as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row and 'ASMID' in row:
                    asmid = row['ASMID']
                    rescue_data[asmid] = row
    except FileNotFoundError:
        print(f"Warning: {rescue_csv} not found")
    return rescue_data

def merge_rescue_into_samples(
    samples_csv: str,
    taxonomy_rescue_csv: str,
    basidiobolus_rescue_csv: str,
    output_csv: str
) -> None:
    """
    Merge rescued taxonomy data into samples.csv.

    Args:
        samples_csv: Path to original samples.csv
        taxonomy_rescue_csv: Path to taxonomy_rescue.csv
        basidiobolus_rescue_csv: Path to basidiobolus_taxonomy_rescue.csv
        output_csv: Path to save updated samples.csv
    """

    # Load rescue data
    print("Loading rescue data...")
    taxonomy_rescue = load_rescue_data(taxonomy_rescue_csv)
    print(f"  Loaded {len(taxonomy_rescue)} entries from taxonomy_rescue.csv")

    basidiobolus_rescue = load_rescue_data(basidiobolus_rescue_csv)
    print(f"  Loaded {len(basidiobolus_rescue)} entries from basidiobolus_taxonomy_rescue.csv")

    # Merge both rescue datasets
    all_rescue = {**taxonomy_rescue, **basidiobolus_rescue}
    print(f"  Total rescue entries: {len(all_rescue)}")

    # Read and update samples.csv
    print(f"\nProcessing {samples_csv}...")
    updated_rows = []
    header = None
    updates_applied = 0
    already_filled = 0

    with open(samples_csv, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        header = reader.fieldnames

        for row in reader:
            asmid = row['ASMID']

            # Check if this ASMID has empty PHYLUM
            if not row.get('PHYLUM', '').strip():
                # Check if we have rescue data for this ASMID
                if asmid in all_rescue:
                    rescue_row = all_rescue[asmid]

                    # Map rescue columns to samples columns
                    taxonomy_cols = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']

                    for col in taxonomy_cols:
                        if col in rescue_row and rescue_row[col].strip():
                            # Update the row with rescue data
                            # For SPECIES, only update if currently empty to preserve existing data
                            if col == 'SPECIES' and row.get('SPECIES', '').strip():
                                continue
                            row[col] = rescue_row[col]

                    updates_applied += 1
            else:
                already_filled += 1

            updated_rows.append(row)

    # Write updated samples.csv
    print(f"\nWriting updated samples to {output_csv}...")
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=header, lineterminator='\n')
        writer.writeheader()
        writer.writerows(updated_rows)

    print(f"\n✓ Merge complete!")
    print(f"  Total rows processed: {len(updated_rows)}")
    print(f"  Rows updated with rescue data: {updates_applied}")
    print(f"  Rows already with taxonomy: {already_filled}")
    print(f"  Output saved to: {output_csv}")

    # Summary of what was rescued
    print(f"\nRescued taxonomy entries:")
    for asmid in sorted(all_rescue.keys()):
        rescue_row = all_rescue[asmid]
        species = rescue_row.get('SPECIES', '')
        phylum = rescue_row.get('PHYLUM', '')
        if phylum:
            print(f"  {asmid}: {species[:50]} → {phylum}")
        else:
            print(f"  {asmid}: {species[:50]} (partial data)")

if __name__ == '__main__':
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent

    samples_csv = repo_root / 'samples.csv'
    taxonomy_rescue_csv = repo_root / 'data' / 'curation' / 'taxonomy_rescue.csv'
    basidiobolus_rescue_csv = repo_root / 'data' / 'curation' / 'basidiobolus_taxonomy_rescue.csv'
    output_csv = repo_root / 'samples_update_rescue.csv'

    if not samples_csv.exists():
        print(f"Error: {samples_csv} not found")
        exit(1)

    merge_rescue_into_samples(
        str(samples_csv),
        str(taxonomy_rescue_csv),
        str(basidiobolus_rescue_csv),
        str(output_csv)
    )
