#!/usr/bin/env python3
"""
Rescue missing taxonomy data from backup_samples.csv for entries with empty PHYLUM in samples.csv.

This script identifies rows in samples.csv where the PHYLUM column is empty,
then attempts to find matching entries in backup_samples.csv that have taxonomy data.
Results are saved to data/curation/taxonomy_rescue.csv for reference and potential future
data restoration.
"""

import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

def load_csv_by_asmid(filepath: str) -> Dict[str, List[str]]:
    """Load CSV and index by ASMID (first column)."""
    data = {}
    with open(filepath, 'r', newline='', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader)
        for row in reader:
            if row:
                asmid = row[0]
                data[asmid] = row
    return data, header

def find_empty_phylum(samples_data: Dict[str, List[str]], header: List[str]) -> List[Tuple[str, List[str]]]:
    """Find all rows where PHYLUM (column 6) is empty."""
    phylum_idx = header.index('PHYLUM')
    empty_phylum = []
    for asmid, row in samples_data.items():
        if phylum_idx < len(row) and row[phylum_idx].strip() == '':
            empty_phylum.append((asmid, row))
    return empty_phylum

def get_taxonomy_columns(header: List[str]) -> List[int]:
    """Get indices of taxonomy-related columns."""
    tax_cols = ['PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES']
    return [header.index(col) for col in tax_cols if col in header]

def rescue_taxonomy(
    samples_csv: str,
    backup_csv: str,
    output_csv: str
) -> None:
    """
    Identify and rescue missing taxonomy data from backup.

    Args:
        samples_csv: Path to current samples.csv
        backup_csv: Path to backup_samples.csv
        output_csv: Path to save rescue report
    """

    print(f"Loading {samples_csv}...")
    samples_data, samples_header = load_csv_by_asmid(samples_csv)
    print(f"  Loaded {len(samples_data)} entries")

    print(f"Loading {backup_csv}...")
    backup_data, backup_header = load_csv_by_asmid(backup_csv)
    print(f"  Loaded {len(backup_data)} entries")

    # Find rows with empty PHYLUM in current samples
    empty_phylum = find_empty_phylum(samples_data, samples_header)
    print(f"\nFound {len(empty_phylum)} rows with empty PHYLUM in samples.csv")

    # Get taxonomy column indices
    tax_cols = get_taxonomy_columns(samples_header)
    tax_col_names = [samples_header[i] for i in tax_cols]

    # Try to rescue from backup
    rescued = []
    not_found = []
    no_tax_data = []

    for asmid, samples_row in empty_phylum:
        if asmid in backup_data:
            backup_row = backup_data[asmid]
            # Check if backup has any taxonomy data
            has_tax = False
            backup_tax_values = []

            for col_idx in tax_cols:
                if col_idx < len(backup_row):
                    val = backup_row[col_idx].strip()
                    backup_tax_values.append(val)
                    if val:
                        has_tax = True
                else:
                    backup_tax_values.append('')

            if has_tax:
                rescued.append({
                    'asmid': asmid,
                    'species': samples_row[samples_header.index('SPECIES_IN')] if 'SPECIES_IN' in samples_header else '',
                    'ncbi_taxonid': samples_row[samples_header.index('NCBI_TAXONID')] if 'NCBI_TAXONID' in samples_header else '',
                    'backup_data': dict(zip(tax_col_names, backup_tax_values))
                })
            else:
                no_tax_data.append((asmid, samples_row[samples_header.index('SPECIES_IN')] if 'SPECIES_IN' in samples_header else ''))
        else:
            not_found.append((asmid, samples_row[samples_header.index('SPECIES_IN')] if 'SPECIES_IN' in samples_header else ''))

    # Write rescue report
    print(f"\nRescue Results:")
    print(f"  ✓ Taxonomy data found in backup and can be rescued: {len(rescued)}")
    print(f"  ✗ Not found in backup: {len(not_found)}")
    print(f"  ✗ Found in backup but no taxonomy data: {len(no_tax_data)}")

    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, lineterminator='\n')
        # Header
        header_row = ['ASMID', 'SPECIES', 'NCBI_TAXONID'] + tax_col_names
        writer.writerow(header_row)

        # Rescued entries
        for entry in rescued:
            row = [
                entry['asmid'],
                entry['species'],
                entry['ncbi_taxonid']
            ]
            row.extend([entry['backup_data'].get(col, '') for col in tax_col_names])
            writer.writerow(row)

    print(f"\nRescue report saved to: {output_csv}")

    # Print summary of rescued entries
    if rescued:
        print(f"\nRescued entries ({len(rescued)} total):")
        for entry in rescued[:10]:  # Show first 10
            print(f"  {entry['asmid']}: {entry['species']}")
        if len(rescued) > 10:
            print(f"  ... and {len(rescued) - 10} more")

    # Print entries that couldn't be rescued
    if not_found:
        print(f"\nEntries not found in backup ({len(not_found)} total):")
        for asmid, species in not_found[:5]:
            print(f"  {asmid}: {species}")
        if len(not_found) > 5:
            print(f"  ... and {len(not_found) - 5} more")

    if no_tax_data:
        print(f"\nEntries in backup but with no taxonomy data ({len(no_tax_data)} total):")
        for asmid, species in no_tax_data[:5]:
            print(f"  {asmid}: {species}")
        if len(no_tax_data) > 5:
            print(f"  ... and {len(no_tax_data) - 5} more")

if __name__ == '__main__':
    # Determine paths relative to script location
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent

    samples_csv = repo_root / 'samples.csv'
    backup_csv = repo_root / 'backup_samples.csv'
    output_csv = repo_root / 'data' / 'curation' / 'taxonomy_rescue.csv'

    # Ensure output directory exists
    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if not samples_csv.exists():
        print(f"Error: {samples_csv} not found", file=sys.stderr)
        sys.exit(1)
    if not backup_csv.exists():
        print(f"Error: {backup_csv} not found", file=sys.stderr)
        sys.exit(1)

    rescue_taxonomy(str(samples_csv), str(backup_csv), str(output_csv))
