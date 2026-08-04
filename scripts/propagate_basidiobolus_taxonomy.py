#!/usr/bin/env python3
"""
Propagate known Basidiobolus taxonomy to rescue entries.

Basidiobolus is a well-characterized fungal genus. All species share the same
higher-level taxonomy. This script fills in the taxonomy for Basidiobolus entries
based on the established taxonomy.
"""

import csv
from pathlib import Path

# Established taxonomy for Basidiobolus genus
BASIDIOBOLUS_TAXONOMY = {
    'PHYLUM': 'Basidiomycota',
    'SUBPHYLUM': 'Agaricomycotina',
    'CLASS': 'Exobasidiomycetes',
    'SUBCLASS': '',
    'ORDER': 'Basidiobolales',
    'FAMILY': 'Basidiobolaceae',
    'GENUS': 'Basidiobolus'
}

def get_basidiobolus_entries(samples_csv: str):
    """Get all Basidiobolus entries with empty phylum from samples.csv."""
    entries = []
    with open(samples_csv, 'r', newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if 'Basidiobolus' in row.get('SPECIES_IN', '') and not row.get('PHYLUM', '').strip():
                entries.append({
                    'asmid': row['ASMID'],
                    'species': row['SPECIES_IN'],
                    'ncbi_taxonid': row['NCBI_TAXONID'],
                    'strain': row.get('STRAIN', '')
                })
    return entries

def save_propagated_taxonomy(entries, output_csv: str) -> None:
    """Save propagated taxonomy to CSV."""
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'ASMID', 'SPECIES', 'NCBI_TAXONID', 'STRAIN',
            'PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS', 'SPECIES_RANK'
        ], lineterminator='\n')
        writer.writeheader()

        for entry in entries:
            # Extract species name from SPECIES_IN (last component)
            species_parts = entry['species'].rsplit(' ', 1)
            species_rank = species_parts[-1] if len(species_parts) > 1 else entry['species']

            writer.writerow({
                'ASMID': entry['asmid'],
                'SPECIES': entry['species'],
                'NCBI_TAXONID': entry['ncbi_taxonid'],
                'STRAIN': entry['strain'],
                'PHYLUM': BASIDIOBOLUS_TAXONOMY['PHYLUM'],
                'SUBPHYLUM': BASIDIOBOLUS_TAXONOMY['SUBPHYLUM'],
                'CLASS': BASIDIOBOLUS_TAXONOMY['CLASS'],
                'SUBCLASS': BASIDIOBOLUS_TAXONOMY['SUBCLASS'],
                'ORDER': BASIDIOBOLUS_TAXONOMY['ORDER'],
                'FAMILY': BASIDIOBOLUS_TAXONOMY['FAMILY'],
                'GENUS': BASIDIOBOLUS_TAXONOMY['GENUS'],
                'SPECIES_RANK': species_rank
            })

if __name__ == '__main__':
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent

    samples_csv = repo_root / 'samples.csv'
    output_csv = repo_root / 'data' / 'curation' / 'basidiobolus_taxonomy_rescue.csv'

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    print(f"Fetching Basidiobolus entries from {samples_csv}...")
    entries = get_basidiobolus_entries(str(samples_csv))
    print(f"Found {len(entries)} Basidiobolus entries\n")

    print(f"Propagating known Basidiobolus taxonomy:")
    print(f"  Phylum: {BASIDIOBOLUS_TAXONOMY['PHYLUM']}")
    print(f"  Subphylum: {BASIDIOBOLUS_TAXONOMY['SUBPHYLUM']}")
    print(f"  Class: {BASIDIOBOLUS_TAXONOMY['CLASS']}")
    print(f"  Order: {BASIDIOBOLUS_TAXONOMY['ORDER']}")
    print(f"  Family: {BASIDIOBOLUS_TAXONOMY['FAMILY']}")
    print(f"  Genus: {BASIDIOBOLUS_TAXONOMY['GENUS']}\n")

    print("Entries to be propagated:")
    for i, entry in enumerate(entries, 1):
        print(f"  {i:2}. {entry['asmid']}: {entry['species']}")

    save_propagated_taxonomy(entries, str(output_csv))
    print(f"\n✓ Saved {len(entries)} propagated entries to {output_csv}")

