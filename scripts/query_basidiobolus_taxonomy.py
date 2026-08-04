#!/usr/bin/env python3
"""
Query NCBI taxonomy for Basidiobolus entries and add to rescue file.

This script queries the NCBI Taxonomy database for Basidiobolus species
using their NCBI taxon IDs and retrieves the full taxonomic lineage.
"""

import csv
import json
import urllib.request
import urllib.error
import time
import sys
from pathlib import Path
from typing import Dict, List, Optional

def query_ncbi_taxonomy(taxon_id: int, retries: int = 3) -> Optional[Dict]:
    """
    Query NCBI Taxonomy database for a given taxon ID.

    Args:
        taxon_id: NCBI taxonomy ID
        retries: Number of retry attempts

    Returns:
        Dictionary with taxonomy lineage or None if failed
    """
    url = f"https://www.ncbi.nlm.nih.gov/taxonomy/?format=json&id={taxon_id}"

    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=5) as response:
                data = json.loads(response.read().decode('utf-8'))
                if 'result' in data and str(taxon_id) in data['result']:
                    return data['result'][str(taxon_id)]
            time.sleep(0.5)  # Rate limiting
            return None
        except urllib.error.HTTPError as e:
            if attempt < retries - 1:
                time.sleep(1)
            else:
                print(f"Warning: Failed to query taxon {taxon_id}: {e}", file=sys.stderr)
                return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(1)
            else:
                print(f"Warning: Error querying taxon {taxon_id}: {e}", file=sys.stderr)
                return None

    return None

def extract_lineage(taxonomy_data: Dict) -> Dict[str, str]:
    """
    Extract taxonomic ranks from NCBI taxonomy data.

    Args:
        taxonomy_data: Taxonomy data from NCBI

    Returns:
        Dictionary with rank names as keys
    """
    ranks = {}
    if 'lineage' in taxonomy_data:
        for item in taxonomy_data['lineage']:
            if 'rank' in item and 'scientificname' in item:
                rank = item['rank'].lower()
                if rank in ['phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'genus']:
                    ranks[rank] = item['scientificname']

    return ranks

def get_basidiobolus_entries(samples_csv: str) -> List[Dict]:
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

def query_and_save_basidiobolus(samples_csv: str, output_csv: str) -> None:
    """
    Query NCBI for Basidiobolus taxonomy and save results.

    Args:
        samples_csv: Path to samples.csv
        output_csv: Path to save results
    """
    print("Fetching Basidiobolus entries from samples.csv...")
    entries = get_basidiobolus_entries(samples_csv)
    print(f"Found {len(entries)} Basidiobolus entries with empty taxonomy")

    rescued = []

    for i, entry in enumerate(entries, 1):
        print(f"[{i}/{len(entries)}] Querying taxon {entry['ncbi_taxonid']} for {entry['species']}...", end=' ', flush=True)

        tax_data = query_ncbi_taxonomy(int(entry['ncbi_taxonid']))

        if tax_data:
            lineage = extract_lineage(tax_data)
            rescued.append({
                'asmid': entry['asmid'],
                'species': entry['species'],
                'ncbi_taxonid': entry['ncbi_taxonid'],
                'phylum': lineage.get('phylum', ''),
                'subphylum': lineage.get('subphylum', ''),
                'class': lineage.get('class', ''),
                'subclass': lineage.get('subclass', ''),
                'order': lineage.get('order', ''),
                'family': lineage.get('family', ''),
                'genus': lineage.get('genus', '')
            })
            print(f"✓ {lineage.get('phylum', 'Unknown')}")
        else:
            print("✗ Query failed")

    # Write results
    with open(output_csv, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=[
            'ASMID', 'SPECIES', 'NCBI_TAXONID',
            'PHYLUM', 'SUBPHYLUM', 'CLASS', 'SUBCLASS', 'ORDER', 'FAMILY', 'GENUS'
        ])
        writer.writeheader()
        for entry in rescued:
            writer.writerow({
                'ASMID': entry['asmid'],
                'SPECIES': entry['species'],
                'NCBI_TAXONID': entry['ncbi_taxonid'],
                'PHYLUM': entry['phylum'],
                'SUBPHYLUM': entry['subphylum'],
                'CLASS': entry['class'],
                'SUBCLASS': entry['subclass'],
                'ORDER': entry['order'],
                'FAMILY': entry['family'],
                'GENUS': entry['genus']
            })

    print(f"\n✓ Saved {len(rescued)}/{len(entries)} successfully queried entries to {output_csv}")

    # Show summary
    if rescued:
        # Check if all have same phylum (for propagation)
        phyla = set(e['phylum'] for e in rescued if e['phylum'])
        if len(phyla) == 1:
            print(f"✓ All entries belong to: {phyla.pop()}")

if __name__ == '__main__':
    script_dir = Path(__file__).parent
    repo_root = script_dir.parent

    samples_csv = repo_root / 'samples.csv'
    output_csv = repo_root / 'data' / 'curation' / 'basidiobolus_taxonomy_rescue.csv'

    output_csv.parent.mkdir(parents=True, exist_ok=True)

    if not samples_csv.exists():
        print(f"Error: {samples_csv} not found", file=sys.stderr)
        sys.exit(1)

    query_and_save_basidiobolus(str(samples_csv), str(output_csv))
