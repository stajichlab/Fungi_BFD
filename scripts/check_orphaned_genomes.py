#!/usr/bin/env python3
"""
Detect genomes on disk that are not in samples.csv.
Run after each samples.csv regeneration to find orphaned files.

Usage:
    python3 scripts/check_orphaned_genomes.py [--remove-log] [--archive-prefix PREFIX]
"""

import csv
import os
import re
import sys
from pathlib import Path
from datetime import datetime
from collections import defaultdict

def load_samples_asmids(samples_file='samples.csv'):
    """Extract all ASMIDs from samples.csv."""
    asmids = set()
    try:
        with open(samples_file) as f:
            reader = csv.DictReader(f)
            for row in reader:
                asmids.add(row['ASMID'])
    except FileNotFoundError:
        print(f"ERROR: {samples_file} not found", file=sys.stderr)
        sys.exit(1)
    return asmids

def extract_asmid(path):
    """Extract ASMID (GCF/GCA_*) from a directory/file path."""
    match = re.search(r'(GC[AF]_\d+\.\d+)', path)
    return match.group(1) if match else None

def find_orphaned_genomes(search_dirs=None, max_depth=3):
    """Search for genome directories/files not in samples.csv."""
    if search_dirs is None:
        search_dirs = [
            'data/genomes',
            'data/raw',
            'data/processed',
            'analysis',
            'misc'
        ]

    samples_asmids = load_samples_asmids()
    print(f"[*] Loaded {len(samples_asmids)} ASMIDs from samples.csv")

    orphaned = defaultdict(list)
    scanned = 0

    for search_dir in search_dirs:
        if not os.path.isdir(search_dir):
            continue

        print(f"[*] Scanning {search_dir}...")

        for root, dirs, files in os.walk(search_dir):
            # Limit depth
            depth = root.count(os.sep) - search_dir.count(os.sep)
            if depth > max_depth:
                dirs.clear()
                continue

            # Check directories
            for d in dirs:
                full_path = os.path.join(root, d)
                asmid = extract_asmid(full_path)

                if asmid and asmid not in samples_asmids:
                    size_mb = get_size(full_path) / (1024 ** 2)
                    orphaned[asmid].append({
                        'path': full_path,
                        'size_mb': size_mb,
                        'type': 'directory'
                    })
                    scanned += 1

            # Check files
            for f in files:
                full_path = os.path.join(root, f)
                asmid = extract_asmid(full_path)

                if asmid and asmid not in samples_asmids:
                    size_mb = os.path.getsize(full_path) / (1024 ** 2)
                    orphaned[asmid].append({
                        'path': full_path,
                        'size_mb': size_mb,
                        'type': 'file'
                    })
                    scanned += 1

    return orphaned, samples_asmids

def get_size(path):
    """Calculate total size of directory."""
    total = 0
    try:
        for dirpath, dirnames, filenames in os.walk(path):
            for f in filenames:
                try:
                    total += os.path.getsize(os.path.join(dirpath, f))
                except OSError:
                    pass
    except OSError:
        pass
    return total

def load_removed_asmids(removed_file='data/curation/removed_asmids.csv'):
    """Load previously recorded removed ASMIDs."""
    removed = {}
    if os.path.exists(removed_file):
        try:
            with open(removed_file) as f:
                reader = csv.DictReader(f)
                for row in reader:
                    removed[row['ASMID']] = row
        except Exception as e:
            print(f"WARNING: Could not read {removed_file}: {e}", file=sys.stderr)
    return removed

def main():
    orphaned, current_asmids = find_orphaned_genomes()
    removed_asmids = load_removed_asmids()

    if not orphaned:
        print("[✓] No orphaned genomes found!")
        return 0

    print(f"\n[!] Found {len(orphaned)} orphaned ASMID groups:")
    print()

    # Analyze orphaned
    total_space = 0
    expected_removals = 0
    unexpected_orphans = []

    for asmid in sorted(orphaned.keys()):
        items = orphaned[asmid]
        total_size = sum(item['size_mb'] for item in items)
        total_space += total_size

        # Check if this was expected
        was_removed = asmid in removed_asmids
        status = "✓ EXPECTED" if was_removed else "✗ UNEXPECTED"

        print(f"{status}: {asmid}")
        print(f"  Size: {total_size:.1f} MB ({len(items)} items)")

        for item in items[:3]:  # Show first 3 files
            print(f"    - {item['path']} ({item['size_mb']:.1f} MB)")

        if len(items) > 3:
            print(f"    ... and {len(items) - 3} more items")

        if was_removed:
            removal_info = removed_asmids[asmid]
            print(f"  Reason: {removal_info.get('REASON', 'N/A')}")
            print(f"  Date removed: {removal_info.get('DATE_REMOVED', 'N/A')}")
            expected_removals += 1
        else:
            unexpected_orphans.append((asmid, total_size, items))

        print()

    # Summary
    print("=" * 70)
    print(f"SUMMARY:")
    print(f"  Total orphaned ASMIDs: {len(orphaned)}")
    print(f"    - Expected removals: {expected_removals}")
    print(f"    - Unexpected orphans: {len(unexpected_orphans)}")
    print(f"  Total space used: {total_space:.1f} MB ({total_space/1024:.1f} GB)")
    print()

    if unexpected_orphans:
        print("UNEXPECTED ORPHANS (may need investigation):")
        for asmid, size, items in unexpected_orphans:
            print(f"  {asmid}: {size:.1f} MB")
            print(f"    Paths: {', '.join(item['path'] for item in items[:2])}")
        print()
        print("ACTION: Check if these should be in removed_asmids.csv")
        return 1

    return 0

if __name__ == '__main__':
    sys.exit(main())
