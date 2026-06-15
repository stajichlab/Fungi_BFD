#!/usr/bin/env python3
"""
Update SRA query cache files with platform information.

For each file in rnaseq_reads/sra_query/ that:
  - has at least one data record (more than just a header line)
  - lacks a 'platform' column

fetches platform info via efetch (ncbi_edirect) and rewrites the file
with a new 'platform' column appended.

Outputs:
  UPDATED_SRA_CACHE.txt  -- filenames that were updated
  UPDATED_WITH_BGI.txt   -- species names where any record used a BGI platform

Usage:
  module load ncbi_edirect
  python nextflow/bin/update_platform_sra_query_cache.py [--sra-dir PATH] [--outdir PATH]
"""

import argparse
import csv
import io
import os
import subprocess
import sys
import time
from pathlib import Path


def find_efetch():
    result = subprocess.run(["which", "efetch"], capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(
            "ERROR: efetch not found in PATH.\n"
            "Run: module load ncbi_edirect\n"
            "then re-run this script."
        )
    return result.stdout.strip()


def fetch_platform_map(accessions, retries=3, sleep_base=2):
    """
    Query efetch runinfo for a list of SRA run accessions.
    Returns dict {accession: platform_string} for all accessions returned.
    """
    id_str = ",".join(accessions)
    for attempt in range(retries):
        try:
            proc = subprocess.run(
                ["efetch", "-db", "sra", "-id", id_str, "-format", "runinfo"],
                capture_output=True, text=True, timeout=120,
            )
            if proc.returncode == 0 and proc.stdout.strip():
                reader = csv.DictReader(io.StringIO(proc.stdout.strip()))
                platform_map = {}
                for row in reader:
                    run = row.get("Run", "").strip()
                    platform = row.get("Platform", "").strip().upper()
                    if run:
                        platform_map[run] = platform
                if platform_map:
                    return platform_map
            # Empty or failed response — wait and retry
            if attempt < retries - 1:
                time.sleep(sleep_base ** attempt)
        except subprocess.TimeoutExpired:
            print(f"  WARNING: efetch timed out for {id_str} (attempt {attempt+1})", file=sys.stderr)
            if attempt < retries - 1:
                time.sleep(sleep_base ** attempt)
        except Exception as exc:
            print(f"  WARNING: efetch error for {id_str}: {exc}", file=sys.stderr)
            if attempt < retries - 1:
                time.sleep(sleep_base ** attempt)
    return {}


def process_file(csv_path, dry_run=False):
    """
    Read one cache CSV file.
    Returns (updated: bool, has_bgi: bool, new_content: str or None).
    updated=False when the file already has a platform column or has no data rows.
    """
    with open(csv_path, newline="") as fh:
        content = fh.read()

    lines = content.splitlines()
    if len(lines) < 2:
        # Header only or empty — nothing to do
        return False, False, None

    reader = csv.DictReader(io.StringIO(content))
    fieldnames = reader.fieldnames or []

    if "platform" in [f.lower() for f in fieldnames]:
        return False, False, None

    rows = list(reader)
    if not rows:
        return False, False, None

    accessions = [r["sra_accession"].strip() for r in rows if r.get("sra_accession", "").strip()]
    if not accessions:
        return False, False, None

    print(f"  Fetching platform for {len(accessions)} accession(s): {', '.join(accessions)}")
    platform_map = fetch_platform_map(accessions)

    if not platform_map:
        print(f"  WARNING: no platform data returned for {csv_path.name}", file=sys.stderr)
        return False, False, None

    # Build updated rows
    new_fieldnames = fieldnames + ["platform"]
    out = io.StringIO()
    writer = csv.DictWriter(out, fieldnames=new_fieldnames, lineterminator="\n")
    writer.writeheader()
    has_bgi = False
    for row in rows:
        acc = row.get("sra_accession", "").strip()
        platform = platform_map.get(acc, "")
        if "BGI" in platform or "BGISEQ" in platform:
            has_bgi = True
        row["platform"] = platform
        writer.writerow(row)

    return True, has_bgi, out.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--sra-dir",
        default="rnaseq_reads/sra_query",
        help="Directory containing *.sra_query.csv files (default: rnaseq_reads/sra_query)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Directory for UPDATED_SRA_CACHE.txt and UPDATED_WITH_BGI.txt (default: .)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Report what would be updated without writing any files",
    )
    args = parser.parse_args()

    find_efetch()

    sra_dir = Path(args.sra_dir)
    if not sra_dir.is_dir():
        sys.exit(f"ERROR: SRA query directory not found: {sra_dir}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cache_files = sorted(sra_dir.glob("*.sra_query.csv"))
    print(f"Found {len(cache_files)} cache files in {sra_dir}")

    updated_files = []
    bgi_species = []

    for csv_path in cache_files:
        print(f"Checking {csv_path.name} ...")
        try:
            updated, has_bgi, new_content = process_file(csv_path, dry_run=args.dry_run)
        except Exception as exc:
            print(f"  ERROR processing {csv_path.name}: {exc}", file=sys.stderr)
            continue

        if updated:
            updated_files.append(csv_path.name)
            if has_bgi:
                # Derive species name from filename (strip .sra_query.csv)
                species = csv_path.name.replace(".sra_query.csv", "")
                bgi_species.append(species)
                print(f"  -> Updated (BGI platform detected)")
            else:
                print(f"  -> Updated")

            if not args.dry_run and new_content:
                csv_path.write_text(new_content)
        else:
            print(f"  -> Skipped (already has platform or no data rows)")

        # Be polite to NCBI — small pause between requests
        time.sleep(0.4)

    updated_path = outdir / "UPDATED_SRA_CACHE.txt"
    bgi_path = outdir / "UPDATED_WITH_BGI.txt"

    if args.dry_run:
        print("\n[dry-run] Would write:")
        print(f"  {updated_path}: {len(updated_files)} entries")
        print(f"  {bgi_path}: {len(bgi_species)} entries")
    else:
        updated_path.write_text("\n".join(updated_files) + ("\n" if updated_files else ""))
        bgi_path.write_text("\n".join(bgi_species) + ("\n" if bgi_species else ""))
        print(f"\nDone. {len(updated_files)} files updated, {len(bgi_species)} with BGI platform.")
        print(f"  Summary: {updated_path}")
        print(f"  BGI list: {bgi_path}")


if __name__ == "__main__":
    main()
