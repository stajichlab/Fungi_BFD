#!/usr/bin/env python3
"""
Sync completed funannotate predict results from Nextflow work directories to
genome_annotation/.

Scans work/funannotate/<2char>/<fullhash>/<SpeciesName>/ for a predict_results/*.gbk file.
If found, rsyncs the entire species folder into genome_annotation/<SpeciesName>/.
Skips work dirs that have no gbk (incomplete or failed runs).
Skips species that already have a gbk in genome_annotation (unless --force).

Usage:
    python scripts/sync_predict_results.py [--dry-run] [--force] [--work-dir PATH] [--dest-dir PATH]
"""

import argparse
import subprocess
import sys
from pathlib import Path
from typing import Optional


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dry-run", action="store_true",
                   help="Print what would be done without copying anything.")
    p.add_argument("--force", action="store_true",
                   help="Overwrite destination even if a gbk already exists there.")
    p.add_argument("--work-dir", type=Path, default=Path("work/funannotate"),
                   help="Nextflow work directory containing 2-level hash subdirectories (default: work/funannotate).")
    p.add_argument("--dest-dir", type=Path, default=Path("genome_annotation"),
                   help="Destination base directory (default: genome_annotation).")
    return p.parse_args()


def iter_work_species_dirs(work_base: Path):
    """Yield all <SpeciesName> dirs found 2 levels deep: work/<2char>/<fullhash>/<SpeciesName>/"""
    for prefix_dir in work_base.iterdir():
        if not prefix_dir.is_dir():
            continue
        for hash_dir in prefix_dir.iterdir():
            if not hash_dir.is_dir():
                continue
            for species_dir in hash_dir.iterdir():
                if species_dir.is_dir() and not species_dir.name.startswith("."):
                    yield species_dir


def find_best_work_dir(work_base: Path, species_name: str) -> Optional[Path]:
    """
    Among all work/<2char>/<hash>/<species_name>/ candidates that contain a gbk,
    return the one with the most recent mtime (most recent successful run).
    """
    candidates = []
    for prefix_dir in work_base.iterdir():
        if not prefix_dir.is_dir():
            continue
        for hash_dir in prefix_dir.iterdir():
            if not hash_dir.is_dir():
                continue
            species_dir = hash_dir / species_name
            if not species_dir.is_dir():
                continue
            predict_results = species_dir / "predict_results"
            gbks = list(predict_results.glob("*.gbk")) if predict_results.is_dir() else []
            if gbks:
                candidates.append((max(g.stat().st_mtime for g in gbks), species_dir))
    if not candidates:
        return None
    candidates.sort(key=lambda x: x[0], reverse=True)
    return candidates[0][1]


def main():
    args = parse_args()
    work_base: Path = args.work_dir
    dest_base: Path = args.dest_dir

    if not work_base.is_dir():
        sys.exit(f"ERROR: work directory not found: {work_base}")
    if not dest_base.is_dir() and not args.dry_run:
        dest_base.mkdir(parents=True)

    # Collect all unique species names across all 2-level hash directories
    species_names = set()
    for species_dir in iter_work_species_dirs(work_base):
        species_names.add(species_dir.name)

    copied = skipped_no_gbk = skipped_exists = 0

    for species in sorted(species_names):
        src_dir = find_best_work_dir(work_base, species)
        if src_dir is None:
            skipped_no_gbk += 1
            continue

        dest_dir = dest_base / species

        # Skip if destination already has a gbk and --force not set
        if dest_dir.is_dir() and not args.force:
            existing_gbks = list((dest_dir / "predict_results").glob("*.gbk")) if (dest_dir / "predict_results").is_dir() else []
            if existing_gbks:
                print(f"SKIP (exists)  {species}")
                skipped_exists += 1
                continue

        action = "DRY-RUN" if args.dry_run else "SYNC"
        print(f"{action}  {src_dir}  ->  {dest_dir}")

        if not args.dry_run:
            dest_dir.mkdir(parents=True, exist_ok=True)
            result = subprocess.run(
                ["rsync", "-a", "--exclude=.command*", "--exclude=.exitcode",
                 "--exclude=.fusion*", "--exclude=.nxf*",
                 f"{src_dir}/", f"{dest_dir}/"],
                capture_output=True, text=True,
            )
            if result.returncode != 0:
                print(f"  ERROR rsync failed:\n{result.stderr}", file=sys.stderr)
            else:
                copied += 1

    print(f"\nDone. Synced: {copied}  Skipped (no gbk): {skipped_no_gbk}  Skipped (already exists): {skipped_exists}")
    if args.dry_run:
        print("(dry-run: no files were copied)")


if __name__ == "__main__":
    main()
