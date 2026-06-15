#!/usr/bin/env python3
"""
Rename genome_stats output files from LOCUSTAG-prefixed names to
SPECIES_STRAIN basenames, using the same sanitisation rules as BFD.nf:

  strain   = first ';'-separated token, single-quotes stripped, ':' → ' '
  basename = f"{species}_{strain}" with whitespace / '/' / '#' → '_'

Subdirectories handled:
  aa_freq/          {LOCUSTAG}.aa_freq.csv.gz
  codon_freq/       {LOCUSTAG}.codon_freq.csv.gz
  gene_stats/       {LOCUSTAG}.gene_*.csv.gz
  intergenic_stats/ {LOCUSTAG}.gene_intergenic_distances.csv.gz

asm_reports/ files are named after ASMID and are left untouched.

Usage:
  python misc_scripts/renamefiles_locustag_to_speciesstrain.py \\
      [--samples samples.csv] \\
      [--genome-stats results/genome_stats] \\
      [--dry-run]
"""

import argparse
import csv
import re
import sys
from pathlib import Path


def make_basename(species: str, strain: str) -> str:
    """Replicate BFD.nf basename construction."""
    strain = strain.split(";")[0].strip().replace("'", "").replace(":", " ")
    parts = [p for p in (species.strip(), strain) if p]
    raw = "_".join(parts)
    return re.sub(r"[\s/#]+", "_", raw)


def build_locustag_map(samples_csv: Path) -> dict[str, str]:
    """Return {locustag: basename} for every row in samples.csv."""
    mapping: dict[str, str] = {}
    with samples_csv.open(newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            locustag = (row.get("LOCUSTAG") or "").strip().strip("\r\n")
            species  = (row.get("SPECIES") or "").strip()
            strain   = (row.get("STRAIN") or "").strip()
            if not locustag:
                continue
            mapping[locustag] = make_basename(species, strain)
    return mapping


LOCUSTAG_RE = re.compile(r"^([0-9A-Fa-f]{8})(\..*)")


def rename_dir(directory: Path, mapping: dict[str, str], dry_run: bool) -> tuple[int, int]:
    """Rename files in *directory* whose names start with a known LOCUSTAG.

    Returns (renamed, skipped) counts.
    """
    renamed = skipped = 0
    for fpath in sorted(directory.iterdir()):
        if not fpath.is_file():
            continue
        m = LOCUSTAG_RE.match(fpath.name)
        if not m:
            skipped += 1
            continue
        locustag, suffix = m.group(1).upper(), m.group(2)
        basename = mapping.get(locustag)
        if basename is None:
            print(f"  WARN  no mapping for {locustag} ({fpath.name})", file=sys.stderr)
            skipped += 1
            continue
        new_path = fpath.parent / f"{basename}{suffix}"
        if new_path == fpath:
            skipped += 1
            continue
        if new_path.exists():
            print(f"  WARN  target exists, skipping: {new_path.name}", file=sys.stderr)
            skipped += 1
            continue
        print(f"  {'DRY ' if dry_run else ''}RENAME  {fpath.name}  →  {new_path.name}")
        if not dry_run:
            fpath.rename(new_path)
        renamed += 1
    return renamed, skipped


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", default="samples.csv",
                        help="Path to samples.csv (default: %(default)s)")
    parser.add_argument("--genome-stats", default="results/genome_stats",
                        help="Path to genome_stats directory (default: %(default)s)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print what would be renamed without making changes")
    args = parser.parse_args()

    samples_csv   = Path(args.samples)
    genome_stats  = Path(args.genome_stats)
    subdirs       = ["aa_freq", "codon_freq", "gene_stats", "intergenic_stats"]

    if not samples_csv.exists():
        sys.exit(f"ERROR: samples file not found: {samples_csv}")
    if not genome_stats.is_dir():
        sys.exit(f"ERROR: genome_stats directory not found: {genome_stats}")

    print(f"Loading sample map from {samples_csv} …")
    mapping = build_locustag_map(samples_csv)
    print(f"  {len(mapping)} locustag → basename entries loaded")

    total_renamed = total_skipped = 0
    for sub in subdirs:
        d = genome_stats / sub
        if not d.is_dir():
            print(f"  SKIP  {d} (not found)")
            continue
        print(f"\nProcessing {d} …")
        r, s = rename_dir(d, mapping, dry_run=args.dry_run)
        print(f"  renamed={r}  skipped/warned={s}")
        total_renamed += r
        total_skipped += s

    print(f"\nDone.  Total renamed={total_renamed}  skipped/warned={total_skipped}")
    if args.dry_run:
        print("(dry-run: no files were actually renamed)")


if __name__ == "__main__":
    main()
