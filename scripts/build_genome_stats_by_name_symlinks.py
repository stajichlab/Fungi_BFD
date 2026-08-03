#!/usr/bin/env python3
"""build_genome_stats_by_name_symlinks.py — human-browsable *_by_name/ symlink trees.

Generates `results/genome_stats_by_name/<GENUS-or-NOGENUS>/<basename>.<type>.<suffix>`
and `results/function_by_name/...` as READ-ONLY, REGENERABLE symlink trees
pointing into the real hash-bucketed canonical store
(`results/genome_stats/*`, `results/function/*`), built entirely from
samples.csv + whatever real bucketed files exist on disk right now.

These trees are NEVER the source of truth -- safe to delete and rebuild at
any time, which is exactly what this script does on every run that isn't a
no-op (full rebuild, not an incremental diff, to avoid ever accumulating
stale symlinks from renamed/removed source data). See
todo/genome_stats_storage_reorg.md §A (T-014), issue #10.

Staleness check
----------------
A stamp file (default `results/.by_name_stamp`) records a hash of:
  - samples.csv's mtime + size
  - every source type directory's own mtime (results/genome_stats/asm_stats,
    results/function/aiupred, etc. -- ~16 stat() calls total)
If unchanged since the last run, this script does nothing (fast no-op).

KNOWN LIMITATION, not a bug: this only detects samples.csv changing or a
BRAND NEW BUCKET subdirectory appearing under a type -- it does NOT detect a
new file being added to an ALREADY-EXISTING bucket (the common case: a new
genome landing in a hash bucket some earlier genome already claimed), because
that only touches the bucket's own mtime, not the parent type directory's.
Checking every individual bucket's mtime would mean up to ~4096 stat() calls
per type (gene_stats/intergenic_stats use 3-hex, 4096 buckets) -- reintroducing
the exact NFS-listing-performance problem this whole reorg exists to fix. Use
--force to bypass the staleness check and always rebuild; a good rule of
thumb is to run this with --force after any batch of new genome_stats/
function output lands, rather than relying on the staleness check alone.

Usage
-----
    # normal run: no-op if nothing changed since the last run
    python3 scripts/build_genome_stats_by_name_symlinks.py

    # always rebuild, ignoring the staleness stamp
    python3 scripts/build_genome_stats_by_name_symlinks.py --force

    # custom paths
    python3 scripts/build_genome_stats_by_name_symlinks.py \
        --samples samples.csv --root results --stamp-file results/.by_name_stamp
"""

import argparse
import csv
import hashlib
import os
import re
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "nextflow" / "bin"))

from genome_stats_paths import hash_bucket_for_type  # noqa: E402

# type -> (target_kind, parent). Mirrors TYPE_CONFIG in
# scripts/one-off/reorg_genome_stats_hash_buckets.py -- duplicated rather than
# imported (that script is a one-off migration tool, this is permanent
# pipeline-supporting tooling; see this repo's CLAUDE.md on scripts/one-off/
# vs scripts/ and the existing make_sample_tag()/cleanStrain() precedent for
# "mirror rather than cross-import" across this codebase).
TYPE_CONFIG = {
    "asm_stats":        ("asmid",    "genome_stats"),
    "BUSCO_genome":     ("asmid",    "genome_stats"),
    "BUSCO_protein":    ("locustag", "genome_stats"),
    "aa_freq":          ("locustag", "genome_stats"),
    "codon_freq":       ("locustag", "genome_stats"),
    "gene_stats":       ("locustag", "genome_stats"),
    "intergenic_stats": ("locustag", "genome_stats"),
    "aiupred":          ("locustag", "function"),
    "merops":           ("locustag", "function"),
    "pfam_hmmscan":     ("locustag", "function"),
    "predgpi":          ("locustag", "function"),
    "signalp":          ("locustag", "function"),
    "targetP":          ("locustag", "function"),
    "tmhmm":            ("locustag", "function"),
    "wolfpsort":        ("locustag", "function"),
}
# cazy already stores one subdirectory per genome (not a hash bucket) and is
# out of scope for the reorg entirely (KNOWN_SUBDIRECTORY_LAYOUT_TYPES in the
# migration script) -- no by-name view generated for it here either.

_PREFIX_DELIMITERS = (".", "_")


def _prefix_lengths(filename: str):
    lengths = [len(filename)]
    lengths.extend(i for i, ch in enumerate(filename) if ch in _PREFIX_DELIMITERS)
    lengths.sort(reverse=True)
    return lengths


def longest_prefix_match(filename: str, candidate_set):
    """Same algorithm as reorg_genome_stats_hash_buckets.py::longest_prefix_match()
    -- O(len(filename)) membership checks, not O(len(candidate_set)). Needed here
    because LOCUSTAG is not fixed-width (8 or 9 hex chars observed in production)."""
    for length in _prefix_lengths(filename):
        candidate = filename[:length]
        if candidate in candidate_set:
            return candidate
    return None


def clean_strain(raw_strain: str) -> str:
    s = (raw_strain or "").strip().replace("'", "").replace('"', "")
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species: str, raw_strain: str) -> str:
    sp = (raw_species or "").strip().replace("'", "").replace('"', "")
    st = clean_strain(raw_strain)
    parts = [p for p in (sp, st) if p]
    tag = "_".join(parts)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def load_samples(samples_path: Path):
    """Return (asmid_set, locustag_set, by_asmid, by_locustag) where
    by_asmid/by_locustag map key -> (genus_or_nogenus, basename)."""
    asmid_set = set()
    locustag_set = set()
    by_asmid = {}
    by_locustag = {}
    with open(samples_path, newline="") as fh:
        for row in csv.DictReader(fh):
            asmid = (row.get("ASMID") or "").strip()
            locustag = (row.get("LOCUSTAG") or "").strip()
            species = (row.get("SPECIES") or "").strip()
            strain = (row.get("STRAIN") or "").strip()
            genus = (row.get("GENUS") or "").strip()
            if not species:
                continue
            genus_dir = re.sub(r"[\s/#\[\]?{}]+", "_", genus) if genus else "NOGENUS"
            basename = make_sample_tag(species, strain)
            if asmid:
                asmid_set.add(asmid)
                by_asmid[asmid] = (genus_dir, basename)
            if locustag:
                locustag_set.add(locustag)
                by_locustag[locustag] = (genus_dir, basename)
    return asmid_set, locustag_set, by_asmid, by_locustag


def compute_signature(samples_path: Path, root: Path) -> str:
    """Cheap staleness signature: samples.csv mtime+size, plus every source
    type directory's own mtime. See module docstring for the known limitation."""
    h = hashlib.sha256()
    st = samples_path.stat()
    h.update(f"samples:{st.st_mtime}:{st.st_size}\n".encode())
    for type_name in sorted(TYPE_CONFIG):
        _kind, parent = TYPE_CONFIG[type_name]
        type_dir = root / parent / type_name
        mtime = type_dir.stat().st_mtime if type_dir.is_dir() else -1
        h.update(f"{type_name}:{mtime}\n".encode())
    return h.hexdigest()


def build_symlink_tree(root: Path, asmid_set, locustag_set, by_asmid, by_locustag):
    """Rebuild every *_by_name/ tree from scratch. Returns (n_linked, n_skipped_collision)."""
    asmid_desc = sorted(asmid_set, key=len, reverse=True)
    locustag_desc = sorted(locustag_set, key=len, reverse=True)

    # Remove and recreate every *_by_name/ root touched by TYPE_CONFIG's parents.
    parents = {parent for _kind, parent in TYPE_CONFIG.values()}
    for parent in parents:
        by_name_root = root / f"{parent}_by_name"
        if by_name_root.exists():
            shutil.rmtree(by_name_root)

    n_linked = 0
    n_skipped = 0
    seen_targets = {}  # symlink_path -> source_path, to detect+skip collisions loudly

    for type_name, (target_kind, parent) in sorted(TYPE_CONFIG.items()):
        src_type_dir = root / parent / type_name
        if not src_type_dir.is_dir():
            continue

        candidate_desc = asmid_desc if target_kind == "asmid" else locustag_desc
        by_key = by_asmid if target_kind == "asmid" else by_locustag

        by_name_root = root / f"{parent}_by_name"

        for bucket_dir in sorted(p for p in src_type_dir.iterdir() if p.is_dir()):
            for entry in sorted(bucket_dir.iterdir()):
                if not entry.is_file():
                    continue
                filename = entry.name
                key = longest_prefix_match(filename, candidate_desc)
                if key is None:
                    print(f"WARN: {entry} does not start with a known "
                          f"{target_kind.upper()} -- skipping (not in samples.csv)",
                          file=sys.stderr)
                    n_skipped += 1
                    continue
                genus_dir, basename = by_key[key]
                suffix = filename[len(key):]
                # type_name is REQUIRED in the symlink filename, not cosmetic:
                # BUSCO_genome and BUSCO_protein both produce
                # "{key}.BUSCO_summary.{lineage}.txt" -- an identical suffix --
                # so without the type embedded, two genuinely different, non-
                # colliding source files (one ASMID-keyed, one LOCUSTAG-keyed)
                # would collapse onto the same by-name target. Confirmed live:
                # this was the majority of collisions on the first production
                # run before this fix.
                link_path = by_name_root / genus_dir / f"{basename}.{type_name}{suffix}"

                if link_path in seen_targets:
                    print(f"WARN: symlink target collision, skipping: {link_path} "
                          f"already points to {seen_targets[link_path]}, "
                          f"also wanted by {entry}", file=sys.stderr)
                    n_skipped += 1
                    continue

                link_path.parent.mkdir(parents=True, exist_ok=True)
                relative_target = os.path.relpath(entry, link_path.parent)
                link_path.symlink_to(relative_target)
                seen_targets[link_path] = entry
                n_linked += 1

    return n_linked, n_skipped


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", default="samples.csv", help="samples.csv path [samples.csv]")
    parser.add_argument("--root", default="results", help="root containing genome_stats/ and function/ [results]")
    parser.add_argument("--stamp-file", default=None, help="staleness stamp path [<root>/.by_name_stamp]")
    parser.add_argument("--force", action="store_true", help="rebuild even if the staleness signature is unchanged")
    args = parser.parse_args()

    root = Path(args.root)
    samples_path = Path(args.samples)
    if not samples_path.exists():
        print(f"ERROR: samples file not found: {samples_path}", file=sys.stderr)
        sys.exit(1)

    stamp_path = Path(args.stamp_file) if args.stamp_file else root / ".by_name_stamp"

    signature = compute_signature(samples_path, root)
    if not args.force and stamp_path.exists() and stamp_path.read_text().strip() == signature:
        print(f"Up to date (signature unchanged since last run) -- nothing to do. "
              f"Use --force to rebuild anyway.")
        return

    asmid_set, locustag_set, by_asmid, by_locustag = load_samples(samples_path)
    print(f"Loaded samples.csv: {len(asmid_set)} ASMIDs, {len(locustag_set)} LOCUSTAGs")

    n_linked, n_skipped = build_symlink_tree(root, asmid_set, locustag_set, by_asmid, by_locustag)
    print(f"Created {n_linked} symlinks, {n_skipped} skipped (unresolvable key or target collision)")

    stamp_path.parent.mkdir(parents=True, exist_ok=True)
    stamp_path.write_text(signature)
    print(f"Stamp updated: {stamp_path}")


if __name__ == "__main__":
    main()
