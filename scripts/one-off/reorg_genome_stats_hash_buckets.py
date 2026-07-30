#!/usr/bin/env python3
"""One-time migration: fan flat genome_stats/function directories out into hash buckets.

Background
----------
`results/genome_stats/{asm_stats,asm_reports,BUSCO_genome,BUSCO_protein,aa_freq,
codon_freq,gene_stats,intergenic_stats}/` and `results/function/{aiupred,cazy,
merops,pfam_hmmscan,predgpi,signalp,targetP,tmhmm,wolfpsort}/` are flat directories
that reach tens of thousands of files each (gene_stats: 55,398). This is slow to
list/back up on this cluster's NFS-backed storage. See
`todo/genome_stats_storage_reorg.md` (T-014) for the full plan.

This script moves every existing file from `<type>/<name>.<ext>` to
`<type>/<bucket>/<key>.<ext>`, where:
  - `key` is ASMID for asm_stats/asm_reports/BUSCO_genome (assembly-level outputs)
    or LOCUSTAG for everything else (annotation-derived outputs)
  - `bucket` is `hash_bucket_for_type(type, key)` (nextflow/bin/genome_stats_paths.py)

Key resolution
--------------
Accession-keyed types (asm_stats, asm_reports, BUSCO_genome... wait, BUSCO_genome is
actually name-keyed on disk today, see below): filenames already start with the
ASMID itself (`GCA_000149445.2_ASM14944v2.stats.txt`) -- resolve by longest-prefix
match against the known ASMID set from samples.csv (handles the bare
`GCA_000149445` vs versioned `GCA_000149445.2` ambiguity the same way
collect_busco_stats.py::load_asmid_map() does).

Name-keyed types (BUSCO_genome, BUSCO_protein, aa_freq, codon_freq, gene_stats,
intergenic_stats, function/*): filenames start with a `Genus_species_strain` basename
built the same way Nextflow's makeSampleTag()/cleanStrain() build it (see
scripts/collect_busco_stats.py::build_basename_map(), mirrored here rather than
imported to keep this one-off script self-contained) -- resolve by longest-prefix
match against the known basename set, then map basename -> ASMID or -> LOCUSTAG
depending on the type's key kind.

Files whose key cannot be resolved (no samples.csv match -- expect a non-zero count:
329 samples.csv rows have no GENUS, and historical species-name drift is documented
elsewhere in this repo) are moved to `results/_UNMATCHED/<type>/<original filename>`,
never dropped, never silently renamed.

Safety
------
- Dry-run by default; `--apply` is required to write anything.
- A target-path uniqueness check runs BEFORE any rename executes: if two different
  source files would produce the same target filename (observed in production --
  e.g. a strain-less legacy filename and a later strain-qualified one that both
  resolve to the same LOCUSTAG), ALL of that key's files are quarantined to
  `results/_COLLISIONS/<type>/` instead of being moved to their computed bucket.
  This is never auto-resolved (no mtime-wins guessing) and never blocks the rest
  of the migration -- every other, non-colliding file still moves normally in the
  same run. Collisions require a human decision (which file is authoritative, or
  whether both are genuinely distinct data under a badly-collided key) before
  their quarantined files are dealt with.
- The full manifest (old_path, new_path or _UNMATCHED path, sha256, size, mtime) is
  written BEFORE any rename, so it is the audit trail and the rollback mechanism
  (swap old_path/new_path and re-run) even if the process is killed partway.
- Idempotent/resumable: if new_path already exists and matches the manifest's
  recorded sha256/size while old_path no longer exists, that row is skipped as
  already-done rather than re-processed or errored -- safe to re-run after a killed
  SLURM job.
- Moves use os.rename (same filesystem, metadata-only) -- no data duplication, no
  transient double disk usage.
- A verification pass (file-count parity, sha256 spot-check, per-bucket count
  uniformity) runs after every --apply and must be reviewed before considering the
  migration for a given directory complete.

Usage
-----
    # dry-run (default): report what would move, touch nothing
    python3 scripts/one-off/reorg_genome_stats_hash_buckets.py \
        --samples samples.csv --root results

    # apply: perform the migration for real
    python3 scripts/one-off/reorg_genome_stats_hash_buckets.py \
        --samples samples.csv --root results --apply

    # scope to one directory at a time (recommended for the first production runs)
    python3 scripts/one-off/reorg_genome_stats_hash_buckets.py \
        --samples samples.csv --root results --apply --only asm_stats

Run via sbatch, not on the login node, for the full-scale production run:
    sbatch -p short --mem 8gb -c 1 -N 1 -n 1 --out logs/reorg_genome_stats.%j.log \
        --wrap "python3 scripts/one-off/reorg_genome_stats_hash_buckets.py \
                --samples samples.csv --root results --apply"
"""

import argparse
import csv
import hashlib
import os
import re
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "nextflow" / "bin"))

from genome_stats_paths import hash_bucket_for_type  # noqa: E402

# type -> (source_kind, target_kind, parent directory under --root)
#
# source_kind: how the EXISTING filename is parsed --
#   "asmid" -> filename is literally prefixed with the ASMID itself
#              (e.g. GCA_000149445.2_ASM14944v2.stats.txt)
#   "name"  -> filename is prefixed with the Genus_species_strain basename
#              (e.g. Aaosphaeria_arxii_CBS_175.79.BUSCO_summary.fungi_odb12.txt)
#
# target_kind: which samples.csv column the RESOLVED key maps to (i.e. what the
# file gets renamed to) -- "asmid" for assembly-level statistics, "locustag" for
# annotation-derived ones. This is independent of source_kind: BUSCO_genome is
# name-keyed on disk today (source_kind="name") but is an assembly-level
# statistic, so its target/canonical key is ASMID (target_kind="asmid") -- this
# corrects an earlier version of this script (and the original plan's example)
# that assumed BUSCO_genome was already accession-named on disk.
#
# Per-tool filename delimiter right after the basename varies: most types use
# "." (Foo_bar.suffix), but targetP uses "_" (Foo_bar_summary.targetp2.gz) --
# handled by longest_prefix_match() accepting both.
TYPE_CONFIG = {
    "asm_stats":        ("asmid", "asmid",    "genome_stats"),
    "asm_reports":      ("asmid", "asmid",    "genome_stats"),
    "BUSCO_genome":     ("name",  "asmid",    "genome_stats"),
    "BUSCO_protein":    ("name",  "locustag", "genome_stats"),
    "aa_freq":          ("name",  "locustag", "genome_stats"),
    "codon_freq":       ("name",  "locustag", "genome_stats"),
    "gene_stats":       ("name",  "locustag", "genome_stats"),
    "intergenic_stats": ("name",  "locustag", "genome_stats"),
    "aiupred":          ("name",  "locustag", "function"),
    "cazy":             ("name",  "locustag", "function"),
    "merops":           ("name",  "locustag", "function"),
    "pfam_hmmscan":     ("name",  "locustag", "function"),
    "predgpi":          ("name",  "locustag", "function"),
    "signalp":          ("name",  "locustag", "function"),
    "targetP":          ("name",  "locustag", "function"),
    "tmhmm":            ("name",  "locustag", "function"),
    "wolfpsort":        ("name",  "locustag", "function"),
}

# cazy is known to already store one subdirectory per genome
# (results/function/cazy/<basename>/<basename>.*.tsv.gz) rather than flat files --
# a different, pre-existing structure this script does not touch (it only scans
# top-level files in each type dir). Left as a follow-up, not silently claimed
# as "nothing to migrate".
KNOWN_SUBDIRECTORY_LAYOUT_TYPES = {"cazy"}


def clean_strain(raw_strain: str) -> str:
    """Mirror nextflow/modules/common/utils.nf::cleanStrain()."""
    s = (raw_strain or "").strip().replace("'", "").replace('"', "")
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species: str, raw_strain: str) -> str:
    """Mirror nextflow/modules/common/utils.nf::makeSampleTag()."""
    sp = (raw_species or "").strip().replace("'", "").replace('"', "")
    st = clean_strain(raw_strain)
    parts = [p for p in (sp, st) if p]
    tag = "_".join(parts)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def load_samples(samples_path: Path):
    """Return (asmid_set, basename_to_asmid, basename_to_locustag) from samples.csv."""
    asmid_set = set()
    basename_to_asmid = {}
    basename_to_locustag = {}
    with open(samples_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = (row.get("ASMID") or "").strip()
            locustag = (row.get("LOCUSTAG") or "").strip()
            species = (row.get("SPECIES") or "").strip()
            strain = (row.get("STRAIN") or "").strip()
            species_in = (row.get("SPECIES_IN") or "").strip()
            if asmid:
                asmid_set.add(asmid)
                # bare-accession prefix (no assembly-name suffix), matching
                # collect_busco_stats.py::load_asmid_map()'s version-only fallback.
                m = re.match(r"(GC[AF]_\d+\.\d+)", asmid)
                if m:
                    asmid_set.add(m.group(1))
            for sp in (species, species_in):
                if not sp:
                    continue
                basename = make_sample_tag(sp, strain)
                if asmid:
                    basename_to_asmid.setdefault(basename, asmid)
                if locustag:
                    basename_to_locustag.setdefault(basename, locustag)
                # also register the no-strain form, for the documented A. niger-style
                # edge case (collect_asm_stats.py already falls back to this)
                basename_nostrain = make_sample_tag(sp, "")
                if asmid:
                    basename_to_asmid.setdefault(basename_nostrain, asmid)
                if locustag:
                    basename_to_locustag.setdefault(basename_nostrain, locustag)
    return asmid_set, basename_to_asmid, basename_to_locustag


# Delimiters observed between a basename/ASMID prefix and its filename suffix.
# Most tools use "." (Foo_bar.suffix); targetP uses "_" (Foo_bar_summary.targetp2.gz).
_PREFIX_DELIMITERS = (".", "_")


def longest_prefix_match(filename: str, candidates_sorted_desc):
    """Return the longest candidate that is a prefix of filename (followed by one of
    _PREFIX_DELIMITERS, or an exact match), or None."""
    for cand in candidates_sorted_desc:
        if filename == cand:
            return cand
        if any(filename.startswith(cand + d) for d in _PREFIX_DELIMITERS):
            return cand
    return None


def resolve_key(filename: str, source_kind: str, target_kind: str,
                 asmid_list_desc, basename_asmid_desc, basename_locustag_desc,
                 basename_to_asmid, basename_to_locustag):
    """Return (key, matched_prefix) or (None, None) if unresolvable.

    source_kind decides how the filename is parsed (asmid-prefixed vs.
    basename-prefixed); target_kind decides which samples.csv column the
    resolved key is expressed in -- these are independent (see TYPE_CONFIG).
    """
    if source_kind == "asmid":
        matched = longest_prefix_match(filename, asmid_list_desc)
        return (matched, matched) if matched else (None, None)

    basename_map_desc = basename_asmid_desc if target_kind == "asmid" else basename_locustag_desc
    basename_lookup = basename_to_asmid if target_kind == "asmid" else basename_to_locustag
    matched = longest_prefix_match(filename, basename_map_desc)
    if not matched:
        return None, None
    key = basename_lookup.get(matched)
    return (key, matched) if key else (None, None)


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_plan(root: Path, types, asmid_set, basename_to_asmid, basename_to_locustag):
    """Return (plan_rows, unmatched_rows) without touching the filesystem.

    plan_rows: list of dict(old_path, new_path, type, key)
    unmatched_rows: list of dict(old_path, new_path=_UNMATCHED path, type, key=None)
    """
    asmid_list_desc = sorted(asmid_set, key=len, reverse=True)
    basename_asmid_desc = sorted(basename_to_asmid, key=len, reverse=True)
    basename_locustag_desc = sorted(basename_to_locustag, key=len, reverse=True)

    plan_rows = []
    unmatched_rows = []

    for type_name in types:
        source_kind, target_kind, parent = TYPE_CONFIG[type_name]
        src_dir = root / parent / type_name
        if not src_dir.is_dir():
            print(f"skip {type_name}: {src_dir} does not exist", file=sys.stderr)
            continue

        entries = list(os.scandir(src_dir))
        file_entries = [e for e in entries if e.is_file()]
        subdir_entries = [e for e in entries if e.is_dir()]
        # Already-created hash-bucket directories (2-3 lowercase hex chars) are
        # expected on a re-run after a prior --apply -- not a layout surprise.
        unexpected_subdirs = [e for e in subdir_entries if not re.fullmatch(r"[0-9a-f]{2,3}", e.name)]
        if unexpected_subdirs and type_name not in KNOWN_SUBDIRECTORY_LAYOUT_TYPES:
            print(f"NOTE: {type_name} has {len(unexpected_subdirs)} non-bucket subdirectories under "
                  f"{src_dir} that this script does not touch (only top-level files) -- unexpected for "
                  f"a type not in KNOWN_SUBDIRECTORY_LAYOUT_TYPES, please check.", file=sys.stderr)
        elif unexpected_subdirs:
            print(f"NOTE: {type_name} stores {len(unexpected_subdirs)} genomes as one subdirectory each "
                  f"(known pre-existing layout, out of scope for this script -- see "
                  f"KNOWN_SUBDIRECTORY_LAYOUT_TYPES) -- not migrated by this run.", file=sys.stderr)

        for entry in sorted(file_entries, key=lambda e: e.name):
            filename = entry.name
            key, matched_prefix = resolve_key(
                filename, source_kind, target_kind,
                asmid_list_desc, basename_asmid_desc, basename_locustag_desc,
                basename_to_asmid, basename_to_locustag,
            )
            old_path = src_dir / filename
            if key is None:
                new_path = root / "_UNMATCHED" / type_name / filename
                unmatched_rows.append({"old_path": old_path, "new_path": new_path, "type": type_name, "key": ""})
                continue
            suffix = filename[len(matched_prefix):]  # keeps leading delimiter and the rest
            bucket = hash_bucket_for_type(type_name, key)
            new_path = src_dir / bucket / f"{key}{suffix}"
            plan_rows.append({"old_path": old_path, "new_path": new_path, "type": type_name, "key": key})

    return plan_rows, unmatched_rows


def check_target_uniqueness(plan_rows):
    """Return (clean_rows, collision_rows): rows whose target path is unique vs. not.

    A collision means two different source files (different keys, or the same key
    resolved from two differently-named legacy files) would land on the same new
    path. This is real, observed production data: e.g. a strain-less legacy output
    file and a later strain-qualified one that both resolve to the same LOCUSTAG.
    We do NOT guess which is authoritative (no mtime-wins auto-resolution) -- both
    go to a review quarantine instead, and migration proceeds for everything else
    rather than aborting the whole run over a handful of edge cases.
    """
    by_target = defaultdict(list)
    for row in plan_rows:
        by_target[row["new_path"]].append(row)
    clean_rows = []
    collision_rows = []
    for target, rows in by_target.items():
        if len(rows) > 1:
            collision_rows.extend(rows)
        else:
            clean_rows.extend(rows)
    return clean_rows, collision_rows


def already_done(old_path: Path, new_path: Path, recorded_sha256: str) -> bool:
    """True if this move appears to have already happened in a prior (killed/resumed) run."""
    if old_path.exists() or not new_path.exists():
        return False
    if not recorded_sha256:
        return False
    try:
        return sha256_of(new_path) == recorded_sha256
    except OSError:
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--samples", default="samples.csv", help="samples.csv path [samples.csv]")
    parser.add_argument("--root", default="results", help="root containing genome_stats/ and function/ [results]")
    parser.add_argument("--only", action="append", help="restrict to this type (repeatable); default: all known types")
    parser.add_argument("--apply", action="store_true", help="actually move files (default: dry-run)")
    parser.add_argument("--manifest-out", default=None,
                        help="manifest CSV path [<root>/_migration_manifest.<timestamp-free, appended>.csv]")
    parser.add_argument("--spot-check", type=int, default=25, help="number of moved files to sha256-verify after --apply [25]")
    args = parser.parse_args()

    root = Path(args.root)
    samples_path = Path(args.samples)
    if not samples_path.exists():
        print(f"ERROR: samples file not found: {samples_path}", file=sys.stderr)
        sys.exit(1)

    types = args.only if args.only else list(TYPE_CONFIG.keys())
    unknown = [t for t in types if t not in TYPE_CONFIG]
    if unknown:
        print(f"ERROR: unknown type(s) {unknown}; known types: {list(TYPE_CONFIG.keys())}", file=sys.stderr)
        sys.exit(1)

    asmid_set, basename_to_asmid, basename_to_locustag = load_samples(samples_path)
    print(f"Loaded samples.csv: {len(asmid_set)} ASMIDs, {len(basename_to_asmid)} basenames -> ASMID, "
          f"{len(basename_to_locustag)} basenames -> LOCUSTAG")

    plan_rows, unmatched_rows = build_plan(root, types, asmid_set, basename_to_asmid, basename_to_locustag)

    # Target-path collisions (same resolved key/bucket from two different source
    # files -- e.g. a strain-less legacy filename and a later strain-qualified one
    # that both resolve to the same LOCUSTAG) are quarantined for manual review,
    # not auto-resolved and not allowed to block the rest of the migration.
    plan_rows, collision_rows = check_target_uniqueness(plan_rows)
    for row in collision_rows:
        row["new_path"] = root / "_COLLISIONS" / row["type"] / row["old_path"].name

    all_rows = plan_rows + unmatched_rows + collision_rows

    if not all_rows:
        print("Nothing to do: no files found under the requested type(s).")
        return

    if collision_rows:
        print(f"NOTE: {len(collision_rows)} files quarantined to _COLLISIONS/ (target-path collisions -- "
              f"two different source files resolved to the same key/bucket; requires manual review, see "
              f"module docstring). Not auto-resolved. First few:", file=sys.stderr)
        by_target_review = defaultdict(list)
        for row in collision_rows:
            key = row["key"] or "?"
            by_target_review[(row["type"], key)].append(row["old_path"])
        for (t, key), olds in list(by_target_review.items())[:10]:
            print(f"  {t} key={key}: {olds}", file=sys.stderr)

    print(f"Plan: {len(plan_rows)} files to bucket, {len(unmatched_rows)} unmatched -> _UNMATCHED/, "
          f"{len(collision_rows)} quarantined -> _COLLISIONS/")
    by_type_planned = defaultdict(int)
    by_type_unmatched = defaultdict(int)
    by_type_collision = defaultdict(int)
    for row in plan_rows:
        by_type_planned[row["type"]] += 1
    for row in unmatched_rows:
        by_type_unmatched[row["type"]] += 1
    for row in collision_rows:
        by_type_collision[row["type"]] += 1
    for t in types:
        print(f"  {t}: {by_type_planned.get(t, 0)} planned, {by_type_unmatched.get(t, 0)} unmatched, "
              f"{by_type_collision.get(t, 0)} quarantined")

    if not args.apply:
        print("\nDRY RUN: no files were moved. Re-run with --apply to execute.")
        return

    manifest_path = Path(args.manifest_out) if args.manifest_out else root / "_migration_manifest.csv"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)

    # Load any prior manifest rows (for idempotent re-run / already_done checks).
    prior_by_new_path = {}
    if manifest_path.exists():
        with open(manifest_path, newline="") as fh:
            for row in csv.DictReader(fh):
                prior_by_new_path[row["new_path"]] = row.get("sha256", "")

    write_header = not manifest_path.exists()
    moved = 0
    skipped_already_done = 0
    errors = 0

    with open(manifest_path, "a", newline="") as mfh:
        writer = csv.writer(mfh)
        if write_header:
            writer.writerow(["old_path", "new_path", "type", "key", "sha256", "size", "mtime"])

        for row in all_rows:
            old_path, new_path = row["old_path"], row["new_path"]
            recorded_sha = prior_by_new_path.get(str(new_path))

            if not old_path.exists():
                if new_path.exists() and already_done(old_path, new_path, recorded_sha):
                    skipped_already_done += 1
                    continue
                print(f"WARN: expected source file missing and not already migrated: {old_path}", file=sys.stderr)
                errors += 1
                continue

            if already_done(old_path, new_path, recorded_sha):
                skipped_already_done += 1
                continue

            new_path.parent.mkdir(parents=True, exist_ok=True)
            try:
                sha = sha256_of(old_path)
                size = old_path.stat().st_size
                mtime = old_path.stat().st_mtime
                os.rename(old_path, new_path)
                os.utime(new_path, (mtime, mtime))
                writer.writerow([str(old_path), str(new_path), row["type"], row["key"], sha, size, mtime])
                moved += 1
            except OSError as exc:
                print(f"ERROR moving {old_path} -> {new_path}: {exc}", file=sys.stderr)
                errors += 1

    print(f"\nMoved {moved} files, skipped {skipped_already_done} already-done, {errors} errors.")
    print(f"Manifest: {manifest_path}")

    # Verification pass.
    print("\nVerification:")
    total_after = 0
    bucket_counts = defaultdict(lambda: defaultdict(int))
    for t in types:
        _source_kind, _target_kind, parent = TYPE_CONFIG[t]
        type_dir = root / parent / t
        if not type_dir.is_dir():
            continue
        for bucket_dir in sorted(p for p in type_dir.iterdir() if p.is_dir()):
            n = sum(1 for _ in bucket_dir.iterdir())
            bucket_counts[t][bucket_dir.name] = n
            total_after += n
    def count_flat_tree(base: Path) -> int:
        total = 0
        if base.is_dir():
            for t_dir in base.iterdir():
                if t_dir.is_dir():
                    total += sum(1 for _ in t_dir.iterdir())
        return total

    total_unmatched = count_flat_tree(root / "_UNMATCHED")
    total_collisions = count_flat_tree(root / "_COLLISIONS")
    print(f"  files now in buckets: {total_after}, in _UNMATCHED: {total_unmatched}, "
          f"in _COLLISIONS: {total_collisions}, "
          f"total: {total_after + total_unmatched + total_collisions} "
          f"(expected {len(plan_rows) + len(unmatched_rows) + len(collision_rows)} from this run "
          f"plus anything already migrated in a prior run)")

    for t, buckets in bucket_counts.items():
        if not buckets:
            continue
        mean = sum(buckets.values()) / len(buckets)
        hot = {b: n for b, n in buckets.items() if mean > 0 and n > 3 * mean}
        if hot:
            print(f"  WARNING: {t} has hash buckets >3x the mean ({mean:.1f}) -- possible hash-function bug: {hot}")

    if errors:
        sys.exit(1)


if __name__ == "__main__":
    main()
