#!/usr/bin/env python3
"""One-time migration: move results/function/cazy/ from one-subdirectory-per-genome
to the LOCUSTAG hash-bucket layout every other function/* type already uses.

    before  results/function/cazy/<out>/<out>.overview.tsv.gz
    after   results/function/cazy/<hash2>/<LOCUSTAG>.overview.tsv.gz

Why
---
`out` is makeSampleTag(SPECIES, STRAIN). It is derived, not stable: a taxonomic
rename rebuilds it, so the next run writes to a NEW directory and the old one is
left behind for ever. The merge glob ("cazy/*/*.overview.tsv.gz") reads whatever
directories exist, so an orphan keeps being merged alongside its replacement --
two copies of one genome under one LOCUSTAG in cazy.overview.parquet. Measured on
2026-11-21: 24 genomes present under two names at once, plus 8 directories for
genomes no longer in samples.csv at all. LOCUSTAG does not move when taxonomy
does, so after this migration a rename rewrites the same path and the failure
cannot recur. Plan: todo/genome_stats_storage_reorg.md §A, which already assigns
LOCUSTAG and width 2 to all of results/function/*.

storeDir, not publishDir
------------------------
nextflow/modules/BFD/CAZY/main.nf declares storeDir, which is the cache: Nextflow
skips the process when the declared outputs already exist there. This migration
and that module must be applied TOGETHER. Either alone leaves storeDir pointing
at a path with no files and re-runs dbCAN for every genome in samples.csv.

What is NOT migrated
--------------------
Nothing is deleted. Directories that must not become the bucketed copy are moved
whole into a quarantine sibling, `cazy.superseded_<DATE>/`, which no glob reads:

  DUPLICATE_LOSER  another directory holds the same LOCUSTAG and carries the
                   current `out` name (verified newer in all 24 live cases)
  PHANTOM          the LOCUSTAG is not in samples.csv at all
  EMPTY            the directory has no overview file to migrate
  UNRESOLVED       the LOCUSTAG could not be determined -- never guessed

Resuming
--------
The plan is rebuilt from what is on disk each run and every move is skipped when
its target already exists, so re-running after an interruption continues rather
than corrupting. --dry-run writes the manifest and changes nothing.

Usage
-----
    migrate_cazy_to_hash_buckets.py --runs-root /path/to/Fungi_BFD_runs --dry-run
    migrate_cazy_to_hash_buckets.py --runs-root /path/to/Fungi_BFD_runs --apply
"""
import argparse
import csv
import datetime
import gzip
import importlib.util
import os
import shutil
import sys

TYPES = ("overview", "cazymes", "substrates")
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


# Imported, never re-implemented: a second copy of either of these drifting from
# the original is exactly the silent path miss this migration exists to remove.
_samples = _load(os.path.join(REPO, "scripts", "create_samples_file.py"), "_cazy_samples")
_paths = _load(os.path.join(REPO, "nextflow", "bin", "genome_stats_paths.py"), "_cazy_paths")
make_sample_tag = _samples.make_sample_tag
hash_bucket_for_type = _paths.hash_bucket_for_type


def load_samples(samples_csv):
    """LOCUSTAG -> current `out`, for every row of samples.csv."""
    out_by_lt = {}
    with open(samples_csv, newline="") as fh:
        for row in csv.DictReader(fh):
            lt = (row.get("LOCUSTAG") or "").strip()
            if lt:
                out_by_lt[lt] = make_sample_tag(row.get("SPECIES", ""), row.get("STRAIN", ""))
    return out_by_lt


def locustag_from_overview(path):
    """Read the LOCUSTAG out of an overview.tsv.gz. Gene IDs are '<LOCUSTAG>_<n>',
    so the first data row names the genome. Returns None rather than guessing."""
    try:
        with gzip.open(path, "rt") as fh:
            next(fh, None)                      # header
            first = next(fh, "").strip()
        return first.split("\t")[0].split("_")[0] or None if first else None
    except (OSError, EOFError, StopIteration):
        return None


def build_plan(cazy_dir, out_by_lt):
    """One row per directory: (dirname, locustag, action, bucket, reason)."""
    lt_by_out = {out: lt for lt, out in out_by_lt.items()}
    dirs = sorted(d for d in os.listdir(cazy_dir) if os.path.isdir(os.path.join(cazy_dir, d)))
    # Buckets are 2 hex chars; a bucket directory is never a genome directory.
    width = _paths.bucket_width_for("cazy")
    dirs = [d for d in dirs if not (len(d) == width and all(c in "0123456789abcdef" for c in d))]

    resolved = {}
    for d in dirs:
        lt = lt_by_out.get(d)
        if lt is None:
            lt = locustag_from_overview(os.path.join(cazy_dir, d, d + ".overview.tsv.gz"))
        resolved[d] = lt

    # A LOCUSTAG claimed by several directories: the one carrying the current
    # `out` name wins; failing that the newest by mtime. Both agreed in all 24
    # live cases on 2026-11-21, so this is a tie-break, not a judgement call.
    claims = {}
    for d, lt in resolved.items():
        if lt:
            claims.setdefault(lt, []).append(d)
    winner = {}
    for lt, ds in claims.items():
        if len(ds) == 1:
            winner[lt] = ds[0]
            continue
        cur = out_by_lt.get(lt)
        if cur in ds:
            winner[lt] = cur
        else:
            winner[lt] = max(ds, key=lambda x: os.path.getmtime(os.path.join(cazy_dir, x)))

    plan = []
    for d in dirs:
        lt = resolved[d]
        if lt is None:
            has_any = any(os.path.exists(os.path.join(cazy_dir, d, "%s.%s.tsv.gz" % (d, t))) for t in TYPES)
            plan.append((d, "", "EMPTY" if not has_any else "UNRESOLVED", "",
                         "no overview file" if not has_any else "LOCUSTAG not determinable"))
        elif lt not in out_by_lt:
            plan.append((d, lt, "PHANTOM", "", "LOCUSTAG not in samples.csv"))
        elif winner[lt] != d:
            plan.append((d, lt, "DUPLICATE_LOSER", "", "superseded by %s" % winner[lt]))
        else:
            plan.append((d, lt, "MIGRATE", hash_bucket_for_type("cazy", lt), ""))
    return plan


def write_manifest(plan, path):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["DIR", "LOCUSTAG", "ACTION", "BUCKET", "REASON"])
        w.writerows(plan)


def apply_plan(plan, cazy_dir, quarantine):
    moved = quarantined = skipped = 0
    for d, lt, action, bucket, _reason in plan:
        src_dir = os.path.join(cazy_dir, d)
        if not os.path.isdir(src_dir):
            skipped += 1                        # already handled by an earlier run
            continue
        if action == "MIGRATE":
            dst_dir = os.path.join(cazy_dir, bucket)
            os.makedirs(dst_dir, exist_ok=True)
            for t in TYPES:
                src = os.path.join(src_dir, "%s.%s.tsv.gz" % (d, t))
                dst = os.path.join(dst_dir, "%s.%s.tsv.gz" % (lt, t))
                if os.path.exists(dst) or not os.path.exists(src):
                    continue
                os.rename(src, dst)
            leftovers = os.listdir(src_dir)
            if leftovers:
                print("  KEEPING %s: %d unexpected file(s) left: %s"
                      % (d, len(leftovers), ", ".join(sorted(leftovers)[:4])), file=sys.stderr)
            else:
                os.rmdir(src_dir)
                moved += 1
        else:
            os.makedirs(quarantine, exist_ok=True)
            dst = os.path.join(quarantine, d)
            if os.path.exists(dst):
                skipped += 1
                continue
            shutil.move(src_dir, dst)
            quarantined += 1
    return moved, quarantined, skipped


def verify(plan, cazy_dir):
    """Every MIGRATE row must have its overview file at the bucketed path, and no
    planned source directory may survive under cazy/."""
    missing, survivors = [], []
    for d, lt, action, bucket, _r in plan:
        if action == "MIGRATE":
            if not os.path.exists(os.path.join(cazy_dir, bucket, "%s.overview.tsv.gz" % lt)):
                missing.append("%s -> %s/%s" % (d, bucket, lt))
        if os.path.isdir(os.path.join(cazy_dir, d)):
            survivors.append(d)
    return missing, survivors


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs-root", required=True, help="Fungi_BFD_runs checkout")
    ap.add_argument("--manifest", default=None, help="where to write the plan (default: alongside cazy/)")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--dry-run", action="store_true", help="write the manifest, change nothing")
    g.add_argument("--apply", action="store_true", help="carry out the plan")
    a = ap.parse_args()

    root = os.path.abspath(a.runs_root)
    cazy_dir = os.path.join(root, "results", "function", "cazy")
    samples_csv = os.path.join(root, "samples.csv")
    for p in (cazy_dir, samples_csv):
        if not os.path.exists(p):
            sys.exit("ERROR: missing %s" % p)

    stamp = datetime.date.today().strftime("%Y%m%d")
    quarantine = os.path.join(root, "results", "function", "cazy.superseded_%s" % stamp)
    manifest = a.manifest or os.path.join(root, "results", "function", "cazy_migration_%s.tsv" % stamp)

    print("cazy dir   : %s" % cazy_dir)
    print("samples.csv: %s" % samples_csv)
    print("quarantine : %s" % quarantine)

    out_by_lt = load_samples(samples_csv)
    print("samples.csv rows with a LOCUSTAG: %d" % len(out_by_lt))

    plan = build_plan(cazy_dir, out_by_lt)
    counts = {}
    for _d, _lt, action, _b, _r in plan:
        counts[action] = counts.get(action, 0) + 1
    print("\nplan over %d directories:" % len(plan))
    for k in ("MIGRATE", "DUPLICATE_LOSER", "PHANTOM", "EMPTY", "UNRESOLVED"):
        print("  %-16s %d" % (k, counts.get(k, 0)))

    write_manifest(plan, manifest)
    print("\nmanifest: %s" % manifest)

    if counts.get("UNRESOLVED"):
        print("\nUNRESOLVED directories (quarantined, never guessed):", file=sys.stderr)
        for d, _lt, action, _b, r in plan:
            if action == "UNRESOLVED":
                print("  %s  (%s)" % (d, r), file=sys.stderr)

    if a.dry_run:
        print("\ndry run: nothing moved")
        return 0

    moved, quarantined, skipped = apply_plan(plan, cazy_dir, quarantine)
    print("\nmigrated %d, quarantined %d, skipped %d (already done)" % (moved, quarantined, skipped))

    missing, survivors = verify(plan, cazy_dir)
    if missing:
        print("\nVERIFY FAILED: %d migrated genome(s) have no file at the target path:" % len(missing),
              file=sys.stderr)
        for m in missing[:20]:
            print("  %s" % m, file=sys.stderr)
        return 1
    if survivors:
        print("\nVERIFY: %d source directory/ies still present (re-run to finish):" % len(survivors),
              file=sys.stderr)
        for s in survivors[:20]:
            print("  %s" % s, file=sys.stderr)
        return 1
    print("verify: every migrated genome is readable at its bucketed path; no source dirs remain")
    return 0


if __name__ == "__main__":
    sys.exit(main())
