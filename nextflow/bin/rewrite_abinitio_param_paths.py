#!/usr/bin/env python3
"""
rewrite_abinitio_param_paths.py — rewrite the `path` entries in the shared
ab-initio parameter stores' parameters.json files between absolute and relative
form, validating that every referenced component actually exists.

Why
---
`funannotate predict -p <parameters.json>` reads a `path` per predictor
(augustus/genemark/snap/glimmerhmm/codingquarry). Historically those were
absolute, which makes a store non-relocatable: move or copy the tree and every
JSON silently points at the old location. Worse, stock funannotate hands a
missing path to lib.copyDirectory(), which swallows the error and leaves an
EMPTY parameter dir behind -- the prediction quietly degrades to untrained
instead of failing.

  --mode relative   Each component is stored beside its parameters.json, so the
                    path becomes a bare filename and the store relocates freely.
                    REQUIRES the patched funannotate that resolves a relative
                    path against the parameters.json directory
                    (resolveTrainingPaths() in funannotate/predict.py). Stock
                    funannotate resolves against the process CWD and will NOT
                    find these.
  --mode absolute   The reverse: rewrite to absolute paths rooted at the store's
                    real on-disk location. Use to repair stores whose absolute
                    paths went stale after a move, or to run against a stock
                    funannotate install.

mtime is preserved by default
-----------------------------
staleSharedParams() in nextflow/modules/funannotate/utils.nf treats a
parameters.json newer than a sibling strain's GBK as "re-predict me". A bulk
rewrite that bumped every mtime would therefore schedule a full re-prediction of
every reuse_eligible genome in the project. Since this rewrite changes only how
a path is spelled and not which trained parameters it names, the original mtime
is restored after writing. Pass --touch if you genuinely want the cascade.

Usage
-----
    # inspect first -- always
    rewrite_abinitio_param_paths.py --shared-root gene_prediction_shared_abinitio \
        --mode relative --dry-run

    # apply
    rewrite_abinitio_param_paths.py --shared-root gene_prediction_shared_abinitio \
        --mode relative
"""

import argparse
import json
import os
import sys
from pathlib import Path

PREDICTORS = ["augustus", "genemark", "codingquarry", "snap", "glimmerhmm"]


def resolve_existing(given: str, store_dir: Path):
    """Find the component `given` names, trying the path as-written, then the
    same basename inside store_dir (which repairs a stale absolute path left
    behind by a relocated store). Returns an existing Path, or None."""
    candidates = []
    p = Path(given)
    if p.is_absolute():
        candidates.append(p)
    else:
        candidates.append(store_dir / given)
        candidates.append(Path(os.path.abspath(given)))
    # basename fallback: the component lives beside the JSON even though the
    # recorded path points somewhere else entirely
    candidates.append(store_dir / p.name)
    for c in candidates:
        if c.exists():
            return c
    return None


def rewrite_one(json_path: Path, mode: str, dry_run: bool, touch: bool):
    """Returns (changed: bool, missing: list[(predictor, given)])."""
    store_dir = json_path.parent
    try:
        data = json.loads(json_path.read_text())
    except (OSError, ValueError) as e:
        print(f"[ERROR] {json_path}: unreadable ({e})", file=sys.stderr)
        return False, [("<file>", str(e))]

    changed = False
    missing = []
    for predictor in PREDICTORS:
        entries = data.get(predictor)
        if not entries or not isinstance(entries, list):
            continue
        given = entries[0].get("path") if isinstance(entries[0], dict) else None
        if not given:
            continue
        found = resolve_existing(given, store_dir)
        if found is None:
            missing.append((predictor, given))
            continue
        if mode == "relative":
            try:
                new = os.path.relpath(found.resolve(), store_dir.resolve())
            except ValueError:  # different drive/mount; leave absolute
                new = str(found.resolve())
        else:
            new = str(found.resolve())
        if new != given:
            entries[0]["path"] = new
            changed = True

    if missing:
        for predictor, given in missing:
            print(f"[WARN] {json_path}: {predictor} component not found: {given}",
                  file=sys.stderr)

    if changed and not dry_run:
        st = json_path.stat()
        tmp = json_path.with_suffix(".json.tmp.%d" % os.getpid())
        tmp.write_text(json.dumps(data, indent=2))
        os.replace(tmp, json_path)
        if not touch:
            os.utime(json_path, (st.st_atime, st.st_mtime))
    return changed, missing


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--shared-root", required=True, type=Path,
                    help="gene_prediction_shared_abinitio directory")
    ap.add_argument("--mode", required=True, choices=["relative", "absolute"])
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--touch", action="store_true",
                    help="do NOT preserve mtime; lets staleSharedParams() "
                         "schedule a re-predict of every reuse_eligible strain")
    args = ap.parse_args()

    if not args.shared_root.is_dir():
        sys.exit(f"[ERROR] not a directory: {args.shared_root}")

    jsons = sorted(args.shared_root.glob("*/parameters.json"))
    if not jsons:
        sys.exit(f"[ERROR] no */parameters.json under {args.shared_root}")

    changed = unchanged = 0
    incomplete = []
    for j in jsons:
        did, missing = rewrite_one(j, args.mode, args.dry_run, args.touch)
        if did:
            changed += 1
        else:
            unchanged += 1
        if missing:
            incomplete.append((j, missing))

    verb = "would rewrite" if args.dry_run else "rewrote"
    print(f"\n{verb} {changed} / {len(jsons)} parameters.json to {args.mode} paths "
          f"({unchanged} already correct)")
    if not args.dry_run and changed and not args.touch:
        print("mtimes preserved -- no re-prediction cascade will be triggered")
    if incomplete:
        print(f"[WARN] {len(incomplete)} store(s) reference components that do not "
              f"exist on disk; those entries were left as-is. With the patched "
              f"funannotate these now warn at predict time and fall back to "
              f"training from scratch rather than silently producing empty "
              f"parameter dirs.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
