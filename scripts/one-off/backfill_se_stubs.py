#!/usr/bin/env python3
"""Backfill zero-byte `_norm_SE.fastq.gz` stubs for genuine paired-end species.

Background
----------
`SRA_FETCH` (nextflow/funannotate.nf) declares a 3-file output tuple per species:
`_norm_R1`, `_norm_R2`, `_norm_SE` (all in `rnaseq_reads/`, via `storeDir`). The
`_norm_SE` output was added on 2026-06-12 (commit f7db039). Species whose reads were
built before that date have R1/R2 but no SE file, so their `storeDir` is incomplete and
Nextflow re-runs the expensive `SRA_FETCH` (32 cpu / 96 GB / ~2 h, re-downloading and
re-normalizing hundreds of MB) purely to regenerate outputs that already exist.

`SRA_FETCH` itself creates the SE output as a zero-byte file for paired-end species
(`: > <tag>_norm_SE.fastq.gz`), and most existing on-disk SE files are zero bytes. So the
correct completion for a paired-end species is simply an empty SE stub.

Safety: only genuine PE species are stubbed
-------------------------------------------
A stub makes `storeDir` consider the species done, which also stops SRA_FETCH_SE from
ever running for it. That is correct for a species that has real paired data (it will
never route to SE fetch anyway), but WRONG for a species with empty R1/R2 (no PE data
found) — those must stay eligible to pick up single-end data from a re-query + SE fetch.
This script therefore stubs ONLY species where both R1 and R2 are non-empty and SE is
missing. Species with empty R1/R2 (and missing SE) are reported and left untouched.

Usage
-----
    python3 scripts/one-off/backfill_se_stubs.py            # dry-run (default)
    python3 scripts/one-off/backfill_se_stubs.py --apply    # create the stubs
    python3 scripts/one-off/backfill_se_stubs.py --dir rnaseq_reads --limit 10
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

R1_SUFFIX = "_norm_R1.fastq.gz"


def classify(reads_dir: str):
    """Return (to_stub, left_empty, left_partial) lists of species tags / info."""
    to_stub, left_empty, left_partial = [], [], []
    for r1 in glob.glob(os.path.join(reads_dir, f"*{R1_SUFFIX}")):
        tag = os.path.basename(r1)[: -len(R1_SUFFIX)]
        r2 = os.path.join(reads_dir, f"{tag}_norm_R2.fastq.gz")
        se = os.path.join(reads_dir, f"{tag}_norm_SE.fastq.gz")
        if os.path.exists(se):
            continue  # store already complete
        r1sz = os.path.getsize(r1)
        r2sz = os.path.getsize(r2) if os.path.exists(r2) else -1
        if r1sz > 0 and r2sz > 0:
            to_stub.append((tag, se))
        elif r1sz == 0 and r2sz == 0:
            left_empty.append(tag)
        else:
            left_partial.append((tag, r1sz, r2sz))
    return to_stub, left_empty, left_partial


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="rnaseq_reads",
                    help="directory of *_norm_*.fastq.gz files (default: %(default)s)")
    ap.add_argument("--apply", action="store_true",
                    help="create the zero-byte SE stubs (default: dry-run)")
    ap.add_argument("--limit", type=int, default=0,
                    help="stub at most N species (0 = all; for testing)")
    args = ap.parse_args()

    if not os.path.isdir(args.dir):
        print(f"[ERROR] not a directory: {args.dir}", file=sys.stderr)
        return 1

    to_stub, left_empty, left_partial = classify(args.dir)
    if args.limit:
        to_stub = to_stub[: args.limit]

    created = 0
    if args.apply:
        for tag, se in to_stub:
            # Match SRA_FETCH's `: > file`: create an empty file, do not clobber a
            # non-empty one that may have appeared since classification.
            if os.path.exists(se) and os.path.getsize(se) > 0:
                continue
            with open(se, "w"):
                pass
            created += 1

    mode = "APPLIED" if args.apply else "DRY-RUN (no files created)"
    print(f"=== backfill_se_stubs [{mode}] ===")
    print(f"reads dir: {args.dir}")
    print(f"stub targets (real PE, SE missing): {len(to_stub)}")
    if args.apply:
        print(f"stubs created:                      {created}")
    print(f"left untouched (empty R1/R2):       {len(left_empty)}   "
          f"<- eligible for SE fetch after re-query")
    if left_partial:
        print(f"\n[WARN] {len(left_partial)} species with exactly one empty mate "
              f"(not stubbed, inspect manually):")
        for tag, a, b in left_partial[:20]:
            print(f"   {tag}: R1={a}B R2={b}B")
    if not args.apply and to_stub:
        print("\nRe-run with --apply to create these stubs.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
