#!/usr/bin/env python3
"""Quarantine empty (header-only) SRA-query CSVs so SRA_QUERY re-searches them.

Background
----------
`SRA_QUERY` / `SRA_QUERY_BATCH` (nextflow/funannotate.nf) cache one CSV per species in
`rnaseq_reads/sra_query/`. Both treat an existing cached CSV as authoritative: SRA_QUERY
uses `storeDir` (skips the process if the output file exists) and SRA_QUERY_BATCH reuses
any file that passes `[ -s "$cached" ]`. A *header-only* CSV (schema header, zero data
rows) means "this species was queried and returned nothing" — and because the header
line makes the file non-empty, it is reused forever and never re-queried.

Many of those empty results are stale: they were produced by the old PAIRED[Layout]-only
query, before the SINGLE[Layout] fallback existed, so the species was never checked for
single-end RNA-seq at all. This script finds the header-only CSVs and moves them into a
quarantine directory (reversible), so the next pipeline run re-queries those species.

After running all data rows are 6-column with no empty fields (see
refresh_sra_query_layout.py), so "empty" here means header-only. As a safety net the
script also flags any file that is not 6-column or has a data row with an empty field.

IMPORTANT — enable_single_end
-----------------------------
The single-end fallback only runs when `params.enable_single_end == true`. In
conf/profile_funannotate.config it is currently **false**, so a bare re-query re-runs
only the PAIRED search (useful only if new paired data was deposited since the original
query). To actually discover single-end data for these species, set
`enable_single_end = true` (and `max_rnaseq_se_runs` as desired) before re-running.

Trigger mechanism
-----------------
Removing a species' CSV from `rnaseq_reads/sra_query/` makes both the storeDir check and
the batch reuse check miss, so the query re-runs. This script *moves* (not deletes) the
files to a timestamped quarantine dir; restore with `mv` if needed. Worst case a re-query
fails and re-writes a header-only CSV — no data is lost, since these files had none.

Usage
-----
    python3 scripts/one-off/requery_empty_sra_query.py                 # dry-run
    python3 scripts/one-off/requery_empty_sra_query.py --apply         # quarantine them
    python3 scripts/one-off/requery_empty_sra_query.py --apply --quarantine-dir <path>
    python3 scripts/one-off/requery_empty_sra_query.py --list species_to_requery.txt
"""
from __future__ import annotations

import argparse
import datetime as dt
import glob
import os
import shutil
import sys

EXPECTED_COLS = 6


def scan(qdir: str):
    """Classify CSVs. Returns (header_only, malformed) as lists of paths."""
    header_only, malformed = [], []
    for path in sorted(glob.glob(os.path.join(qdir, "*.sra_query.csv"))):
        lines = open(path).read().splitlines()
        if not lines:
            malformed.append((path, "empty-file"))
            continue
        if lines[0].count(",") + 1 != EXPECTED_COLS:
            malformed.append((path, f"{lines[0].count(',') + 1}col-header"))
        data = [l for l in lines[1:] if l.strip()]
        if not data:
            header_only.append(path)
            continue
        for l in data:
            fields = l.split(",")
            if len(fields) != EXPECTED_COLS or any(x.strip() == "" for x in fields):
                malformed.append((path, "empty-field-in-data-row"))
                break
    return header_only, malformed


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dir", default="rnaseq_reads/sra_query",
                    help="sra_query cache dir (default: %(default)s)")
    ap.add_argument("--apply", action="store_true",
                    help="move header-only CSVs to the quarantine dir (default: dry-run)")
    ap.add_argument("--quarantine-dir", default=None,
                    help="destination for quarantined CSVs "
                         "(default: <dir>/../sra_query_requery_quarantine_<UTCdate>)")
    ap.add_argument("--list", metavar="FILE", default=None,
                    help="also write the species tags to re-query to FILE")
    ap.add_argument("--limit", type=int, default=0,
                    help="quarantine at most N files (0 = all; for testing)")
    args = ap.parse_args()

    if not os.path.isdir(args.dir):
        print(f"[ERROR] not a directory: {args.dir}", file=sys.stderr)
        return 1

    header_only, malformed = scan(args.dir)
    targets = header_only[: args.limit] if args.limit else header_only

    qdir = args.quarantine_dir or os.path.join(
        os.path.dirname(os.path.abspath(args.dir)),
        f"sra_query_requery_quarantine_{dt.datetime.utcnow():%Y%m%d}")

    moved = 0
    if args.apply and targets:
        os.makedirs(qdir, exist_ok=True)
        for path in targets:
            dest = os.path.join(qdir, os.path.basename(path))
            shutil.move(path, dest)
            moved += 1

    if args.list:
        tags = [os.path.basename(p)[: -len(".sra_query.csv")] for p in targets]
        with open(args.list, "w") as fh:
            fh.write("\n".join(tags) + ("\n" if tags else ""))

    mode = "APPLIED" if args.apply else "DRY-RUN (no files moved)"
    print(f"=== requery_empty_sra_query [{mode}] ===")
    print(f"cache dir:                {args.dir}")
    print(f"header-only (re-query):   {len(header_only)}")
    print(f"malformed (flagged only): {len(malformed)}")
    if args.apply:
        print(f"moved to quarantine:      {moved}")
        print(f"quarantine dir:           {qdir}")
    if args.list:
        print(f"species list written to:  {args.list}")
    if malformed:
        print("\n[WARN] malformed files (NOT moved; inspect):")
        for path, why in malformed[:20]:
            print(f"   {why:24s} {os.path.basename(path)}")
        if len(malformed) > 20:
            print(f"   ... and {len(malformed) - 20} more")
    print("\nNOTE: the single-end fallback only fires when params.enable_single_end=true")
    print("      (set in conf/profile_funannotate.config). Confirm it is true before")
    print("      re-running if the goal is to discover single-end RNA-seq.")
    if not args.apply and header_only:
        print("\nRe-run with --apply to quarantine these and trigger re-query.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
