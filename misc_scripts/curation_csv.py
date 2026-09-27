#!/usr/bin/env python3
"""Append rows to, and check, the hand-curated RNA-seq CSV files.

Why this exists
---------------
rnaseq_skip.csv, rnaseq_force_no_rnaseq.csv and similar files carry a
free-text reason or evidence column. Rows were appended by hand with unquoted
commas inside that text, so a real CSV parser saw 3-7 columns instead of 2 or 4
(found 2026-09-27: all 7 rnaseq_skip.csv rows, 2 of 4
rnaseq_force_no_rnaseq.csv rows). The Nextflow loaders only split on ',' and
read the leading columns, so they kept working, but csv.DictReader readers
(e.g. Funannotate_benchmarking/scripts/tag_trinity_failures.py) got shifted
fields. Use `append` instead of echo/printf so text fields are quoted, and
`check` after any manual edit.

Quoting keeps the leading key columns unquoted, so the Nextflow split(',')
loaders still read the same values.

Usage
-----
  curation_csv.py append rnaseq_skip.csv out=Foo_bar_X1 reason="text, with commas"
  curation_csv.py append rnaseq_force_no_rnaseq.csv species_tag=Foo_bar \\
      reason="..." date=2026-09-27 evidence="..."
  curation_csv.py check rnaseq_skip.csv rnaseq_force_no_rnaseq.csv
  curation_csv.py check --columns 4 --allow-extra rnaseq_blacklist.csv

Files whose first line starts with '#' (rnaseq_blacklist.csv) have no header
row; pass --columns for them. '#' comment lines are skipped.
"""
import argparse
import csv
import io
import sys
from pathlib import Path


def read_rows(path: Path):
    """Return (header or None, [(line_no, row), ...]) skipping blank/# lines."""
    text = path.read_text()
    rows = []
    header = None
    for n, row in enumerate(csv.reader(io.StringIO(text)), start=1):
        if not row or not "".join(row).strip() or row[0].lstrip().startswith("#"):
            continue
        if header is None and n == 1:
            header = row
            continue
        rows.append((n, row))
    return header, rows


def check(paths, columns=None, allow_extra=False) -> int:
    bad = 0
    for p in paths:
        path = Path(p)
        header, rows = read_rows(path)
        want = columns or (len(header) if header else None)
        if want is None:
            print(f"[ERROR] {path}: no header row; pass --columns", file=sys.stderr)
            bad += 1
            continue
        errs = [(n, len(r)) for n, r in rows
                if len(r) < want or (len(r) > want and not allow_extra)]
        for n, k in errs:
            print(f"[ERROR] {path}:{n}: {k} columns, expected {want}", file=sys.stderr)
        bad += len(errs)
        print(f"{path}: {len(rows)} rows, {len(errs)} bad (expected {want} columns)")
    return 1 if bad else 0


def append(path: Path, pairs) -> int:
    header, rows = read_rows(path)
    if header is None:
        print(f"[ERROR] {path}: no header row; append by hand and run check", file=sys.stderr)
        return 1
    given = {}
    for p in pairs:
        if "=" not in p:
            print(f"[ERROR] expected column=value, got {p!r}", file=sys.stderr)
            return 1
        k, v = p.split("=", 1)
        given[k] = v
    missing = [c for c in header if c not in given]
    unknown = [k for k in given if k not in header]
    if missing or unknown:
        print(f"[ERROR] {path} columns are {header}; missing {missing}, unknown {unknown}",
              file=sys.stderr)
        return 1
    if any("\n" in v for v in given.values()):
        print("[ERROR] values must not contain newlines", file=sys.stderr)
        return 1
    text = path.read_text()
    with open(path, "a", newline="") as fh:
        if text and not text.endswith("\n"):
            fh.write("\n")
        csv.writer(fh, lineterminator="\n").writerow([given[c] for c in header])
    return check([path])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("append", help="append one row (column=value ...)")
    a.add_argument("file", type=Path)
    a.add_argument("pairs", nargs="+")
    c = sub.add_parser("check", help="report rows with the wrong column count")
    c.add_argument("files", nargs="+")
    c.add_argument("--columns", type=int, help="expected columns when the file has no header")
    c.add_argument("--allow-extra", action="store_true",
                   help="accept trailing extra columns (rnaseq_blacklist.csv info column)")
    args = ap.parse_args()
    if args.cmd == "append":
        sys.exit(append(args.file, args.pairs))
    sys.exit(check(args.files, args.columns, args.allow_extra))


if __name__ == "__main__":
    main()
