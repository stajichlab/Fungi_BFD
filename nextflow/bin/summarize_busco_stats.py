#!/usr/bin/env python3
"""Parse BUSCO_GENOME `*.BUSCO_summary.*.txt` reports, write a merged TSV.

Inputs are supplied either as a manifest file (--manifest, one TAB-delimited
`<path>\\t<mtime>\\t<size>` record per line, as written by BFD.nf toManifest) or
as a directory to glob (--reportdir). Output is written to -o; if it ends in
.gz it is gzip-compressed.

BUSCO_GENOME is an assembly-level statistic (nextflow/modules/BFD/BUSCO_GENOME/
main.nf writes {bucket}/{ASMID}.BUSCO_summary.{lineage}.txt), so the key IS the
ASMID directly -- no samples.csv join needed, unlike summarize_asm_stats.py.
"""

import argparse
import csv
import glob
import gzip
import os
import re
import sys

_COMPLETE_RE = re.compile(
    r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)")

COLUMNS = ["ASMID", "complete_pct", "single_pct", "duplicated_pct",
           "fragmented_pct", "missing_pct", "n_markers", "lineage"]


def parse_summary(path):
    """Parse a BUSCO_summary.*.txt file; return dict of metrics or None."""
    text = open(path, errors="replace").read()
    m = _COMPLETE_RE.search(text)
    if not m:
        return None
    return {
        "complete_pct": float(m.group(1)),
        "single_pct": float(m.group(2)),
        "duplicated_pct": float(m.group(3)),
        "fragmented_pct": float(m.group(4)),
        "missing_pct": float(m.group(5)),
        "n_markers": int(m.group(6)),
    }


def read_manifest(path):
    """Read a toManifest file; return the list of report paths (first TAB field)."""
    files = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            files.append(line.split("\t", 1)[0])
    return files


def collect_report_files(args):
    """Resolve the list of BUSCO_summary.*.txt files from --manifest or --reportdir."""
    if args.manifest:
        files = read_manifest(args.manifest)
    else:
        files = glob.glob(os.path.join(args.reportdir, "*.BUSCO_summary.*.txt"))
    return sorted(files)


def open_out(path):
    """Open path for text writing, gzip-compressed if it ends in .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--manifest", help="TAB-delimited manifest of report paths (path\\tmtime\\tsize)")
    src.add_argument("--reportdir", help="Directory to glob for *.BUSCO_summary.*.txt")
    p.add_argument("-o", "--output", required=True, help="output TSV (.gz -> gzip-compressed)")
    args = p.parse_args()

    report_files = collect_report_files(args)
    if not report_files:
        src_desc = args.manifest or args.reportdir
        print(f"No BUSCO summary files found ({src_desc})", file=sys.stderr)
        sys.exit(1)

    rows = []
    unparsed = []
    for path in report_files:
        # Filename: {ASMID}.BUSCO_summary.{lineage}.txt -- ASMID may itself
        # contain dots, so split on the fixed ".BUSCO_summary." marker rather
        # than on the first ".".
        base = os.path.basename(path)
        marker = ".BUSCO_summary."
        if marker not in base:
            unparsed.append(base)
            continue
        asmid, rest = base.split(marker, 1)
        lineage = rest[:-len(".txt")] if rest.endswith(".txt") else rest

        stats = parse_summary(path)
        if stats is None:
            unparsed.append(base)
            continue

        row = {"ASMID": asmid, "lineage": lineage}
        row.update(stats)
        rows.append(row)

    if unparsed:
        print(f"WARNING: {len(unparsed)} files could not be parsed:", file=sys.stderr)
        for u in unparsed[:10]:
            print(f"  {u}", file=sys.stderr)
        if len(unparsed) > 10:
            print(f"  ... and {len(unparsed) - 10} more", file=sys.stderr)

    with open_out(args.output) as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
