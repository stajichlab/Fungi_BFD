#!/usr/bin/env python3
"""Write samples.csv (optionally subset to matched/processed rows) for the run manifest.

The BFD.nf MERGE_SAMPLES step passes the full samples.csv here. Per T-014 §D.2,
MERGE_SAMPLES always builds the full, unscoped table now -- --taxon runs still
only *process* that clade (via taxonRowFilter() at the genome channel), they no
longer produce a separate taxon-restricted samples.csv.gz. --matched is kept as
an optional filter (mirrors build_species_table.py) for any future caller that
still wants a restricted subset; omit it to keep every row.
"""

import argparse
import csv
import gzip
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Write samples.csv, optionally subset to rows whose key was matched.")
    parser.add_argument("--samples", required=True,
                        help="full samples CSV (with header)")
    parser.add_argument("--matched",
                        help="file with one matched key per line; if omitted, keep all rows")
    parser.add_argument("--key", default="LOCUSTAG",
                        help="samples column to match against --matched values [LOCUSTAG]")
    parser.add_argument("-o", "--outfile", default="tables/samples.csv.gz",
                        help="output CSV; gzipped when name ends in .gz [tables/samples.csv.gz]")
    args = parser.parse_args()

    keep = None
    if args.matched:
        with open(args.matched) as fh:
            keep = {line.strip() for line in fh if line.strip()}

    outdir = os.path.dirname(args.outfile)
    if outdir:
        os.makedirs(outdir, exist_ok=True)

    opener = gzip.open if args.outfile.endswith(".gz") else open
    n = 0
    with open(args.samples, newline="") as ifh, opener(args.outfile, "wt", newline="") as ofh:
        reader = csv.DictReader(ifh)
        if args.key not in (reader.fieldnames or []):
            print(f"ERROR: key column '{args.key}' not in {args.samples} "
                  f"(have: {reader.fieldnames})", file=sys.stderr)
            sys.exit(1)
        writer = csv.DictWriter(ofh, fieldnames=reader.fieldnames)
        writer.writeheader()
        for row in reader:
            if keep is not None and (row.get(args.key) or "").strip() not in keep:
                continue
            writer.writerow(row)
            n += 1

    if keep is not None:
        print(f"Wrote {n} of {len(keep)} matched rows to {args.outfile}")
    else:
        print(f"Wrote {n} rows to {args.outfile}")


if __name__ == "__main__":
    main()
