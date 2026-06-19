#!/usr/bin/env python3
"""Write a samples CSV containing only the rows that were matched/processed in a run.

The BFD.nf MERGE_SAMPLES step passes the full samples.csv plus a file listing the
keys (one LOCUSTAG per line) of the genomes that survived the taxonomy/pattern
filters and were processed. The output keeps the original header and column order
and is gzipped, landing next to the other merged tables (e.g.
tables/<SUBSET>/samples.csv.gz). For a --taxon run this naturally contains only
the matched clade because the caller supplies only that clade's keys.
"""

import argparse
import csv
import gzip
import os
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Subset a samples CSV to the rows whose key was matched/processed.")
    parser.add_argument("--samples", required=True,
                        help="full samples CSV (with header)")
    parser.add_argument("--matched", required=True,
                        help="file with one matched key per line (default key column LOCUSTAG)")
    parser.add_argument("--key", default="LOCUSTAG",
                        help="samples column to match against --matched values [LOCUSTAG]")
    parser.add_argument("-o", "--outfile", default="tables/All_Taxa/samples.csv.gz",
                        help="output CSV; gzipped when name ends in .gz [tables/All_Taxa/samples.csv.gz]")
    args = parser.parse_args()

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
            if (row.get(args.key) or "").strip() in keep:
                writer.writerow(row)
                n += 1

    print(f"Wrote {n} of {len(keep)} matched rows to {args.outfile}")


if __name__ == "__main__":
    main()
