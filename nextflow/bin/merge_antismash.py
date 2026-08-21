#!/usr/bin/env python3
"""Merge per-genome antiSMASH region JSON into a tables-loadable CSV.

antiSMASH's own *.json output nests BGC calls under records[].areas[]
(recent versions) or records[].regions[] (older versions) -- both are
handled here since the exact key has shifted across antiSMASH releases.
"""

import csv
import gzip
import json
import sys
import argparse


def iter_areas(record):
    for key in ("areas", "regions"):
        for area in record.get(key, []):
            yield area


def main():
    parser = argparse.ArgumentParser(description="Merge antiSMASH region JSON results")
    parser.add_argument("jsons", nargs="+", help="*.antismash.json.gz files")
    parser.add_argument("-o", "--outfile", default="antismash.csv")
    args = parser.parse_args()

    with open(args.outfile, "w", newline="") as of:
        w = csv.writer(of)
        w.writerow(['species_prefix', 'record_id', 'cluster_num',
                    'product', 'contig_edge', 'start', 'end'])
        for f in sorted(args.jsons):
            species_prefix = f.split("/")[-1].split(".")[0]
            with gzip.open(f, "rt") as fh:
                try:
                    data = json.load(fh)
                except json.JSONDecodeError:
                    continue
            for record in data.get("records", []):
                record_id = record.get("id", "")
                for i, area in enumerate(iter_areas(record), start=1):
                    products = area.get("products") or [area.get("product", "")]
                    w.writerow([
                        species_prefix, record_id, i,
                        ";".join(p for p in products if p),
                        area.get("contig_edge", ""),
                        area.get("start", ""), area.get("end", ""),
                    ])

    print(f"Written: {args.outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()
