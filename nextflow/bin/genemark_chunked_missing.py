#!/usr/bin/env python3
"""List genomes that took the chunked GeneMark route but have no result in genemark_store.

genemark_chunked.tsv (written by GENEMARK_TRAIN, one row per chunk) records every genome that
entered the scatter/gather route. With errorStrategy 'ignore' a permanently failed chunk or merge
leaves no GTF, and nothing else records that, so compare against the store after each run.

A genome counts as missing when <store>/<out>/<out>.genemark.gtf does not exist or is empty.
Usage: genemark_chunked_missing.py --chunked-tsv <target>/genemark_chunked.tsv --store <genemark_store> > missing.tsv
Reads plain or .gz input.
"""
import argparse
import csv
import gzip
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--chunked-tsv", required=True)
    ap.add_argument("--store", required=True)
    a = ap.parse_args()
    opener = gzip.open if a.chunked_tsv.endswith(".gz") else open
    seen = {}
    with opener(a.chunked_tsv, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            seen[r["out"]] = r  # last row per genome wins (a rerun appends new rows)
    w = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    w.writerow(["out", "asmid", "route", "n_chunks", "genome_contigs", "store_status"])
    missing = 0
    for out, r in sorted(seen.items()):
        gtf = os.path.join(a.store, out, out + ".genemark.gtf")
        if not os.path.exists(gtf):
            status = "missing"
        elif os.path.getsize(gtf) == 0:
            status = "empty"
        else:
            continue
        w.writerow([out, r["asmid"], r["route"], r["n_chunks"], r.get("genome_contigs") or r.get("contigs", ""), status])
        missing += 1
    print(f"{len(seen)} chunked genomes, {missing} without a GeneMark result in {a.store}", file=sys.stderr)


main()
