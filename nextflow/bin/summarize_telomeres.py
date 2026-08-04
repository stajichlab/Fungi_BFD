#!/usr/bin/env python3
"""Aggregate per-genome telomere feature tables into a merged summary.

Inputs are supplied as a manifest file (--manifest, one TAB-delimited
``<path>\\t<mtime>\\t<size>`` record per line) or as a directory to glob
(--reportdir).  Output is written to -o; if it ends in .gz it is gzip-compressed.

Per-genome aggregates include total telomeric tract length, total repeat count,
number of scaffolds with at least one telomere, counts by scaffold end and
strand, and the most frequently observed monomer.
"""

import argparse
import csv
import glob
import gzip
import os
import sys
from collections import Counter


COLUMNS = [
    "ASMID", "SPECIES", "STRAIN",
    "telomere_scaffolds",
    "telomere_tracts",
    "telomere_total_length_bp",
    "telomere_total_repeats",
    "telomere_5prime_count",
    "telomere_3prime_count",
    "telomere_plus_count",
    "telomere_minus_count",
    "telomere_terminal_count",
    "telomere_internal_count",
    "telomere_top_monomer",
    "telomere_monomers",
]


def parse_int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return 0


def parse_bool(value):
    return value.strip().lower() in ("true", "1", "yes", "t")


def load_samples(path):
    """Return a dict mapping ASMID -> {SPECIES, STRAIN}."""
    samples = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = row["ASMID"]
            samples[asmid] = {
                "SPECIES": row.get("SPECIES_IN", row.get("SPECIES", "")),
                "STRAIN": row.get("STRAIN", ""),
            }
    return samples


def read_manifest(path):
    """Read a toManifest file; return the list of report paths."""
    files = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            files.append(line.split("\t", 1)[0])
    return files


def open_in(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "r", newline="")


def open_out(path):
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def aggregate_genome(path):
    """Aggregate a single genome's telomere TSV into a dict."""
    scaffolds = set()
    total_length = 0
    total_repeats = 0
    end_5 = 0
    end_3 = 0
    strand_plus = 0
    strand_minus = 0
    terminal = 0
    internal = 0
    monomer_counts = Counter()

    with open_in(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            scaffolds.add(row["scaffold"])
            total_length += parse_int(row["tract_length"])
            total_repeats += parse_int(row["repeat_count"])
            if row["end"] == "5prime":
                end_5 += 1
            elif row["end"] == "3prime":
                end_3 += 1
            if row["strand"] == "+":
                strand_plus += 1
            elif row["strand"] == "-":
                strand_minus += 1
            if parse_bool(row.get("terminal", "True")):
                terminal += 1
            else:
                internal += 1
            monomer_counts[row["monomer"]] += 1

    top_monomer = monomer_counts.most_common(1)[0][0] if monomer_counts else ""
    return {
        "telomere_scaffolds": len(scaffolds),
        "telomere_tracts": sum(monomer_counts.values()),
        "telomere_total_length_bp": total_length,
        "telomere_total_repeats": total_repeats,
        "telomere_5prime_count": end_5,
        "telomere_3prime_count": end_3,
        "telomere_plus_count": strand_plus,
        "telomere_minus_count": strand_minus,
        "telomere_terminal_count": terminal,
        "telomere_internal_count": internal,
        "telomere_top_monomer": top_monomer,
        "telomere_monomers": ";".join(sorted(monomer_counts)),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--manifest", help="TAB-delimited manifest of report paths")
    src.add_argument("--reportdir", help="Directory to glob for *.telomeres.tsv.gz")
    p.add_argument("--samples", required=True, help="samples CSV with ASMID/SPECIES_IN/STRAIN")
    p.add_argument("-o", "--output", required=True, help="output TSV (.gz -> gzip-compressed)")
    args = p.parse_args()

    samples = load_samples(args.samples)

    if args.manifest:
        report_files = read_manifest(args.manifest)
    else:
        report_files = glob.glob(os.path.join(args.reportdir, "**/*.telomeres.tsv.gz"), recursive=True)
    report_files = sorted(report_files)

    if not report_files:
        src_desc = args.manifest or args.reportdir
        print(f"No telomere report files found ({src_desc})", file=sys.stderr)
        sys.exit(1)

    rows = []
    missing_meta = []
    for path in report_files:
        stem = os.path.basename(path).replace(".telomeres.tsv.gz", "").replace(".telomeres.tsv", "")
        agg = aggregate_genome(path)

        meta = samples.get(stem)
        if meta is None:
            missing_meta.append(stem)
            species = ""
            strain = ""
        else:
            species = meta["SPECIES"]
            strain = meta["STRAIN"]

        row = {"ASMID": stem, "SPECIES": species, "STRAIN": strain}
        row.update(agg)
        rows.append(row)

    if missing_meta:
        print(f"WARNING: {len(missing_meta)} ASMIDs not found in {args.samples}:", file=sys.stderr)
        for m in missing_meta[:10]:
            print(f"  {m}", file=sys.stderr)
        if len(missing_meta) > 10:
            print(f"  ... and {len(missing_meta) - 10} more", file=sys.stderr)

    with open_out(args.output) as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()
