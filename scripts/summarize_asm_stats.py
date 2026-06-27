#!/usr/bin/env python3
"""Parse AAFTF `assess` *.stats.txt reports, join with samples.csv, write a merged TSV.

Inputs are supplied either as a manifest file (--manifest, one TAB-delimited
`<path>\\t<mtime>\\t<size>` record per line, as written by BFD.nf toManifest) or
as a directory to glob (--reportdir).  Output is written to -o; if it ends in
.gz it is gzip-compressed.
"""

import argparse
import csv
import glob
import gzip
import os
import sys

FIELD_MAP = {
    "CONTIG COUNT": "contig_count",
    "TOTAL LENGTH": "total_length_bp",
    "MIN": "min_contig_bp",
    "MAX": "max_contig_bp",
    "MEDIAN": "median_contig_bp",
    "MEAN": "mean_contig_bp",
    "L50": "L50",
    "N50": "N50_bp",
    "L90": "L90",
    "N90": "N90_bp",
    "GC%": "gc_pct",
    "N GAP COUNT": "n_gap_count",
    "TOTAL N BASES": "total_n_bases",
    "BASES MASKED": "masked_bases",
    "PERCENT MASKED": "masked_pct",
    "T2T SCAFFOLDS": "t2t_scaffolds",
    "TELOMERE FWD": "telomere_fwd",
    "TELOMERE REV": "telomere_rev",
}

COLUMNS = [
    "ASMID", "SPECIES", "STRAIN",
    "contig_count", "total_length_bp",
    "min_contig_bp", "max_contig_bp", "median_contig_bp", "mean_contig_bp",
    "L50", "N50_bp", "L90", "N90_bp",
    "gc_pct", "n_gap_count", "total_n_bases",
    "masked_bases", "masked_pct",
    "t2t_scaffolds", "telomere_fwd", "telomere_rev",
]


def parse_stats(path):
    """Parse a .stats.txt file and return a dict of normalized field names to values."""
    stats = {}
    with open(path) as fh:
        for line in fh:
            if "=" not in line:
                continue
            key, _, val = line.partition("=")
            key = key.strip()
            val = val.strip()
            if key in FIELD_MAP:
                stats[FIELD_MAP[key]] = val
    return stats


def load_samples(path):
    """Return a dict mapping ASMID → {SPECIES, STRAIN} from the samples CSV at path."""
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
    """Resolve the list of .stats.txt report files from --manifest or --reportdir."""
    if args.manifest:
        files = read_manifest(args.manifest)
    else:
        files = glob.glob(os.path.join(args.reportdir, "*.stats.txt"))
    return sorted(files)


def open_out(path):
    """Open path for text writing, gzip-compressed if it ends in .gz."""
    if path.endswith(".gz"):
        return gzip.open(path, "wt", newline="")
    return open(path, "w", newline="")


def main():
    """Join asm stats reports with samples metadata and write a merged TSV."""
    p = argparse.ArgumentParser(description=__doc__)
    src = p.add_mutually_exclusive_group(required=True)
    src.add_argument("--manifest", help="TAB-delimited manifest of report paths (path\\tmtime\\tsize)")
    src.add_argument("--reportdir", help="Directory to glob for *.stats.txt")
    p.add_argument("--samples", required=True, help="samples CSV with ASMID/SPECIES_IN/STRAIN")
    p.add_argument("-o", "--output", required=True, help="output TSV (.gz → gzip-compressed)")
    args = p.parse_args()

    samples = load_samples(args.samples)

    report_files = collect_report_files(args)
    if not report_files:
        src_desc = args.manifest or args.reportdir
        print(f"No stats files found ({src_desc})", file=sys.stderr)
        sys.exit(1)

    rows = []
    missing_meta = []
    for path in report_files:
        stem = os.path.basename(path).replace(".stats.txt", "")
        stats = parse_stats(path)

        meta = samples.get(stem)
        if meta is None:
            missing_meta.append(stem)
            species = ""
            strain = ""
        else:
            species = meta["SPECIES"]
            strain = meta["STRAIN"]

        row = {"ASMID": stem, "SPECIES": species, "STRAIN": strain}
        row.update(stats)
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
