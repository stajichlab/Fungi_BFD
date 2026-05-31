#!/usr/bin/env python3
"""Parse results/asm_reports/*.stats.txt, join with samples.csv, write tables/asm_stats.tsv."""

import csv
import os
import sys
import glob

REPORT_DIR = "results/asm_reports"
SAMPLES_CSV = "samples.csv"
OUT_FILE = "tables/asm_stats.tsv"

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
                "SPECIES": row["SPECIES_IN"],
                "STRAIN": row["STRAIN"],
            }
    return samples


def main():
    """Join asm stats files with samples metadata and write a TSV summary to tables/asm_stats.tsv."""
    os.makedirs("tables", exist_ok=True)

    samples = load_samples(SAMPLES_CSV)

    report_files = sorted(glob.glob(os.path.join(REPORT_DIR, "*.stats.txt")))
    if not report_files:
        print(f"No stats files found in {REPORT_DIR}", file=sys.stderr)
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
        print(f"WARNING: {len(missing_meta)} ASMIDs not found in {SAMPLES_CSV}:", file=sys.stderr)
        for m in missing_meta[:10]:
            print(f"  {m}", file=sys.stderr)
        if len(missing_meta) > 10:
            print(f"  ... and {len(missing_meta) - 10} more", file=sys.stderr)

    with open(OUT_FILE, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {OUT_FILE}")


if __name__ == "__main__":
    main()
