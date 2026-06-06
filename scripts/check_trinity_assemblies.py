#!/usr/bin/env python3
"""
Identify failed/low-quality Trinity genome-guided assemblies and summarize stats
for all successful ones.

Outputs:
  TrinityGG_failed.tsv   — species with missing, empty, or low-transcript assemblies
  TrinityGG_summary.tsv  — per-species assembly stats for all non-empty assemblies
"""

import os
import re
import sys
import math
import argparse
from pathlib import Path
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reads-dir",   default="rnaseq_reads",
                   help="Directory containing paired-end FASTQ files")
    p.add_argument("--asm-dir",     default="rnaseq_data",
                   help="Directory containing Trinity GG FASTA assemblies")
    p.add_argument("--failed-out",  default="TrinityGG_failed.tsv")
    p.add_argument("--summary-out", default="TrinityGG_summary.tsv")
    p.add_argument("--min-transcripts", type=int, default=2000,
                   help="Minimum transcript count; below this is flagged as LOW_COUNT")
    p.add_argument("--write-all", action="store_true",
                   help="Write all failed entries including those with empty R1 reads "
                        "(default: skip species whose R1 file is 0 bytes)")
    p.add_argument("--show-path", action="store_true",
                   help="Include full paths in R1, R2, and FASTA columns "
                        "(default: basename only)")
    return p.parse_args()


# ── length extraction ──────────────────────────────────────────────────────────

LEN_RE = re.compile(r"\blen=(\d+)\b")

def transcript_lengths_from_fasta(fasta_path):
    """
    Extract transcript lengths from a Trinity FASTA.
    Trinity embeds 'len=NNN' in every header, so we never read sequence lines.
    Falls back to counting sequence characters when the tag is absent.
    Returns a list of ints (one per transcript).
    """
    lengths = []
    current_len = None
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                if current_len is not None:
                    lengths.append(current_len)
                    current_len = None
                m = LEN_RE.search(line)
                if m:
                    lengths.append(int(m.group(1)))
                else:
                    current_len = 0   # will accumulate sequence
            else:
                if current_len is not None:
                    current_len += len(line)
        if current_len is not None:
            lengths.append(current_len)
    return lengths


# ── assembly statistics ────────────────────────────────────────────────────────

def n50(sorted_lengths_desc, total):
    """Compute N50 from lengths already sorted descending."""
    half = total / 2
    running = 0
    for l in sorted_lengths_desc:
        running += l
        if running >= half:
            return l
    return 0


def assembly_stats(lengths):
    """Return a dict of summary statistics for a list of transcript lengths."""
    if not lengths:
        return None
    lengths_sorted = sorted(lengths, reverse=True)
    n       = len(lengths_sorted)
    total   = sum(lengths_sorted)
    mean    = total / n
    median  = (lengths_sorted[n // 2] if n % 2 == 1
               else (lengths_sorted[n // 2 - 1] + lengths_sorted[n // 2]) / 2)
    n50_val = n50(lengths_sorted, total)

    # GC / composition not available without sequences; skip.
    # Additional useful stats:
    #   N90  — length at which 90 % of bases are covered
    #   ExN50 (expression-weighted) not computed here (needs abundance)
    n90_val = _nx(lengths_sorted, total, 0.90)
    longest  = lengths_sorted[0]
    shortest = lengths_sorted[-1]

    # Length bins
    bins = {">= 1 kb": 0, ">= 2 kb": 0, ">= 5 kb": 0}
    for l in lengths:
        if l >= 1000: bins[">= 1 kb"] += 1
        if l >= 2000: bins[">= 2 kb"] += 1
        if l >= 5000: bins[">= 5 kb"] += 1

    return {
        "NUM_TRANSCRIPTS":   n,
        "TOTAL_LENGTH":      total,
        "MEAN_LENGTH":       round(mean, 2),
        "MEDIAN_LENGTH":     median,
        "N50":               n50_val,
        "N90":               n90_val,
        "MIN_LENGTH":        shortest,
        "MAX_LENGTH":        longest,
        "TRANSCRIPTS_GE1KB": bins[">= 1 kb"],
        "TRANSCRIPTS_GE2KB": bins[">= 2 kb"],
        "TRANSCRIPTS_GE5KB": bins[">= 5 kb"],
    }


def _nx(sorted_desc, total, fraction):
    target = total * fraction
    running = 0
    for l in sorted_desc:
        running += l
        if running >= target:
            return l
    return 0


# ── reads discovery ───────────────────────────────────────────────────────────

READS_SUFFIXES = [
    ("_norm_R1.fastq.gz",    "_norm_R2.fastq.gz"),
    ("_trimmed_R1.fastq.gz", "_trimmed_R2.fastq.gz"),
    ("_R1.fastq.gz",         "_R2.fastq.gz"),
]

def discover_species_reads(reads_dir):
    """
    Return dict: species_name -> (r1_path, r2_path)
    When a species has both norm and trimmed reads, prefer norm.
    """
    files = set(os.listdir(reads_dir))
    species = {}
    for fname in sorted(files):
        for s1, s2 in READS_SUFFIXES:
            if fname.endswith(s1):
                name = fname[: -len(s1)]
                r2   = name + s2
                if r2 in files:
                    if name not in species:
                        species[name] = (
                            os.path.join(reads_dir, fname),
                            os.path.join(reads_dir, r2),
                        )
                    break  # stop trying other suffix patterns for this file
    return species


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    args = parse_args()

    reads_dir = Path(args.reads_dir)
    asm_dir   = Path(args.asm_dir)

    print(f"Scanning reads in:     {reads_dir}", file=sys.stderr)
    print(f"Scanning assemblies in: {asm_dir}", file=sys.stderr)

    species_reads = discover_species_reads(reads_dir)
    print(f"Found {len(species_reads)} species with read pairs", file=sys.stderr)

    failed_rows  = []
    summary_rows = []

    skipped_empty_reads = 0

    for i, (sp, (r1, r2)) in enumerate(sorted(species_reads.items()), 1):
        if i % 500 == 0:
            print(f"  ... processed {i}/{len(species_reads)}", file=sys.stderr)

        # ── skip species with empty R1 unless --write-all ──────────────────
        if not args.write_all and Path(r1).stat().st_size == 0:
            skipped_empty_reads += 1
            continue

        fasta = asm_dir / f"{sp}.trinity-GG.fasta"

        # ── assembly missing entirely ──────────────────────────────────────
        if not fasta.exists():
            failed_rows.append({
                "SPECIES":          sp,
                "R1":               r1,
                "R2":               r2,
                "FASTA":            str(fasta),
                "REASON":           "MISSING",
                "NUM_TRANSCRIPTS":  0,
            })
            continue

        # ── assembly empty ─────────────────────────────────────────────────
        if fasta.stat().st_size == 0:
            failed_rows.append({
                "SPECIES":          sp,
                "R1":               r1,
                "R2":               r2,
                "FASTA":            str(fasta),
                "REASON":           "EMPTY",
                "NUM_TRANSCRIPTS":  0,
            })
            continue

        # ── parse transcript lengths ───────────────────────────────────────
        lengths = transcript_lengths_from_fasta(fasta)
        n_tx    = len(lengths)

        if n_tx < args.min_transcripts:
            failed_rows.append({
                "SPECIES":          sp,
                "R1":               r1,
                "R2":               r2,
                "FASTA":            str(fasta),
                "REASON":           "LOW_COUNT",
                "NUM_TRANSCRIPTS":  n_tx,
            })
        # Include in summary regardless of count (skip empty only)
        stats = assembly_stats(lengths)
        if stats:
            stats["SPECIES"] = sp
            stats["FASTA"]   = str(fasta)
            summary_rows.append(stats)

    if skipped_empty_reads:
        print(f"Skipped {skipped_empty_reads} species with empty R1 files "
              f"(use --write-all to include them)", file=sys.stderr)

    # ── path formatting helper ─────────────────────────────────────────────────
    def fmt(path_str):
        return path_str if args.show_path else os.path.basename(path_str)

    # ── write failed TSV ───────────────────────────────────────────────────────
    failed_cols = ["SPECIES", "R1", "R2", "FASTA", "REASON", "NUM_TRANSCRIPTS"]
    with open(args.failed_out, "w") as fh:
        fh.write("\t".join(failed_cols) + "\n")
        for row in failed_rows:
            fh.write("\t".join([
                row["SPECIES"],
                fmt(row["R1"]),
                fmt(row["R2"]),
                fmt(row["FASTA"]),
                row["REASON"],
                str(row["NUM_TRANSCRIPTS"]),
            ]) + "\n")

    print(f"\nWrote {len(failed_rows)} failed entries to {args.failed_out}", file=sys.stderr)

    # ── write summary TSV ──────────────────────────────────────────────────────
    summary_cols = [
        "SPECIES", "FASTA",
        "NUM_TRANSCRIPTS", "TOTAL_LENGTH", "MEAN_LENGTH", "MEDIAN_LENGTH",
        "N50", "N90",
        "MIN_LENGTH", "MAX_LENGTH",
        "TRANSCRIPTS_GE1KB", "TRANSCRIPTS_GE2KB", "TRANSCRIPTS_GE5KB",
    ]
    with open(args.summary_out, "w") as fh:
        fh.write("\t".join(summary_cols) + "\n")
        for row in sorted(summary_rows, key=lambda r: r["SPECIES"]):
            fh.write("\t".join(
                fmt(str(row[c])) if c == "FASTA" else str(row[c])
                for c in summary_cols
            ) + "\n")

    print(f"Wrote {len(summary_rows)} assembly summaries to {args.summary_out}", file=sys.stderr)

    # ── breakdown ──────────────────────────────────────────────────────────────
    reasons = defaultdict(int)
    for r in failed_rows:
        reasons[r["REASON"]] += 1
    print("\nFailure breakdown:", file=sys.stderr)
    for reason, count in sorted(reasons.items()):
        print(f"  {reason:12s} {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
