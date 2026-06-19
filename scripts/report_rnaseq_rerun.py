#!/usr/bin/env python3
"""Report species whose RNA-seq reads need to be re-fetched.

A species "needs a re-run" when its SRA query found usable accessions (the
per-species sra_query CSV has at least one data row beyond the header) but the
normalized read files in the reads directory do not actually carry data.

SRA_FETCH / SRA_FETCH_SE / WRITE_EMPTY_READS always emit a 4-tuple of files per
species, but only the relevant ones are populated:
  - paired-end fetch  -> <tag>_norm_R1.fastq.gz + <tag>_norm_R2.fastq.gz hold data,
                         <tag>_norm_SE.fastq.gz is a zero-byte stub
  - single-end fetch  -> <tag>_norm_SE.fastq.gz holds data, R1/R2 are zero-byte stubs
  - no SRA data       -> all three are zero-byte stubs (this is correct, not a failure)

So a species is considered OK when (R1 and R2 both have data) OR (SE has data).
Anything else, for a species that *did* have accessions queued, is flagged:
  - all_empty   : accessions were found but every read file is empty/missing
  - broken_pair : exactly one of R1/R2 has data and SE is empty (mismatched pair)
  - no_files    : accessions were found but no *_norm_* files exist at all

A real fastq.gz with reads is many KB; empty stubs are 0 bytes (an empty gzip is
~20-30 bytes). Files at or below --min-bytes (default 50) are treated as empty.

Inputs:
  --query-dir  directory of <tag>.sra_query.csv files (default rnaseq_reads/sra_query)
  --reads-dir  directory of <tag>_norm_R1/R2/SE.fastq.gz files (default rnaseq_reads)

Output:
  A TSV report (stdout, or --out) of only the species needing a re-run, with the
  accession counts, file sizes, and the reason. A one-line summary goes to stderr.

Usage:
  python scripts/report_rnaseq_rerun.py
  python scripts/report_rnaseq_rerun.py --query-dir rnaseq_reads/sra_query \\
      --reads-dir rnaseq_reads --out rnaseq_rerun_report.tsv
"""
import argparse
import csv
import os
import sys
from glob import glob

QUERY_SUFFIX = ".sra_query.csv"


def size_or_missing(path):
    """Return file size in bytes, or -1 if the file does not exist."""
    try:
        return os.path.getsize(path)
    except OSError:
        return -1


def load_skip_accessions(blacklist_path):
    """Return the set of accessions marked action=skip in rnaseq_blacklist.csv.

    Blacklist column order is sra_accession,species_tag,taxonid,action. SE_trinity
    entries are NOT skips: they are re-routed to single-end fetch, so their reads
    should still appear (in the SE file) and they remain countable.
    """
    skip = set()
    if not blacklist_path or not os.path.exists(blacklist_path):
        return skip
    with open(blacklist_path, newline="") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split(",")
            if len(fields) >= 4 and fields[3].strip() == "skip":
                skip.add(fields[0].strip())
    return skip


def count_accessions(csv_path, skip_accessions):
    """Return (n_total, n_paired, n_single) usable data rows in a sra_query CSV.

    Accessions listed as skip in the blacklist are excluded entirely (they are
    never downloaded, so they don't make a species a re-run candidate).
    Handles both the 5-column legacy header (no layout) and the 6-column current
    header (trailing 'layout' = PAIRED|SINGLE). Rows whose layout is unknown are
    counted as paired, matching the SRA default the pipeline assumes.
    """
    n_total = n_paired = n_single = 0
    with open(csv_path, newline="") as fh:
        reader = csv.reader(fh)
        for i, row in enumerate(reader):
            if i == 0 or not row or not row[0].strip():
                continue  # header or blank line
            acc = row[2].strip() if len(row) >= 3 else ""
            if acc in skip_accessions:
                continue  # blacklisted skip -> not expected to produce reads
            n_total += 1
            layout = row[5].strip().upper() if len(row) >= 6 else ""
            if layout == "SINGLE":
                n_single += 1
            else:
                n_paired += 1
    return n_total, n_paired, n_single


def classify(r1, r2, se, min_bytes):
    """Return (ok, reason) given the three read-file sizes (-1 == missing)."""
    has_r1 = r1 > min_bytes
    has_r2 = r2 > min_bytes
    has_se = se > min_bytes

    if (has_r1 and has_r2) or has_se:
        return True, "ok"
    if r1 < 0 and r2 < 0 and se < 0:
        return False, "no_files"
    if (has_r1 != has_r2) and not has_se:
        return False, "broken_pair"
    return False, "all_empty"


def fmt_size(n):
    return "NA" if n < 0 else str(n)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--query-dir", default="rnaseq_reads/sra_query",
                    help="directory of <tag>.sra_query.csv files")
    ap.add_argument("--reads-dir", default="rnaseq_reads",
                    help="directory of <tag>_norm_R1/R2/SE.fastq.gz files")
    ap.add_argument("--min-bytes", type=int, default=50,
                    help="files at or below this size are treated as empty (default 50)")
    ap.add_argument("--blacklist", default="rnaseq_blacklist.csv",
                    help="rnaseq_blacklist.csv; accessions marked skip are excluded "
                         "(set to '' to disable)")
    ap.add_argument("--out", default="-",
                    help="output TSV path (default stdout)")
    args = ap.parse_args()

    query_csvs = sorted(glob(os.path.join(args.query_dir, "*" + QUERY_SUFFIX)))
    if not query_csvs:
        sys.exit(f"[ERROR] no {QUERY_SUFFIX} files found in {args.query_dir}")

    skip_accessions = load_skip_accessions(args.blacklist)

    out = sys.stdout if args.out == "-" else open(args.out, "w", newline="")
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["species_tag", "n_accessions", "n_paired", "n_single",
                     "r1_bytes", "r2_bytes", "se_bytes", "status"])

    n_species = n_with_data = n_flagged = 0
    for csv_path in query_csvs:
        n_species += 1
        tag = os.path.basename(csv_path)[: -len(QUERY_SUFFIX)]
        n_total, n_paired, n_single = count_accessions(csv_path, skip_accessions)
        if n_total == 0:
            continue  # no accessions found -> empty reads are correct, not a re-run
        n_with_data += 1

        r1 = size_or_missing(os.path.join(args.reads_dir, f"{tag}_norm_R1.fastq.gz"))
        r2 = size_or_missing(os.path.join(args.reads_dir, f"{tag}_norm_R2.fastq.gz"))
        se = size_or_missing(os.path.join(args.reads_dir, f"{tag}_norm_SE.fastq.gz"))

        ok, reason = classify(r1, r2, se, args.min_bytes)
        if not ok:
            n_flagged += 1
            writer.writerow([tag, n_total, n_paired, n_single,
                             fmt_size(r1), fmt_size(r2), fmt_size(se), reason])

    if out is not sys.stdout:
        out.close()

    print(f"[INFO] {n_species} species queried; {n_with_data} have SRA accessions; "
          f"{n_flagged} need a re-run.", file=sys.stderr)


if __name__ == "__main__":
    main()
