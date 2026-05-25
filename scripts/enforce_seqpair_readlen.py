#!/usr/bin/env python3
"""
Enforce equal read lengths in paired-end FASTQ files.

For each read pair (matched by name, R1 has /1 suffix, R2 has /2 suffix):
  - If either read is below --min-len, discard the pair entirely.
  - If the reads differ in length, truncate the longer to match the shorter.

Usage (BBTools-style key=value parameters):
  enforce_pair_readlen_for_normalization.py in=R1.fastq.gz in2=R2.fastq.gz \
      out=R1_trunc.fastq.gz out2=R2_trunc.fastq.gz minlen=75
"""

import argparse
import gzip
import sys
import os
from itertools import zip_longest


def open_fastq(path, mode="r"):
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")
    return open(path, mode)


def fastq_records(fh):
    """Yield (name, seq, plus, qual) tuples from an open FASTQ file handle."""
    while True:
        name = fh.readline().rstrip("\n")
        if not name:
            break
        seq  = fh.readline().rstrip("\n")
        plus = fh.readline().rstrip("\n")
        qual = fh.readline().rstrip("\n")
        yield name, seq, plus, qual


def base_name(read_name):
    """Strip /1 or /2 (and leading @) to get the canonical pair key."""
    n = read_name.lstrip("@")
    if n.endswith("/1") or n.endswith("/2"):
        n = n[:-2]
    return n


def write_record(fh, name, seq, plus, qual):
    fh.write(f"{name}\n{seq}\n{plus}\n{qual}\n")


def parse_kv_args(argv):
    """
    Accept both --flag style and key=value style arguments.
    Returns a namespace with: in1, in2, out1, out2, minlen, stats.
    """
    # Pre-process argv to convert key=value to --key value
    converted = []
    for token in argv:
        if "=" in token and not token.startswith("-"):
            k, v = token.split("=", 1)
            # map BBTools aliases
            if k == "in":
                if any(a.startswith("in=") or a == "--in" for a in converted):
                    k = "in2"  # second 'in=' becomes in2
                else:
                    k = "in"
            if k == "out":
                if any(a.startswith("--out") and not a.startswith("--out2") for a in converted):
                    k = "out2"
                else:
                    k = "out"
            converted.extend([f"--{k}", v])
        else:
            converted.append(token)

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--in",  dest="in1",  required=True,  help="Forward (R1) FASTQ input")
    parser.add_argument("--in2", dest="in2",  required=True,  help="Reverse (R2) FASTQ input")
    parser.add_argument("--out", dest="out1", default=None,   help="Forward (R1) FASTQ output")
    parser.add_argument("--out2",dest="out2", default=None,   help="Reverse (R2) FASTQ output")
    parser.add_argument("--minlen", dest="minlen", type=int, default=75,
                        help="Minimum read length; discard pair if either read is shorter (default: 75)")
    parser.add_argument("--stats", dest="stats", default=None,
                        help="Write summary statistics to this file (default: stderr)")
    return parser.parse_args(converted)


def default_output(path, suffix="_trunc"):
    """Derive a default output path by inserting a suffix before the extension."""
    base = path
    ext = ""
    for gz in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if base.endswith(gz):
            base = base[: -len(gz)]
            ext = gz
            break
    return base + suffix + ext


def main():
    args = parse_kv_args(sys.argv[1:])

    if args.out1 is None:
        args.out1 = default_output(args.in1)
    if args.out2 is None:
        args.out2 = default_output(args.in2)

    n_pairs = 0
    n_kept = 0
    n_too_short = 0
    n_truncated = 0

    with open_fastq(args.in1) as fh1, \
         open_fastq(args.in2) as fh2, \
         open_fastq(args.out1, "w") as oh1, \
         open_fastq(args.out2, "w") as oh2:

        for (r1, r2) in zip_longest(fastq_records(fh1), fastq_records(fh2)):
            if r1 is None or r2 is None:
                sys.exit(
                    "ERROR: R1 and R2 files have different numbers of records. "
                    "Ensure files are properly paired."
                )

            n1, s1, p1, q1 = r1
            n2, s2, p2, q2 = r2

            if base_name(n1) != base_name(n2):
                sys.exit(
                    f"ERROR: Read name mismatch at pair {n_pairs + 1}:\n"
                    f"  R1: {n1}\n  R2: {n2}\n"
                    "Ensure files are sorted in the same order."
                )

            n_pairs += 1
            len1, len2 = len(s1), len(s2)
            min_len = min(len1, len2)

            # Discard if either read is below minimum length
            if min_len < args.minlen:
                n_too_short += 1
                continue

            # Truncate to the shorter length if they differ
            if len1 != len2:
                n_truncated += 1
                s1 = s1[:min_len]
                q1 = q1[:min_len]
                s2 = s2[:min_len]
                q2 = q2[:min_len]

            write_record(oh1, n1, s1, p1, q1)
            write_record(oh2, n2, s2, p2, q2)
            n_kept += 1

    summary = (
        f"Pairs processed : {n_pairs}\n"
        f"Pairs kept      : {n_kept}\n"
        f"Pairs discarded (too short, min={args.minlen}): {n_too_short}\n"
        f"Pairs truncated : {n_truncated}\n"
        f"Output R1       : {args.out1}\n"
        f"Output R2       : {args.out2}\n"
    )

    if args.stats:
        with open(args.stats, "w") as sf:
            sf.write(summary)
    else:
        sys.stderr.write(summary)


if __name__ == "__main__":
    main()
