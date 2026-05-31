#!/usr/bin/env python3
"""Report scaffold count, assembly size, and softmasked (lowercase) repeat content of a FASTA file."""

import sys
import argparse
from pathlib import Path


def count_mask_stats(fasta_path):
    """Return (num_scaffolds, total_bases, masked_bases) for a softmasked FASTA file."""
    num_scaffolds = 0
    total_bases = 0
    masked_bases = 0

    with open(fasta_path) as fh:
        seq_chars = []
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq_chars:
                    seq = "".join(seq_chars)
                    total_bases += len(seq)
                    masked_bases += sum(1 for c in seq if c.islower())
                    seq_chars = []
                num_scaffolds += 1
            else:
                seq_chars.append(line)
        if seq_chars:
            seq = "".join(seq_chars)
            total_bases += len(seq)
            masked_bases += sum(1 for c in seq if c.islower())

    return num_scaffolds, total_bases, masked_bases


def fmt(n):
    """Format an integer with thousands separators."""
    return f"{n:,}"


def main():
    """Print scaffold count, assembly size, and softmasked repeat percentage for a FASTA file."""
    parser = argparse.ArgumentParser(
        description="Report softmask (lowercase) statistics for a FASTA file."
    )
    parser.add_argument("fasta", help="Input FASTA file (softmasked)")
    args = parser.parse_args()

    path = Path(args.fasta)
    if not path.exists():
        sys.exit(f"Error: file not found: {args.fasta}")

    num_scaffolds, total_bases, masked_bases = count_mask_stats(path)

    pct = (masked_bases / total_bases * 100) if total_bases else 0.0

    print(f"num scaffolds: {fmt(num_scaffolds)}")
    print(f"assembly size: {fmt(total_bases)} bp")
    print(f"masked repeats: {fmt(masked_bases)} bp ({pct:.2f}%)")


if __name__ == "__main__":
    main()
