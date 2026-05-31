#!/usr/bin/env python3
"""Fix FASTQ read headers: strip description and append /1 or /2."""

import argparse
import gzip
import re
import sys


def open_fastq(path):
    """Open a plain or gzipped FASTQ file for reading, or return stdin if path is None."""
    if path is None:
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def fix_headers(infile, read_num):
    """Strip description from each FASTQ header in infile and append /{read_num}, writing to stdout."""
    suffix = f"/{read_num}"
    fh = open_fastq(infile)
    try:
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            # Replace header: keep only @identifier, drop description
            header = re.sub(r'^(@\S+).+', r'\1' + suffix, header.rstrip()) + '\n'
            sys.stdout.write(header + seq + "+\n" + qual)
    finally:
        if fh is not sys.stdin:
            fh.close()


def main():
    """Parse arguments and run fix_headers to rewrite FASTQ headers with /1 or /2 suffixes."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fastq", nargs="?", default=None, help="Input FASTQ file (plain or gzipped); omit to read from stdin")
    parser.add_argument(
        "read_num",
        choices=["1", "2"],
        help="Read number to append (/1 or /2)",
    )
    args = parser.parse_args()
    fix_headers(args.fastq, args.read_num)


if __name__ == "__main__":
    main()
