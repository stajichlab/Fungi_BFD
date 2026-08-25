#!/usr/bin/env python3
"""Cheap total-bp / contig-count / N50 stats for one FASTA, plus a
small-and-fragmented verdict.

Shared preflight guard for GENEMARK_RUN and FUNANNOTATE_PREDICT (both call
this instead of duplicating the same awk pipeline -- see
nextflow/docs/GENEMARK_RUN_DESIGN.md's "Known gap" section for why GENEMARK_RUN
needs the identical policy to FUNANNOTATE_PREDICT's own guard, upstream of it
in the DAG). Assemblies that are both small AND fragmented cannot yield
funannotate's required 30 training models, and starve GeneMark-ES/ET's own
training-contig selection (--min_contig 10000, after masking) down to nothing
usable -- both fail slowly (predict runs for hours; GeneMark burns a full
--ES/--ET attempt) instead of being skipped up front.

N50: sort contig lengths descending, walk until the cumulative sum reaches
half the assembly length, report that contig's length -- the standard
definition (matches seqkit stats / AAFTF assess).

Usage:
    asm_preflight_stats.py GENOME.fa[.gz] --min-bp N --max-n50 N --max-contigs N

Prints one TSV line to stdout: total_bp<TAB>contigs<TAB>n50<TAB>verdict
verdict is "small_fragmented" only when BOTH gates trip (small AND
fragmented) -- a complete small genome (e.g. Malassezia) is not flagged.
--min-bp 0 disables the guard entirely (verdict is always "ok").
"""
import argparse
import gzip
import sys


def contig_lengths(path):
    opener = gzip.open if path.endswith(".gz") else open
    lengths = []
    length = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if length:
                    lengths.append(length)
                length = 0
            else:
                length += len(line.strip())
    if length:
        lengths.append(length)
    return lengths


def n50_of(lengths_desc, total_bp):
    half = total_bp / 2
    running = 0
    for length in lengths_desc:
        running += length
        if running >= half:
            return length
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("genome", help="FASTA path, optionally gzip-compressed")
    ap.add_argument("--min-bp", type=int, default=0,
                     help="total assembled bp below this = 'small' (0 disables the whole guard)")
    ap.add_argument("--max-n50", type=int, default=0,
                     help="N50 below this = 'fragmented' (0 disables this gate)")
    ap.add_argument("--max-contigs", type=int, default=0,
                     help="contig count above this = 'fragmented' (0 disables this gate)")
    args = ap.parse_args()

    lengths = sorted(contig_lengths(args.genome), reverse=True)
    total_bp = sum(lengths)
    contigs = len(lengths)
    n50 = n50_of(lengths, total_bp) if lengths else 0

    verdict = "ok"
    if args.min_bp > 0:
        small = total_bp < args.min_bp
        fragmented = (args.max_n50 > 0 and n50 < args.max_n50) or (
            args.max_contigs > 0 and contigs > args.max_contigs
        )
        if small and fragmented:
            verdict = "small_fragmented"

    print(f"{total_bp}\t{contigs}\t{n50}\t{verdict}")


if __name__ == "__main__":
    sys.exit(main())
