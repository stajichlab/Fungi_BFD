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
    asm_preflight_stats.py GENOME.fa[.gz] --report-repeat-pct   # adds a 5th column

Prints one TSV line to stdout: total_bp<TAB>contigs<TAB>n50<TAB>verdict
[<TAB>repeat_pct with --report-repeat-pct]
verdict is "small_fragmented" only when BOTH gates trip (small AND
fragmented) -- a complete small genome (e.g. Malassezia) is not flagged.
--min-bp 0 disables the guard entirely (verdict is always "ok").

repeat_pct (only with --report-repeat-pct, kept off by default so
GENEMARK_RUN's existing 4-variable `read` doesn't silently swallow a 5th
field into ASM_VERDICT): percentage of bases that are soft-masked
(lowercase acgtn) in the same single pass used for contig lengths, so
FUNANNOTATE_PREDICT can gate EVM's --repeats2evm / wider
--evm-partition-interval on assemblies this repeat-dense without a second
full-genome scan. See analysis of Austropuccinia_psidii (GCA_902702905.1,
65% masked, 368K raw ab-initio models) and GCA_003724095.1 EVM failures.
"""
import argparse
import gzip
import sys


def contig_stats(path):
    """Return (lengths, lowercase_base_count) in one pass over the FASTA."""
    opener = gzip.open if path.endswith(".gz") else open
    lengths = []
    length = 0
    lower = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if length:
                    lengths.append(length)
                length = 0
            else:
                seq = line.strip()
                length += len(seq)
                # str.count() is C-level; much faster than a per-char Python loop
                # over a multi-GB genome.
                lower += sum(seq.count(c) for c in "acgtn")
    if length:
        lengths.append(length)
    return lengths, lower


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
    ap.add_argument("--report-repeat-pct", action="store_true",
                     help="append a 5th column: pct of bases soft-masked (lowercase)")
    args = ap.parse_args()

    lengths, lower = contig_stats(args.genome)
    lengths.sort(reverse=True)
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

    line = f"{total_bp}\t{contigs}\t{n50}\t{verdict}"
    if args.report_repeat_pct:
        repeat_pct = (100.0 * lower / total_bp) if total_bp else 0.0
        line += f"\t{repeat_pct:.2f}"
    print(line)


if __name__ == "__main__":
    sys.exit(main())
