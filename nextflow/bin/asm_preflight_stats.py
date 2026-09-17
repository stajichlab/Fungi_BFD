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
    asm_preflight_stats.py GENOME.fa[.gz] --min-contig-len N --min-training-contigs N
    asm_preflight_stats.py GENOME.fa[.gz] --report-repeat-pct   # adds a 5th column

Prints one TSV line to stdout: total_bp<TAB>contigs<TAB>n50<TAB>verdict
[<TAB>repeat_pct with --report-repeat-pct]
verdict is "small_fragmented" only when BOTH gates trip (small AND
fragmented) -- a complete small genome (e.g. Malassezia) is not flagged.
--min-bp 0 disables the guard entirely (verdict is always "ok").

verdict is "too_small" when total_bp is below --abs-min-bp, an UNCONDITIONAL
floor checked on its own regardless of fragmentation. Unlike --min-bp above,
this gate does not require the assembly to also be fragmented: no real
fungal genome assembly is a handful of large contigs totaling under this
floor, so a few large-ish contigs (i.e. "not fragmented" by the --max-n50/
--max-contigs gates) do not exempt it. Added after GCA_026108285.1
(Diplocarpon coronariae GW01-2013), whose own deposited assembly is
genuinely 3 contigs / ~106kb total -- small but not fragmented (N50 77.9kb),
so it cleared --min-bp's small_fragmented gate as "ok" and funannotate
predict ran to 0 valid BUSCO models before failing. --abs-min-bp 0 disables
this gate.

verdict is "no_training_contigs" when fewer than --min-training-contigs
contigs reach --min-contig-len bp. This is a SEPARATE, unconditional gate --
it does not require "small" to also be true. It exists because total
assembled bp says nothing about the contig-length distribution: a 15Mb
assembly shredded into 5000 contigs averaging ~3kb is not "small" by the
--min-bp gate, but has zero contigs at GeneMark-ES's own hardcoded
--min_contig 10000 training-split cutoff, so training.fna is never built and
GeneMark-ES/ET fails outright (not just a weak model). N50 being low is a
correlate of this but doesn't guarantee it -- this gate checks the literal
count directly. (Xylaria_striata_RK1-1 / GCA_002749545.1_ASM274954v1,
2026-09-16: 15,197,463 bp, 5293 contigs, N50 2832, 0 contigs >=10kb --
total_bp cleared the old --min-bp 8000000 gate so verdict was "ok" and
GeneMark-ES burned a full attempt before failing.)
--min-contig-len 0 or --min-training-contigs 0 disables this gate.

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
    ap.add_argument("--abs-min-bp", type=int, default=0,
                     help="total assembled bp below this = 'too_small', unconditionally (does not "
                          "require fragmentation too; 0 disables this gate)")
    ap.add_argument("--max-n50", type=int, default=0,
                     help="N50 below this = 'fragmented' (0 disables this gate)")
    ap.add_argument("--max-contigs", type=int, default=0,
                     help="contig count above this = 'fragmented' (0 disables this gate)")
    ap.add_argument("--min-contig-len", type=int, default=0,
                     help="length threshold for a 'training-length' contig, e.g. 10000 to match "
                          "GeneMark-ES's own --min_contig cutoff (0 disables this gate)")
    ap.add_argument("--min-training-contigs", type=int, default=0,
                     help="minimum number of contigs >= --min-contig-len required, else verdict "
                          "'no_training_contigs' (0 disables this gate)")
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

    if args.min_contig_len > 0 and args.min_training_contigs > 0:
        training_contigs = sum(1 for l in lengths if l >= args.min_contig_len)
        if training_contigs < args.min_training_contigs:
            verdict = "no_training_contigs"

    if args.abs_min_bp > 0 and total_bp < args.abs_min_bp:
        verdict = "too_small"

    line = f"{total_bp}\t{contigs}\t{n50}\t{verdict}"
    if args.report_repeat_pct:
        repeat_pct = (100.0 * lower / total_bp) if total_bp else 0.0
        line += f"\t{repeat_pct:.2f}"
    print(line)


if __name__ == "__main__":
    sys.exit(main())
