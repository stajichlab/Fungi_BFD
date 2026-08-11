#!/usr/bin/env python3
"""Filter a (multi-)FASTA by minimum sequence length and alphabet size.

Reads FASTA from stdin (or --input) and writes to stdout (or --output) only the
records whose sequence length is >= --len. Used by funannotate.nf:

    pigz -dc genome.fna.gz | clean_genome_fa.py --len 2000 > genome.fa
    cat purged.fasta       | clean_genome_fa.py --len 2000 > asmid.fa

Additionally, contigs whose sequence contains fewer than --min-alphabet (default
4) distinct nucleotide characters are dropped with a warning. This mirrors
funannotate's analyzeAssembly() check ("Found N bad contigs, where alphabet is
less than 4 [this should not happen]"), which abort predict unless --force is
given. A contig of only A/G/T (no C), for example, is flagged exactly like
funannotate would be.

Stdlib-only (no Biopython) so it runs in any environment. The header is reduced
to its first whitespace-delimited token (">ABCDE This is more desc" -> "ABCDE",
matching Biopython's record.id); output is wrapped at --width columns
(0 = no wrapping).
"""
import argparse
import sys


def fasta_records(handle):
    """Yield (header, sequence) tuples from an open FASTA handle."""
    header = None
    chunks = []
    for line in handle:
        line = line.rstrip("\n")
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(chunks)
            # Keep only the first whitespace-delimited token as the ID, e.g.
            # ">ABCDE This is more desc" -> "ABCDE" (matches Biopython record.id).
            parts = line[1:].split(None, 1)
            header = parts[0] if parts else ""
            chunks = []
        elif line:
            chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def alphabet_counts(seq):
    """Return {uppercase_char: count} of distinct characters observed.

    Mirrors funannotate's analyzeAssembly()/Sortbysize() character census (count
    every distinct character, not just canonical ACGT), so contigs funannotate
    would flag as suspect are flagged identically here.
    """
    chars = {}
    for nuc in seq:
        nuc = nuc.upper()
        chars[nuc] = chars.get(nuc, 0) + 1
    return chars


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--len", "-l", type=int, default=1,
                    help="minimum contig length to keep (default: 1)")
    ap.add_argument("--min-alphabet", type=int, default=4,
                    help="minimum distinct nucleotide characters a contig must "
                         "contain to be kept; mirrors funannotate's "
                         "'alphabet is less than 4' suspect check (default: 4, "
                         "set to 1 to disable)")
    ap.add_argument("--input", "-i", default="-",
                    help="input FASTA (default: stdin)")
    ap.add_argument("--output", "-o", default="-",
                    help="output FASTA (default: stdout)")
    ap.add_argument("--width", "-w", type=int, default=60,
                    help="line-wrap width for sequence; 0 disables wrapping (default: 60)")
    args = ap.parse_args()

    fin = sys.stdin if args.input == "-" else open(args.input)
    fout = sys.stdout if args.output == "-" else open(args.output, "w")

    kept = 0
    dropped = 0
    dropped_alphabet = 0
    low_alphabet = []
    try:
        for header, seq in fasta_records(fin):
            bad_alphabet = False
            if args.min_alphabet > 1:
                chars = alphabet_counts(seq)
                if len(chars) < args.min_alphabet:
                    bad_alphabet = True
                    dropped_alphabet += 1
                    low_alphabet.append((header, len(seq), chars))
            if bad_alphabet:
                continue
            if len(seq) >= args.len:
                kept += 1
                fout.write(f">{header}\n")
                if args.width and args.width > 0:
                    for i in range(0, len(seq), args.width):
                        fout.write(seq[i:i + args.width] + "\n")
                else:
                    fout.write(seq + "\n")
            else:
                dropped += 1
    finally:
        if fin is not sys.stdin:
            fin.close()
        if fout is not sys.stdout:
            fout.close()

    for header, length, chars in low_alphabet:
        comp = ", ".join(f"{n}:{c:,}" for n, c in sorted(chars.items()))
        sys.stderr.write(
            f"[clean_genome_fa] WARNING: dropping {header} (len {length:,} bp): "
            f"only {len(chars)} distinct nucleotide characters "
            f"(< {args.min_alphabet}): {comp}. funannotate predict would reject "
            f"this contig ('alphabet is less than 4') unless --force is passed\n")

    sys.stderr.write(
        f"[clean_genome_fa] kept {kept} contigs >= {args.len} bp; "
        f"dropped {dropped} shorter contigs; "
        f"dropped {dropped_alphabet} contigs with < {args.min_alphabet} "
        f"distinct nucleotide characters\n")


if __name__ == "__main__":
    main()
