#!/usr/bin/env python3
"""Flexible telomere repeat finder for genome assemblies.

This tool generalizes the regex-based approach of find_telomeres.py (Hiltunen
2021) to support multiple canonical and user-supplied monomer patterns,
including IUPAC degenerate bases and regex quantifiers such as ``[C]+``.

For every scaffold the tool examines the terminal ``--search-window`` bp and
reports telomere tracts that reach the scaffold terminus (``terminal=True``)
and satisfy the repeat-count and length thresholds.  By default it follows the
canonical telomere orientation: the supplied monomer is searched at the 5' end
and its reverse complement at the 3' end.  Use ``--both-ends`` to search both
monomer orientations at both scaffold ends.

Output is a tab-separated table with one row per telomere tract.  Recorded
features include repeat monomer, repeat count, total tract length, scaffold
coordinates, strand, and an inward flanking sequence.

Typical usage::

    find_telomeres.py genome.scaffolds.fa -o telomeres.tsv
    find_telomeres.py genome.scaffolds.fa --fuzzy --max-mismatch 1 -o telomeres.tsv
    find_telomeres.py genome.scaffolds.fa --patterns CCCTAA TTAGGG -o telomeres.tsv
"""

import argparse
import csv
import gzip
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Tuple

from Bio import SeqIO

DEFAULT_MONOMERS = [
    # Fungal / eukaryote telomere monomers. Listed as the forward-strand
    # 5'->3' monomer; the reverse complement is searched automatically at the
    # 3' scaffold end. Sources: Hiltunen 2021 (Rhodotorula), TeloBase
    # (Lycka et al. 2023), canonical yeast/fungal motifs.
    "TAA[C]+",          # Candida / Rhodotorula-like (Hiltunen default)
    "TTAGGG",           # Vertebrate-like; Neurospora, many ascomycetes
    "TTAGG",            # Insect-like; also some fungi
    "TTAGGGG",          # Saccharomyces-like
    "TTAGGGGT",         # Variant
    "TTAGGGTCAACA",     # Aspergillus section Flavi (TeloBase)
    "TTATTAGGG",        # Aspergillus restricti / parasitoid wasps (TeloBase)
    "TTTATTAGGG",       # Aspergillus transcarpathicus / Chrysobalanaceae (TeloBase)
    "TTATTGGGG",        # Bombus / some wasp variants (TeloBase)
    "TACAAGG",          # Some basidiomycete telomeres
    "CCCTAA",           # Reverse of TTAGGG (explicit for symmetry)
]

# IUPAC nucleotide ambiguity codes -> regex character class.
IUPAC_TO_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "U",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

# DNA complement table for regex-aware reverse complement.
_COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K", "B": "V", "D": "H",
    "H": "D", "V": "B", "N": "N",
}


@dataclass(frozen=True)
class TelomereHit:
    """One telomere tract hit."""
    scaffold: str
    end: str              # '5prime' or '3prime'
    strand: str           # '+' or '-'
    monomer: str          # User-supplied monomer (5'->3' canonical)
    repeat_count: int
    tract_length: int
    start: int            # 0-based inclusive
    end_coord: int        # 0-based exclusive
    tract_seq: str
    flank_seq: str        # Sequence inward from the telomere tract
    terminal: bool        # True if tract reaches the scaffold terminus


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genome", help="Genome FASTA file (may be gzip-compressed).")
    p.add_argument("-o", "--output", default="-",
                   help="Output TSV file (default: stdout).")
    p.add_argument("-p", "--patterns", nargs="+", default=None,
                   help="One or more monomer patterns. Defaults to a curated "
                        "fungal set. Patterns may contain IUPAC ambiguity "
                        "codes and '+' quantifiers (e.g. 'TAA[C]+').")
    p.add_argument("-n", "--min-repeats", type=int, default=2,
                   help="Minimum number of contiguous monomer repeats [2].")
    p.add_argument("-l", "--min-length", type=int, default=12,
                   help="Minimum telomere tract length in bp [12].")
    p.add_argument("-w", "--flank-window", type=int, default=500,
                   help="Length of inward flanking sequence to report [500].")
    p.add_argument("-s", "--search-window", type=int, default=10000,
                   help="Number of terminal bp to search for telomere tracts "
                        "[10000]. Increase for long subtelomeric / repeat "
                        "regions, decrease for speed.")
    p.add_argument("--fuzzy", action="store_true",
                   help="Enable fuzzy monomer alignment (slower but more "
                        "sensitive for divergent repeats).")
    p.add_argument("--max-mismatch", type=int, default=1,
                   help="Maximum mismatches per monomer in fuzzy mode [1].")
    p.add_argument("--max-indel", type=int, default=0,
                   help="Maximum indels per monomer in fuzzy mode [0].")
    p.add_argument("--allow-internal", action="store_true",
                   help="Also report telomere-like tracts that do not reach "
                        "the scaffold terminus (default: terminal-only).")
    p.add_argument("--both-ends", action="store_true",
                   help="Search both monomer orientations at both scaffold "
                        "ends (default: canonical orientation only).")
    p.add_argument("--no-rc", action="store_true",
                   help="Disable automatic reverse-complement search at the "
                        "3' end (only useful with --both-ends).")
    p.add_argument("--tsv-sequence", action="store_true", default=True,
                   help="Include tract and flank sequences in TSV output "
                        "(default: True).")
    p.add_argument("--no-tsv-sequence", dest="tsv_sequence",
                   action="store_false",
                   help="Omit tract and flank sequences from TSV output.")
    p.add_argument("--version", action="version", version="find_telomeres 1.0")
    return p.parse_args(argv)


def iupac_to_regex(pattern: str) -> str:
    """Convert a pattern containing IUPAC codes to a regex string."""
    out = []
    for ch in pattern.upper():
        if ch in IUPAC_TO_REGEX:
            out.append(IUPAC_TO_REGEX[ch])
        elif ch in "+*?[](){}|^$\\":
            out.append(ch)
        else:
            raise ValueError(f"Invalid character in pattern '{pattern}': {ch}")
    return "".join(out)


def reverse_complement_pattern(pattern: str) -> str:
    """Reverse complement a monomer pattern that may contain regex quantifiers.

    Handles IUPAC ambiguity codes and the simple quantifier form ``[BASE]+",
    as used by ``TAA[C]+``.  More complex regexes are rejected rather than
    silently mangled.
    """
    # Tokenize into: plain base, bracketed class, quantifier attached to prior.
    tokens: List[str] = []
    i = 0
    n = len(pattern)
    while i < n:
        ch = pattern[i]
        if ch == "[":
            j = pattern.find("]", i)
            if j == -1:
                raise ValueError(f"Unclosed '[' in pattern '{pattern}'")
            tokens.append(pattern[i:j + 1])
            i = j + 1
        elif ch in "+*?{}":
            if not tokens:
                raise ValueError(f"Quantifier '{ch}' without target in '{pattern}'")
            tokens[-1] = tokens[-1] + ch
            i += 1
        else:
            tokens.append(ch)
            i += 1

    # Reverse complement each token and reverse the list.
    rc_tokens = []
    for tok in tokens[::-1]:
        if tok.startswith("[") and "]" in tok:
            body = tok[1:tok.index("]")]
            quant = tok[tok.index("]") + 1:]
            rc_body = "".join(_COMPLEMENT.get(b, "N") for b in body)
            rc_tokens.append(f"[{rc_body}]{quant}")
        else:
            # Plain base possibly followed by a quantifier.
            base = tok[0]
            quant = tok[1:]
            rc_tokens.append(f"{_COMPLEMENT.get(base, 'N')}{quant}")
    return "".join(rc_tokens)


def compile_pattern(pattern: str) -> re.Pattern:
    """Compile a single monomer pattern to a regex matching one or more repeats."""
    regex = iupac_to_regex(pattern)
    try:
        return re.compile(f"({regex})+", re.IGNORECASE)
    except re.error as exc:
        raise ValueError(f"Cannot compile pattern '{pattern}' -> '{regex}': {exc}")


def extract_monomer_unit(tract_seq: str, monomer: str,
                         fuzzy: bool = False,
                         max_mismatch: int = 1,
                         max_indel: int = 0) -> Tuple[int, str]:
    """Infer the number of monomer repeats in a tract and the consensus unit."""
    if not fuzzy:
        regex = re.compile(iupac_to_regex(monomer), re.IGNORECASE)
        count = len(regex.findall(tract_seq))
        return max(count, 1), monomer

    mlen = len(monomer)
    if mlen == 0 or len(tract_seq) < mlen:
        return 1, monomer

    count = round(len(tract_seq) / mlen)
    units = [tract_seq[i:i + mlen]
             for i in range(0, len(tract_seq) - mlen + 1, mlen)]
    if not units:
        return 1, monomer

    consensus = []
    for pos in range(mlen):
        votes: Dict[str, int] = {}
        for unit in units:
            if pos < len(unit):
                base = unit[pos].upper()
                votes[base] = votes.get(base, 0) + 1
        consensus.append(max(votes.items(), key=lambda kv: kv[1])[0]
                         if votes else "N")
    return count, "".join(consensus)


def _unit_edits(unit: str, monomer: str, max_indel: int) -> int:
    """Bounded edit distance between a unit and monomer."""
    if max_indel == 0:
        return sum(1 for a, b in zip(unit, monomer) if a != b)
    m, n = len(unit), len(monomer)
    max_diff = abs(m - n)
    if max_diff > max_indel:
        return max_diff
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        curr = [i] + [0] * n
        for j in range(1, n + 1):
            cost = 0 if unit[i - 1] == monomer[j - 1] else 1
            curr[j] = min(curr[j - 1] + 1, prev[j] + 1, prev[j - 1] + cost)
        prev = curr
    return prev[n]


def fuzzy_find_telomere(seq: str, monomer: str,
                        min_repeats: int,
                        max_mismatch: int,
                        max_indel: int) -> List[Tuple[int, int, str]]:
    """Find maximal telomere tracts by sliding a fuzzy monomer window."""
    mlen = len(monomer)
    if mlen == 0:
        return []

    monomer_upper = monomer.upper()
    seq_upper = seq.upper()
    n = len(seq)
    matches: List[Tuple[int, int, str]] = []
    i = 0

    while i <= n - mlen * min_repeats:
        run_start = i
        run_end = i
        total_edits = 0
        total_units = 0

        while run_end + mlen <= n:
            unit = seq_upper[run_end:run_end + mlen]
            edits = _unit_edits(unit, monomer_upper, max_indel)
            if edits <= max_mismatch + max_indel:
                total_edits += edits
                total_units += 1
                run_end += mlen
                if total_edits / total_units > max_mismatch:
                    break
            else:
                break

        if total_units >= min_repeats:
            matches.append((run_start, run_end, seq[run_start:run_end]))
            i = run_end
        else:
            i += 1

    return matches


def regex_terminal_hits(scaffold: str, seq: str, pattern: re.Pattern,
                        monomer: str, strand: str, end: str,
                        min_repeats: int, min_length: int,
                        flank_window: int, allow_internal: bool,
                        scaffold_len: int, search_window: int,
                        ) -> Iterable[TelomereHit]:
    """Yield regex-based TelomereHit objects for one scaffold end/strand."""
    is_3prime = end == "3prime"
    search_seq = seq[::-1][:search_window] if is_3prime else seq[:search_window]

    for m in pattern.finditer(search_seq):
        tract_len = m.end() - m.start()
        if tract_len < min_length:
            continue

        tract_seq = search_seq[m.start():m.end()]
        repeat_count, inferred_unit = extract_monomer_unit(tract_seq, monomer)
        if repeat_count < min_repeats:
            continue

        if is_3prime:
            start_coord = scaffold_len - m.end()
            end_coord = scaffold_len - m.start()
            tract_seq_orig = tract_seq[::-1]
            flank_start = max(0, start_coord - flank_window)
            flank_seq = seq[flank_start:start_coord]
            terminal = end_coord == scaffold_len
        else:
            start_coord = m.start()
            end_coord = m.end()
            tract_seq_orig = tract_seq
            flank_end = min(scaffold_len, end_coord + flank_window)
            flank_seq = seq[end_coord:flank_end]
            terminal = start_coord == 0

        if not terminal and not allow_internal:
            continue

        yield TelomereHit(
            scaffold=scaffold, end=end, strand=strand,
            monomer=monomer, repeat_count=repeat_count,
            tract_length=tract_len, start=start_coord,
            end_coord=end_coord, tract_seq=tract_seq_orig,
            flank_seq=flank_seq, terminal=terminal,
        )


def fuzzy_terminal_hits(scaffold: str, seq: str, monomer: str,
                        strand: str, end: str,
                        min_repeats: int, min_length: int,
                        flank_window: int, allow_internal: bool,
                        scaffold_len: int, search_window: int,
                        max_mismatch: int, max_indel: int,
                        ) -> Iterable[TelomereHit]:
    """Yield fuzzy-alignment TelomereHit objects for one scaffold end/strand."""
    is_3prime = end == "3prime"
    search_seq = seq[::-1][:search_window] if is_3prime else seq[:search_window]

    matches = fuzzy_find_telomere(
        search_seq, monomer, min_repeats, max_mismatch, max_indel
    )

    for start, end_pos, tract_seq in matches:
        tract_len = end_pos - start
        if tract_len < min_length:
            continue
        repeat_count, inferred_unit = extract_monomer_unit(
            tract_seq, monomer, fuzzy=True,
            max_mismatch=max_mismatch, max_indel=max_indel
        )
        if repeat_count < min_repeats:
            continue

        if is_3prime:
            start_coord = scaffold_len - end_pos
            end_coord = scaffold_len - start
            tract_seq_orig = tract_seq[::-1]
            flank_start = max(0, start_coord - flank_window)
            flank_seq = seq[flank_start:start_coord]
            terminal = end_coord == scaffold_len
        else:
            start_coord = start
            end_coord = end_pos
            tract_seq_orig = tract_seq
            flank_end = min(scaffold_len, end_coord + flank_window)
            flank_seq = seq[end_coord:flank_end]
            terminal = start_coord == 0

        if not terminal and not allow_internal:
            continue

        yield TelomereHit(
            scaffold=scaffold, end=end, strand=strand,
            monomer=monomer, repeat_count=repeat_count,
            tract_length=tract_len, start=start_coord,
            end_coord=end_coord, tract_seq=tract_seq_orig,
            flank_seq=flank_seq, terminal=terminal,
        )


def open_handle(path: str):
    """Open a possibly-gzipped FASTA file."""
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def build_pattern_set(patterns: List[str], both_ends: bool, no_rc: bool
                      ) -> List[Tuple[str, str, str]]:
    """Return list of (canonical_monomer, search_monomer, strand) tuples.

    canonical_monomer: the pattern as supplied by the user.
    search_monomer: the monomer to actually search at this end.
    strand: '+' if search_monomer is the canonical monomer, '-' if it is the
            reverse complement.
    """
    result: List[Tuple[str, str, str]] = []
    seen: set = set()
    for pat in patterns:
        pat = pat.upper()
        rc = reverse_complement_pattern(pat)
        if no_rc:
            result.append((pat, pat, "+"))
            seen.add((pat, pat, "+"))
            if both_ends:
                result.append((pat, pat, "-"))
            continue

        # Canonical orientation: forward at 5', RC at 3'.
        key = (pat, pat, "+")
        if key not in seen:
            result.append(key)
            seen.add(key)
        key = (pat, rc, "-")
        if key not in seen:
            result.append(key)
            seen.add(key)

        if both_ends:
            # Also search RC at 5' and forward at 3'.
            key = (pat, rc, "+")
            if key not in seen:
                result.append(key)
                seen.add(key)
            key = (pat, pat, "-")
            if key not in seen:
                result.append(key)
                seen.add(key)
    return result


def main(argv: Optional[List[str]] = None) -> int:
    args = parse_args(argv)

    patterns = args.patterns if args.patterns else DEFAULT_MONOMERS
    pattern_set = build_pattern_set(patterns, args.both_ends, args.no_rc)

    # Pre-compile regexes for exact mode.
    regexes: Dict[str, re.Pattern] = {}
    for _canon, search_monomer, _strand in pattern_set:
        if search_monomer not in regexes:
            regexes[search_monomer] = compile_pattern(search_monomer)

    out_fh = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    writer = csv.writer(out_fh, delimiter="\t", lineterminator="\n")
    columns = [
        "scaffold", "end", "strand", "monomer", "repeat_count",
        "tract_length", "start", "end_coord", "terminal",
    ]
    if args.tsv_sequence:
        columns += ["tract_seq", "flank_seq"]
    writer.writerow(columns)

    total_hits = 0
    with open_handle(args.genome) as fh:
        for record in SeqIO.parse(fh, "fasta"):
            scaffold_id = record.id
            seq = str(record.seq)
            slen = len(seq)
            if slen < args.min_length:
                continue

            hits: List[TelomereHit] = []
            for canonical, search_monomer, strand in pattern_set:
                end = "5prime" if strand == "+" else "3prime"
                if args.fuzzy:
                    hits.extend(fuzzy_terminal_hits(
                        scaffold_id, seq, search_monomer,
                        strand, end,
                        min_repeats=args.min_repeats,
                        min_length=args.min_length,
                        flank_window=args.flank_window,
                        allow_internal=args.allow_internal,
                        scaffold_len=slen,
                        search_window=args.search_window,
                        max_mismatch=args.max_mismatch,
                        max_indel=args.max_indel,
                    ))
                else:
                    pattern = regexes[search_monomer]
                    hits.extend(regex_terminal_hits(
                        scaffold_id, seq, pattern,
                        canonical, strand, end,
                        min_repeats=args.min_repeats,
                        min_length=args.min_length,
                        flank_window=args.flank_window,
                        allow_internal=args.allow_internal,
                        scaffold_len=slen,
                        search_window=args.search_window,
                    ))

            for hit in hits:
                row = [
                    hit.scaffold, hit.end, hit.strand, hit.monomer,
                    hit.repeat_count, hit.tract_length,
                    hit.start, hit.end_coord, hit.terminal,
                ]
                if args.tsv_sequence:
                    row += [hit.tract_seq, hit.flank_seq]
                writer.writerow(row)
                total_hits += 1

    if args.output != "-":
        out_fh.close()

    print(f"# find_telomeres: {total_hits} telomere tract(s) reported",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
