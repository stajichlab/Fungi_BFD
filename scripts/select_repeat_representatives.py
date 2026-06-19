#!/usr/bin/env python3
"""Select one representative genome per species for EarlGrey repeat-library construction.

A species qualifies for the expensive EarlGrey treatment only if its chosen
representative genome exceeds a size cutoff (default 200 Mb). The representative is
the best assembly in the species, preferring a RefSeq (GCF_) accession, then the
highest N50, then the fewest contigs.

For each qualifying species the output lists the representative ASMID plus all
conspecific strain ASMIDs that have an existing clean genome in --genome-dir; the
EarlGrey library built on the representative is later applied to every member.

Inputs:
  - samples.csv         (ASMID -> SPECIES)
  - tables/asm_stats.tsv.gz  (per-ASMID total_length_bp, N50_bp, contig_count)

Output (CSV): SPECIES,REP_ASMID,REP_SIZE_MB,N_MEMBERS,MEMBER_ASMIDS
  MEMBER_ASMIDS is a ';'-joined list of non-representative strain ASMIDs.
"""

import argparse
import csv
import gzip
import os
import sys
from collections import defaultdict


def load_species_map(samples_path):
    """Return {ASMID: SPECIES} from samples.csv."""
    mapping = {}
    with open(samples_path, newline="") as fh:
        for row in csv.DictReader(fh):
            asmid = (row.get("ASMID") or "").strip()
            species = (row.get("SPECIES") or "").strip()
            if asmid and species:
                mapping[asmid] = species
    return mapping


def load_asm_stats(stats_path):
    """Return {ASMID: (size_bp, n50_bp, contig_count)} from asm_stats.tsv(.gz)."""
    opener = gzip.open if stats_path.endswith(".gz") else open
    stats = {}
    with opener(stats_path, "rt", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asmid = (row.get("ASMID") or "").strip()
            if not asmid:
                continue
            try:
                size = int(row["total_length_bp"])
                n50 = int(row["N50_bp"])
                contigs = int(row["contig_count"])
            except (KeyError, ValueError):
                continue
            stats[asmid] = (size, n50, contigs)
    return stats


def rep_sort_key(entry):
    """Rank assemblies: RefSeq first, then highest N50, then fewest contigs.

    entry = (asmid, size, n50, contigs). Returned tuple is sorted descending,
    so larger is better: is_refseq (1/0), n50, -contigs.
    """
    asmid, _size, n50, contigs = entry
    is_refseq = 1 if asmid.startswith("GCF_") else 0
    return (is_refseq, n50, -contigs)


def main():
    parser = argparse.ArgumentParser(
        description="Pick per-species representative genomes (>cutoff) for EarlGrey masking",
    )
    parser.add_argument("--samples", default="samples.csv", help="samples.csv [samples.csv]")
    parser.add_argument("--asm-stats", default="tables/asm_stats.tsv.gz",
                        help="per-ASMID assembly stats TSV(.gz) [tables/asm_stats.tsv.gz]")
    parser.add_argument("--genome-dir", default="input_clean_genomes",
                        help="dir with clean <ASMID><suffix> genomes [input_clean_genomes]")
    parser.add_argument("--genome-suffix", default=".fa",
                        help="clean genome filename suffix [.fa]")
    parser.add_argument("--cutoff-mb", type=float, default=200.0,
                        help="representative size cutoff in Mb [200]")
    parser.add_argument("-o", "--output", default="misc/repeat_representatives.csv",
                        help="output CSV [misc/repeat_representatives.csv]")
    parser.add_argument("-v", "--debug", action="store_true", help="verbose stderr logging")
    args = parser.parse_args()

    cutoff_bp = args.cutoff_mb * 1_000_000

    species_of = load_species_map(args.samples)
    stats = load_asm_stats(args.asm_stats)

    # Group assemblies (with stats) by species.
    by_species = defaultdict(list)
    for asmid, species in species_of.items():
        if asmid in stats:
            size, n50, contigs = stats[asmid]
            by_species[species].append((asmid, size, n50, contigs))
        elif args.debug:
            print(f"[skip] {asmid}: no entry in {args.asm_stats}", file=sys.stderr)

    def clean_exists(asmid):
        return os.path.isfile(os.path.join(args.genome_dir, f"{asmid}{args.genome_suffix}"))

    rows = []
    n_species_over = 0
    for species, entries in by_species.items():
        rep = sorted(entries, key=rep_sort_key, reverse=True)[0]
        rep_asmid, rep_size = rep[0], rep[1]
        if rep_size <= cutoff_bp:
            continue
        n_species_over += 1

        if not clean_exists(rep_asmid):
            print(f"[warn] {species}: representative {rep_asmid} has no clean genome "
                  f"in {args.genome_dir} — skipping species", file=sys.stderr)
            continue

        members = []
        for asmid, _size, _n50, _contigs in entries:
            if asmid == rep_asmid:
                continue
            if clean_exists(asmid):
                members.append(asmid)
            elif args.debug:
                print(f"[skip] {species}: member {asmid} has no clean genome", file=sys.stderr)

        rows.append({
            "SPECIES": species,
            "REP_ASMID": rep_asmid,
            "REP_SIZE_MB": f"{rep_size / 1_000_000:.1f}",
            "N_MEMBERS": str(len(members)),
            "MEMBER_ASMIDS": ";".join(members),
        })

    rows.sort(key=lambda r: r["SPECIES"])

    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["SPECIES", "REP_ASMID", "REP_SIZE_MB", "N_MEMBERS", "MEMBER_ASMIDS"])
        writer.writeheader()
        writer.writerows(rows)

    total_strains = len(rows) + sum(int(r["N_MEMBERS"]) for r in rows)
    print(f"Species with representative > {args.cutoff_mb:g} Mb: {n_species_over}", file=sys.stderr)
    print(f"Species written (rep clean genome present): {len(rows)}", file=sys.stderr)
    print(f"Total strains (reps + members) to mask: {total_strains}", file=sys.stderr)
    print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main()
