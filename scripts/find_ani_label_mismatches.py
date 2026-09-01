#!/usr/bin/env python3
"""
Proactive, dataset-wide scan for genomes whose species label likely does not
match the ANI evidence -- the KCTC_13826BP/MRD-KRBAY class of failure
(nextflow/docs/DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md, Option 7a).

Every genome compared by `compare_ani` is pre-grouped by its own claimed
SPECIES before any ANI is computed, so a pairwise ANI table can never show two
genomes with *different* claimed species -- the interesting signal is instead
a claimed-species group whose ANI graph is NOT a single connected component at
the pipeline's own cluster threshold. A group that splits into a dominant
majority cluster plus one or more small, ANI-isolated minority clusters is
exactly what KCTC_13826BP + MRD-KRBAY looked like inside "Saccharomyces
cerevisiae": Cluster 2, N=2, <90% ANI to the 1,305-genome majority cluster.

This reuses ANI data the pipeline has already computed (`compare_ani`'s merged
pairwise table) -- no new Nextflow compute, no new sketching/comparison.

Usage:
    python scripts/find_ani_label_mismatches.py \\
        --ani-merged results/ANI/skani/SPECIES/all_pairs_merged.tsv \\
        --samples-csv samples.csv \\
        [--cluster-threshold 95.0] [--min-group-size 3] [--outfile FILE]

Output (TSV, to --outfile or stdout):
    species  group_size  n_components  majority_size  minority_asmids
    minority_size  max_ani_minority_to_majority

Exit code is the number of flagged species groups (0 = nothing found).
"""

import argparse
import csv
import os
import sys
from collections import defaultdict


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--ani-merged", default="results/ANI/skani/SPECIES/all_pairs_merged.tsv",
                   help="Merged pairwise ANI table (query, ref, ANI columns)")
    p.add_argument("--samples-csv", default="samples.csv",
                   help="samples.csv, for ASMID -> declared SPECIES lookup")
    p.add_argument("--cluster-threshold", type=float, default=95.0,
                   help="ANI%% at/above which two genomes are same-cluster (default: 95.0, matches ani_cluster_threshold)")
    p.add_argument("--outlier-threshold", type=float, default=90.0,
                   help="ANI%% below which a minority cluster is HIGH confidence rather than REVIEW "
                        "(default: 90.0, matches profile_ANI.config's ani_outlier_threshold / README_ANI.md's "
                        "documented 'potential misidentifications, novel lineages, or contaminants' convention)")
    p.add_argument("--min-group-size", type=int, default=3,
                   help="Skip species groups with fewer than this many genomes (default: 3 -- a 2-genome group splitting is not informative)")
    p.add_argument("--outfile", default=None, help="Write TSV here instead of stdout")
    return p.parse_args()


def asmid_from_genome_name(name):
    """'query/GCA_003277715.1_ASM327771v1.fa.gz' -> 'GCA_003277715.1_ASM327771v1'.

    skani's query/ref columns may carry a directory prefix (e.g. a 'query/'
    staging dir) ahead of the genome filename; strip it via basename before
    stripping the FASTA suffix, so ASMID lookup against samples.csv works
    regardless of the path shape a given skani invocation uses.
    """
    base = os.path.basename(name)
    for suffix in (".fa.gz", ".fasta.gz", ".fa", ".fasta"):
        if base.endswith(suffix):
            return base[: -len(suffix)]
    return base


def load_species_map(samples_csv):
    m = {}
    with open(samples_csv, newline="") as f:
        for row in csv.DictReader(f):
            m[row["ASMID"]] = row["SPECIES"]
    return m


class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def main():
    args = parse_args()
    species_of = load_species_map(args.samples_csv)

    # groups[species] = {"members": set(asmid), "edges": [(a, b, ani)], "cross_max": {}}
    groups = defaultdict(lambda: {"members": set(), "edges": []})

    n_lines = n_no_species = n_species_mismatch = 0
    with open(args.ani_merged, newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        assert header[:3] == ["query", "ref", "ANI"], f"unexpected header: {header}"
        for row in reader:
            n_lines += 1
            q_asm = asmid_from_genome_name(row[0])
            r_asm = asmid_from_genome_name(row[1])
            ani = float(row[2])
            q_sp = species_of.get(q_asm)
            r_sp = species_of.get(r_asm)
            if q_sp is None or r_sp is None:
                n_no_species += 1
                continue
            if q_sp != r_sp:
                # Should not happen -- compare_ani groups by claimed species
                # before running ANI. Surface it rather than silently drop it.
                n_species_mismatch += 1
                continue
            g = groups[q_sp]
            g["members"].add(q_asm)
            g["members"].add(r_asm)
            g["edges"].append((q_asm, r_asm, ani))

    if n_species_mismatch:
        print(f"[warn] {n_species_mismatch} pairs had query/ref under different "
              f"claimed SPECIES -- compare_ani's own grouping invariant broke; "
              f"investigate before trusting this report", file=sys.stderr)
    print(f"[info] {n_lines} pairs, {n_no_species} skipped (ASMID not in "
          f"samples.csv), {len(groups)} species groups with >=1 pair", file=sys.stderr)

    flagged = []
    for species, g in groups.items():
        members = g["members"]
        if len(members) < args.min_group_size:
            continue
        uf = UnionFind(members)
        for a, b, ani in g["edges"]:
            if ani >= args.cluster_threshold:
                uf.union(a, b)
        components = defaultdict(set)
        for m in members:
            components[uf.find(m)].add(m)
        if len(components) < 2:
            continue  # single connected component -- no split, nothing to flag
        comps_sorted = sorted(components.values(), key=len, reverse=True)
        majority = comps_sorted[0]
        for minority in comps_sorted[1:]:
            max_ani = 0.0
            for a, b, ani in g["edges"]:
                if (a in minority and b in majority) or (b in minority and a in majority):
                    max_ani = max(max_ani, ani)
            flagged.append({
                "species": species,
                "group_size": len(members),
                "n_components": len(components),
                "majority_size": len(majority),
                "minority_asmids": ",".join(sorted(minority)),
                "minority_size": len(minority),
                "max_ani_minority_to_majority": f"{max_ani:.2f}",
                "confidence": "HIGH" if max_ani < args.outlier_threshold else "REVIEW",
            })

    out = open(args.outfile, "w", newline="\n") if args.outfile else sys.stdout
    fieldnames = ["species", "group_size", "n_components", "majority_size",
                  "minority_asmids", "minority_size", "max_ani_minority_to_majority",
                  "confidence"]
    writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in sorted(flagged, key=lambda r: (r["confidence"] != "HIGH", -r["group_size"])):
        writer.writerow(row)
    if args.outfile:
        out.close()

    n_high = sum(1 for f in flagged if f["confidence"] == "HIGH")
    print(f"[info] flagged {len(flagged)} minority clusters across "
          f"{len({f['species'] for f in flagged})} species groups "
          f"({n_high} HIGH confidence, {len(flagged) - n_high} REVIEW)", file=sys.stderr)
    return len(flagged)


if __name__ == "__main__":
    sys.exit(main())
