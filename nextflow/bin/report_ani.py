#!/usr/bin/env python3
"""
report_ani.py — Parse an all-vs-all ANI TSV and write a clustering + outlier report.

The input is a tab-separated file whose first three columns are:
  query  reference  ANI

This format is emitted by every --ani_method in compare_ANI.nf:
  - fastani   : native fastANI output (cols 4/5 = mappings, ignored here)
  - skani     : normalized from `skani triangle --sparse`
  - mash      : normalized from `mash dist`           (ANI = 100*(1-distance))
  - sourmash  : normalized from `sourmash compare --ani` matrix

Sparse methods (skani/sourmash/mash-prefilter) and fastANI both *omit* pairs
below their internal mapping/identity floor, so a fully divergent genome may have
NO line at all in the TSV. To count and flag such genomes as outliers we take the
authoritative genome universe from the --names file (one row per genome) when it
is provided, and fall back to the genomes seen in the pairs otherwise.

Two ANI thresholds are used:
  --cluster-threshold  (default 95.0)  Minimum ANI to place two genomes in the same cluster
  --outlier-threshold  (default 90.0)  Max ANI below which a genome is flagged as an outlier
"""

import argparse
import os
import sys
from collections import defaultdict
from datetime import date
from statistics import median


# ── Union-Find ────────────────────────────────────────────────────────────────

class UnionFind:
    def __init__(self, items):
        self.parent = {x: x for x in items}
        self.rank   = {x: 0  for x in items}

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x, y):
        px, py = self.find(x), self.find(y)
        if px == py:
            return
        if self.rank[px] < self.rank[py]:
            px, py = py, px
        self.parent[py] = px
        if self.rank[px] == self.rank[py]:
            self.rank[px] += 1

    def clusters(self):
        groups = defaultdict(list)
        for item in self.parent:
            groups[self.find(item)].append(item)
        return list(groups.values())


# ── Parsing ───────────────────────────────────────────────────────────────────

def parse_ani_tsv(path):
    """Return dict[(query_base, ref_base)] = max_ANI, excluding self-hits."""
    pairs = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            q   = os.path.basename(parts[0])
            r   = os.path.basename(parts[1])
            try:
                ani = float(parts[2])
            except ValueError:
                # header line (e.g. skani 'ANI') or malformed row — skip
                continue
            if q == r:
                continue
            key = tuple(sorted([q, r]))
            if key not in pairs or pairs[key] < ani:
                pairs[key] = ani
    return pairs


def parse_names_tsv(path):
    """Return dict[filename] = species_name."""
    names = {}
    if not path or not os.path.exists(path):
        return names
    with open(path) as fh:
        next(fh, None)  # skip header
        for line in fh:
            parts = line.rstrip('\n').split('\t', 1)
            if len(parts) == 2:
                names[parts[0]] = parts[1]
    return names


# ── Report ────────────────────────────────────────────────────────────────────

def write_report(args, pairs, names, genomes):
    cluster_thresh  = args.cluster_threshold
    outlier_thresh  = args.outlier_threshold

    # Build cluster membership via Union-Find
    uf = UnionFind(genomes)
    for (q, r), ani in pairs.items():
        if ani >= cluster_thresh:
            uf.union(q, r)

    raw_clusters = uf.clusters()
    # Sort clusters by size descending
    clusters = sorted(raw_clusters, key=len, reverse=True)

    # Identify outliers: max ANI to any other genome < outlier_thresh
    max_ani_for = {g: 0.0 for g in genomes}
    best_match  = {g: '' for g in genomes}
    for (q, r), ani in pairs.items():
        if ani > max_ani_for[q]:
            max_ani_for[q] = ani
            best_match[q]  = r
        if ani > max_ani_for[r]:
            max_ani_for[r] = ani
            best_match[r]  = q

    outliers = [g for g in genomes if max_ani_for[g] < outlier_thresh]

    with open(args.output, 'w') as out:
        # ── Header ────────────────────────────────────────────────────────────
        out.write(f"=== ANI Report: {args.group} ===\n")
        out.write(f"Compare level       : {args.level}\n")
        out.write(f"Genomes             : {len(genomes)}\n")
        out.write(f"Cluster threshold   : {cluster_thresh}%\n")
        out.write(f"Outlier threshold   : {outlier_thresh}%\n")
        out.write(f"Generated           : {date.today()}\n")
        out.write("\n")

        # ── Clusters ──────────────────────────────────────────────────────────
        out.write(f"--- Clusters (ANI >= {cluster_thresh}%) ---\n\n")

        for idx, cluster in enumerate(clusters, 1):
            # Compute within-cluster ANI statistics
            within = [
                ani for (q, r), ani in pairs.items()
                if q in cluster and r in cluster
            ]
            med_ani = f"{median(within):.2f}%" if within else "n/a"
            min_ani = f"{min(within):.2f}%"    if within else "n/a"

            out.write(f"Cluster {idx}  N={len(cluster)}  median_ANI={med_ani}  min_ANI={min_ani}\n")
            for g in sorted(cluster):
                sp = names.get(g, '')
                out.write(f"  {g:<70}  {sp}\n")
            out.write("\n")

        # ── Outliers ──────────────────────────────────────────────────────────
        out.write(f"--- Outliers (best ANI < {outlier_thresh}%) ---\n\n")

        if outliers:
            for g in sorted(outliers):
                sp    = names.get(g, '')
                bani  = f"{max_ani_for[g]:.2f}%"
                bmatch = best_match[g]
                out.write(f"  {g:<70}  {sp}\n")
                out.write(f"    best_ANI={bani}  best_match={bmatch}\n")
        else:
            out.write("  (none)\n")

    print(f"Report written to {args.output}", file=sys.stderr)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input',              required=True,  help='fastANI TSV output')
    ap.add_argument('--names',              default=None,   help='filename→species name TSV')
    ap.add_argument('--group',              required=True,  help='Group name (e.g. genus name)')
    ap.add_argument('--level',              default='GENUS',help='Compare level label')
    ap.add_argument('--cluster-threshold',  type=float, default=95.0,
                    help='ANI%% cutoff for cluster formation (default 95.0)')
    ap.add_argument('--outlier-threshold',  type=float, default=90.0,
                    help='ANI%% below which a genome is an outlier (default 90.0)')
    ap.add_argument('--output',             required=True,  help='Output report file')
    args = ap.parse_args()

    if args.outlier_threshold > args.cluster_threshold:
        print(
            f"Warning: outlier-threshold ({args.outlier_threshold}) > "
            f"cluster-threshold ({args.cluster_threshold}); "
            "all outliers will also be singletons.",
            file=sys.stderr
        )

    pairs = parse_ani_tsv(args.input)
    names = parse_names_tsv(args.names)

    # Genome universe: prefer the names file (one row per genome in the group) so
    # that genomes with no surviving pair — true outliers under sparse methods —
    # are still counted and flagged. Fall back to the genomes seen in the pairs.
    genomes_set = set(names.keys())
    for q, r in pairs:
        genomes_set.add(q)
        genomes_set.add(r)
    genomes = sorted(genomes_set)

    if not genomes:
        with open(args.output, 'w') as out:
            out.write(f"=== ANI Report: {args.group} ===\n")
            out.write("No genomes found (empty names file and no pairs).\n")
        print("Warning: empty input — no genomes or pairs.", file=sys.stderr)
        return

    write_report(args, pairs, names, genomes)


if __name__ == '__main__':
    main()
