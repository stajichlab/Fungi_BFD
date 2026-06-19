#!/usr/bin/env python3
"""
mash_components.py — Pre-cluster genomes into connected components from a mash
all-vs-all distance table, for the fastANI prefilter cascade in compare_ANI.nf.

mash dist output columns:
  reference  query  distance  p-value  shared-hashes

Two genomes are joined by an edge when their mash-estimated ANI
(100 * (1 - distance)) is >= --ani. Connected components are found with union-find.
Only components with at least --min-size members are emitted; singletons are
dropped (downstream they surface as outliers via the names file, never having
been compared by fastANI).

Output (TSV, no header) to stdout — one line per genome in a kept component:
  component_id  <TAB>  genome_filename

This lets the workflow run fastANI all-vs-all *within* each component only,
turning a clade-wide O(N^2) comparison into sum_i O(n_i^2).
"""

import argparse
import sys
from collections import defaultdict


class UnionFind:
    def __init__(self):
        self.parent = {}
        self.rank = {}

    def add(self, x):
        if x not in self.parent:
            self.parent[x] = x
            self.rank[x] = 0

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


def basename(path):
    return path.rsplit('/', 1)[-1]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--input', required=True, help='mash dist TSV (ref query dist p shared)')
    ap.add_argument('--ani', type=float, default=80.0,
                    help='mash-ANI%% floor to join two genomes (default 80.0)')
    ap.add_argument('--min-size', type=int, default=2,
                    help='minimum members for a component to be emitted (default 2)')
    args = ap.parse_args()

    uf = UnionFind()
    with open(args.input) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            ref = basename(parts[0])
            qry = basename(parts[1])
            try:
                dist = float(parts[2])
            except ValueError:
                continue  # skip any header/garbage line
            # Every genome that appears is a node, even if it joins nothing.
            uf.add(ref)
            uf.add(qry)
            if ref == qry:
                continue
            ani = 100.0 * (1.0 - dist)
            if ani >= args.ani:
                uf.union(ref, qry)

    comps = defaultdict(list)
    for node in uf.parent:
        comps[uf.find(node)].append(node)

    # Stable, size-descending component ids
    kept = sorted(
        (members for members in comps.values() if len(members) >= args.min_size),
        key=lambda m: (-len(m), sorted(m)[0]),
    )

    n_emitted = 0
    for cid, members in enumerate(kept, 1):
        for genome in sorted(members):
            sys.stdout.write(f"{cid}\t{genome}\n")
            n_emitted += 1

    print(
        f"mash_components: {len(uf.parent)} genomes -> "
        f"{len(kept)} components (>= {args.min_size}), {n_emitted} genomes kept "
        f"at ANI >= {args.ani}%",
        file=sys.stderr,
    )


if __name__ == '__main__':
    main()
