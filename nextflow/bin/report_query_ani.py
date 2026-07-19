#!/usr/bin/env python3
"""
report_query_ani.py — Classify taxonomically-unassigned ("query") genomes against
their named ("reference") relatives from a single query_ANI.nf group.

Input is the asymmetric ANI TSV produced by `skani dist` (query<TAB>reference<TAB>ANI)
and the group's names TSV (filename, asmid, genus, species, strain, role), where
role is 'query' (missing --query_rank, e.g. no GENUS) or 'reference'.

Because skani dist only emits pairs above --skani_min_af, a query genome with no
line at all is genuinely unmatched at that alignment-fraction floor, not just
below the ANI thresholds — it is reported as tier 'no_alignment'.

Tiers (using the same --cluster-threshold / --outlier-threshold as report_ani.py):
  same_genus_high_confidence  best_ani >= cluster_threshold
  closely_related_review      outlier_threshold <= best_ani < cluster_threshold
  no_close_match              0 < best_ani < outlier_threshold
  no_alignment                no reference passed --skani_min_af at all
"""

import argparse
import os
import sys
from collections import defaultdict


def parse_ani_tsv(path):
    """Return dict[query_filename] = list of (ref_filename, ani), best-first."""
    hits = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            q, r = os.path.basename(parts[0]), os.path.basename(parts[1])
            try:
                ani = float(parts[2])
            except ValueError:
                continue
            if q == r:
                continue
            hits[q].append((r, ani))
    for q in hits:
        hits[q].sort(key=lambda x: x[1], reverse=True)
    return hits


def parse_names_tsv(path):
    """Return (queries: dict[filename]->info, refs: dict[filename]->info)."""
    queries, refs = {}, {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        cols = {name.strip().lower(): i for i, name in enumerate(header)}

        def get(parts, name):
            i = cols.get(name)
            return parts[i].strip() if i is not None and i < len(parts) else ''

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts or not parts[0]:
                continue
            info = {
                'asmid':   get(parts, 'asmid'),
                'genus':   get(parts, 'genus'),
                'species': get(parts, 'species'),
                'strain':  get(parts, 'strain'),
            }
            (queries if get(parts, 'role') == 'query' else refs)[parts[0]] = info
    return queries, refs


def classify(best_ani, cluster_thresh, outlier_thresh):
    if best_ani is None:
        return 'no_alignment'
    if best_ani >= cluster_thresh:
        return 'same_genus_high_confidence'
    if best_ani >= outlier_thresh:
        return 'closely_related_review'
    return 'no_close_match'


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True, help='query_ANI.nf asymmetric ANI TSV')
    ap.add_argument('--names', required=True, help='group names TSV (filename, asmid, genus, species, strain, role)')
    ap.add_argument('--group', required=True)
    ap.add_argument('--level', default='CLASS', help='Compare level label')
    ap.add_argument('--cluster-threshold', type=float, default=95.0)
    ap.add_argument('--outlier-threshold', type=float, default=90.0)
    ap.add_argument('--top-n', type=int, default=10, help='Reference hits shown per query in the text report')
    ap.add_argument('--report', required=True, help='Human-readable report output')
    ap.add_argument('--calls', required=True, help='Machine-readable per-query classification TSV')
    args = ap.parse_args()

    hits = parse_ani_tsv(args.input)
    queries, refs = parse_names_tsv(args.names)

    calls = []
    for qfile in sorted(queries):
        qinfo = queries[qfile]
        qhits = hits.get(qfile, [])
        best_ref, best_ani = (qhits[0] if qhits else (None, None))
        tier = classify(best_ani, args.cluster_threshold, args.outlier_threshold)
        n_ge_outlier = sum(1 for _, a in qhits if a >= args.outlier_threshold)
        rinfo = refs.get(best_ref, {}) if best_ref else {}
        calls.append({
            'group': args.group,
            'level': args.level,
            'query_asmid': qinfo['asmid'],
            'query_species': qinfo['species'],
            'query_strain': qinfo['strain'],
            'best_ref_asmid': rinfo.get('asmid', ''),
            'best_ref_genus': rinfo.get('genus', ''),
            'best_ref_species': rinfo.get('species', ''),
            'best_ani': f"{best_ani:.2f}" if best_ani is not None else '',
            'n_refs_ge_outlier': n_ge_outlier,
            'tier': tier,
        })

    with open(args.calls, 'w') as out:
        cols = ['group', 'level', 'query_asmid', 'query_species', 'query_strain',
                'best_ref_asmid', 'best_ref_genus', 'best_ref_species',
                'best_ani', 'n_refs_ge_outlier', 'tier']
        out.write('\t'.join(cols) + '\n')
        for c in calls:
            out.write('\t'.join(str(c[k]) for k in cols) + '\n')

    with open(args.report, 'w') as out:
        out.write(f"=== Query ANI Report: {args.group} ===\n")
        out.write(f"Compare level      : {args.level}\n")
        out.write(f"Query genomes      : {len(queries)}\n")
        out.write(f"Reference genomes  : {len(refs)}\n")
        out.write(f"Cluster threshold  : {args.cluster_threshold}%\n")
        out.write(f"Outlier threshold  : {args.outlier_threshold}%\n\n")

        for qfile in sorted(queries):
            qinfo = queries[qfile]
            qhits = hits.get(qfile, [])
            call = next(c for c in calls if c['query_asmid'] == qinfo['asmid'])
            label = ' '.join(p for p in (qinfo['species'], qinfo['strain']) if p)
            out.write(f"Query: {label}  [{qinfo['asmid']}]  -> tier={call['tier']}\n")
            if not qhits:
                out.write("  (no reference passed the alignment-fraction floor)\n\n")
                continue
            for rfile, ani in qhits[:args.top_n]:
                rinfo = refs.get(rfile, {})
                rlabel = ' '.join(p for p in (rinfo.get('species', ''), rinfo.get('strain', '')) if p)
                out.write(f"  {ani:6.2f}%  {rlabel}  [{rinfo.get('asmid', '')}]  genus={rinfo.get('genus', '')}\n")
            out.write("\n")

    print(f"Wrote {args.report} and {args.calls}", file=sys.stderr)


if __name__ == '__main__':
    main()
