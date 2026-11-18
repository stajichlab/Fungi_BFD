#!/usr/bin/env python3
"""One-way average amino-acid identity (AAI) from DIAMOND best hits, per query x reference genome.

Input: a directory of DIAMOND outfmt-6 files named <QUERY>__<REF>.tsv with columns
  qseqid sseqid pident length qlen slen evalue bitscore
(run with --max-target-seqs 1, so one row per query protein that has a hit; extra HSP rows of the
same pair are ignored by keeping the first row per qseqid), plus the panel TSV from select_panel.py.

A query protein counts when evalue <= 1e-5 and the alignment covers >= 50% of the query length.
AAI = mean pident over those proteins. frac_hit = counted proteins / proteins in the query proteome.

Output TSV, one row per (query, reference), sorted by query then AAI descending, with the reference's
taxonomy, and a per-query rank.

Usage: aai_summary.py --panel panel.tsv --hits hits/ --nprot nprot.tsv --out aai.tsv
"""
import argparse
import csv
import os
import sys

EVALUE_MAX = 1e-5
QCOV_MIN = 0.5


def aai_from_hits(lines):
    """Return (n_counted, mean_pident) from DIAMOND outfmt-6 lines (first row per qseqid)."""
    seen, ids = set(), []
    for line in lines:
        f = line.rstrip('\n').split('\t')
        if len(f) < 8:
            continue
        q = f[0]
        if q in seen:
            continue
        seen.add(q)
        pident, length, qlen, evalue = float(f[2]), int(f[3]), int(f[4]), float(f[6])
        if evalue <= EVALUE_MAX and length >= QCOV_MIN * qlen:
            ids.append(pident)
    if not ids:
        return 0, float('nan')
    return len(ids), sum(ids) / len(ids)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--panel', required=True)
    ap.add_argument('--hits', required=True)
    ap.add_argument('--nprot', required=True, help='TSV LOCUSTAG<TAB>n_proteins')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    with open(a.panel) as fh:
        panel = {r['LOCUSTAG']: r for r in csv.DictReader(fh, delimiter='\t')}
    with open(a.nprot) as fh:
        nprot = {l.split('\t')[0]: int(l.split('\t')[1]) for l in fh if l.strip()}
    rows = []
    for fn in sorted(os.listdir(a.hits)):
        if not fn.endswith('.tsv') or '__' not in fn:
            continue
        q, ref = fn[:-4].split('__', 1)
        if q == ref:
            continue
        with open(os.path.join(a.hits, fn)) as fh:
            n, aai = aai_from_hits(fh)
        p = panel[ref]
        rows.append({'query': q, 'query_species': panel[q]['SPECIES'], 'ref': ref, 'ref_role': p['role'],
                     'ref_species': p['SPECIES'], 'ref_class': p['CLASS'], 'ref_order': p['ORDER'],
                     'ref_family': p['FAMILY'], 'ref_genus': p['GENUS'], 'n_counted': n,
                     'frac_hit': round(n / nprot[q], 4), 'aai': round(aai, 2)})
    if not rows:
        sys.exit('no hit files found')
    rows.sort(key=lambda r: (r['query'], -(r['aai'] if r['aai'] == r['aai'] else -1)))
    rank = {}
    for r in rows:
        rank[r['query']] = rank.get(r['query'], 0) + 1
        r['rank'] = rank[r['query']]
    cols = ['query', 'query_species', 'rank', 'ref', 'ref_role', 'ref_species', 'ref_class', 'ref_order',
            'ref_family', 'ref_genus', 'n_counted', 'frac_hit', 'aai']
    with open(a.out, 'w', newline='') as o:
        w = csv.DictWriter(o, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print(f'{len(rows)} query x reference rows -> {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
