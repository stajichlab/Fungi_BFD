#!/usr/bin/env python3
"""Report files under Fungi_BFD_runs that do not belong to any row of a samples.csv.

READ ONLY. It never moves or deletes anything. It writes a TSV report; a person
reviews the report, and a separate step (find_replaced_assemblies.py --script for
replaced assemblies, or a script built from a reviewed report) moves files into a
dated quarantine folder.

What is checked, against the samples.csv given in --samples (normally the
regenerated candidate, i.e. the state after install):
    input_clean_genomes/<ASMID>.*          keyed by ASMID
    genome_annotation/<out>                keyed by species+strain (`out`)
    genemark_store/<out>
    genome_annotation_training/<out>

Why each item is stale (CATEGORY). The script looks the item up in older
samples.csv files (--history), the curation lists, and the new samples.csv:
    replaced               an assembly with the same accession number (other version)
                           is in the new samples.csv
    renamed                the same ASMID is in the new samples.csv under a different `out`
    suppressed_or_excluded every ASMID that used this key is in suppress.txt or exclude_asmids.txt
    dropped_other          the ASMID was in an old samples.csv, and is not in the new one,
                           and none of the reasons above apply
    no_history             the key never appears in any samples.csv we have: investigate by hand

ACTION_HINT is advice only:
    genemark_store folders hold GeneMark models that predict reuses and that are
    slow to rebuild, so keep them unless the genome is replaced or renamed.

Do not run `ls` or `find` on genome_annotation yourself: it has about 23,600
folders and takes over a minute. This script lists each store once.
"""
import argparse
import collections
import csv
import glob
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from create_samples_file import make_sample_tag  # noqa: E402

DEF_RUNS = os.path.abspath(os.path.join(HERE, '..', '..', 'Fungi_BFD_runs'))
CODE = os.path.abspath(os.path.join(HERE, '..'))
OUT_STORES = ('genome_annotation', 'genemark_store', 'genome_annotation_training')
# Everything after the ASMID: '.purge...' (purge copies and their FCS-GX reports),
# '.masked.fasta.gz', '.fa.gz', '.fasta.gz'. ASMIDs themselves contain dots (GCA_000.1_name).
INPUT_SUFFIX = re.compile(r'(\.purge.*|\.masked\.fasta\.gz|\.fa\.gz|\.fasta\.gz)$')


def acc_base(asmid):
    m = re.match(r'^(GC[AF]_\d+)\.\d+', asmid)
    return m.group(1) if m else None


def read_rows(path):
    try:
        with open(path, newline='') as f:
            return list(csv.DictReader(f))
    except Exception as e:  # unreadable history file: skip it
        print(f'[stale] skipping {path}: {e}', file=sys.stderr)
        return []


def load_ids(path):
    ids = set()
    if path and os.path.exists(path):
        with open(path) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    ids.add(re.split(r'[\t,]', line)[0].strip())
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--samples', required=True, help='samples.csv that defines what is current')
    ap.add_argument('--runs', default=DEF_RUNS)
    ap.add_argument('--history', nargs='*', default=None,
                    help='older samples.csv files (default: runs samples.csv plus samples*.csv in Fungi_BFD)')
    ap.add_argument('--suppress', default=os.path.join(CODE, 'data/curation/suppress.txt'))
    ap.add_argument('--exclude', default=os.path.join(CODE, 'data/curation/exclude_asmids.txt'))
    ap.add_argument('--stores', nargs='*', default=['input_clean_genomes', *OUT_STORES])
    ap.add_argument('--report', default='stale_report.tsv')
    args = ap.parse_args()

    new = {r['ASMID']: r for r in read_rows(args.samples)}
    new_outs = collections.defaultdict(list)
    for a, r in new.items():
        new_outs[make_sample_tag(r['SPECIES'], r['STRAIN'])].append(a)
    new_by_base = collections.defaultdict(list)
    for a in new:
        new_by_base[acc_base(a)].append(a)

    hist_files = args.history
    if hist_files is None:
        hist_files = [os.path.join(args.runs, 'samples.csv')]
        for pat in ('samples*.csv', 'samples.csv.*', 'backup_samples.csv', 'new_samples.csv'):
            hist_files += glob.glob(os.path.join(CODE, pat))
    hist_files = [p for p in dict.fromkeys(hist_files) if os.path.isfile(p)]
    hist_out = collections.defaultdict(set)    # out -> ASMIDs seen in any history file
    hist_asm = set()
    for p in hist_files:
        for r in read_rows(p):
            if r.get('ASMID') and r.get('SPECIES') is not None:
                hist_out[make_sample_tag(r['SPECIES'], r.get('STRAIN', ''))].add(r['ASMID'])
                hist_asm.add(r['ASMID'])
    print(f'[stale] new samples {len(new)}; history files {len(hist_files)}; '
          f'history ASMIDs {len(hist_asm)}', file=sys.stderr)

    suppressed = load_ids(args.suppress)
    excluded = load_ids(args.exclude)
    listed = suppressed | excluded

    def asm_reason(a):
        """Category for one ASMID that is not in the new samples.csv."""
        reps = [x for x in new_by_base.get(acc_base(a), []) if x != a]
        if reps:
            return 'replaced', ','.join(reps)
        if a in listed:
            return 'suppressed_or_excluded', 'in suppress/exclude list'
        if a in hist_asm:
            return 'dropped_other', 'in an old samples.csv, not in the new one'
        return 'no_history', ''

    rows = []
    for store in args.stores:
        path = os.path.join(args.runs, store)
        if not os.path.isdir(path):
            print(f'[stale] {store}: not found, skipped', file=sys.stderr)
            continue
        if store == 'input_clean_genomes':
            seen = set()
            for e in os.scandir(path):
                if not e.name.startswith('GC') or e.is_dir():
                    continue
                a = INPUT_SUFFIX.sub('', e.name)
                if a in new or a in seen:
                    continue
                seen.add(a)
                cat, detail = asm_reason(a)
                rows.append((store, a, cat, detail, 'quarantine after review'))
            # clean/ subfolder holds purge copies keyed by ASMID
            cdir = os.path.join(path, 'clean')
            if os.path.isdir(cdir):
                for e in os.scandir(cdir):
                    a = INPUT_SUFFIX.sub('', e.name)
                    if e.name.startswith('GC') and a not in new and a not in seen:
                        seen.add(a)
                        cat, detail = asm_reason(a)
                        rows.append((store + '/clean', a, cat, detail, 'quarantine after review'))
            continue
        for e in os.scandir(path):
            if not e.is_dir() or e.name.startswith(('_', '.')) or e.name in new_outs:
                continue
            hist = sorted(hist_out.get(e.name, ()))
            if not hist:
                cat, detail = 'no_history', ''
            else:
                moved = [(a, [o for o, aa in new_outs.items() if a in aa]) for a in hist if a in new]
                if moved:
                    cat, detail = 'renamed', '; '.join(f'{a} -> {o[0]}' for a, o in moved if o)
                else:
                    reasons = [asm_reason(a) for a in hist]
                    cats = {c for c, _ in reasons}
                    if 'replaced' in cats:
                        cat = 'replaced'
                    elif cats == {'suppressed_or_excluded'}:
                        cat = 'suppressed_or_excluded'
                    else:
                        cat = 'dropped_other'
                    detail = '; '.join(f'{a}:{c}' for a, (c, _) in zip(hist, reasons))[:300]
            hint = 'quarantine after review'
            if store == 'genemark_store' and cat not in ('replaced', 'renamed'):
                hint = 'keep (models reusable) unless genome is dropped for good'
            if cat == 'no_history':
                hint = 'investigate by hand'
            rows.append((store, e.name, cat, detail, hint))

    rows.sort()
    with open(args.report, 'w', newline='') as f:
        w = csv.writer(f, delimiter='\t', lineterminator='\n')
        w.writerow(['STORE', 'KEY', 'CATEGORY', 'DETAIL', 'ACTION_HINT'])
        w.writerows(rows)
    tally = collections.Counter((r[0], r[2]) for r in rows)
    print(f'[stale] {len(rows)} stale items -> {args.report}', file=sys.stderr)
    for (store, cat), n in sorted(tally.items()):
        print(f'[stale]   {store:28s} {cat:24s} {n}', file=sys.stderr)


if __name__ == '__main__':
    main()
