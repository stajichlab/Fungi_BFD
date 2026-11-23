#!/usr/bin/env python3
"""Find assemblies in samples.csv that NCBI no longer lists, and list the files to remove.

Run this after refreshing ncbi_accessions.csv and BEFORE you regenerate samples.csv.
It compares the current samples.csv with the new accessions table.

An ASMID in samples.csv that is missing from the accessions table is "gone".
Each gone ASMID gets one class:
    version_update    same accession base (GCA_000000000) with a new version or name
    replaced_by_other a new assembly has the same species and strain
    removed_upstream  no replacement found (NCBI suppressed or withdrew it)

For each gone ASMID the script lists the paths it would clean:
    <runs>/input_clean_genomes/<ASMID>.*              (keyed by ASMID)
    <runs>/input_clean_genomes/clean/<ASMID>.*        (keyed by ASMID)
    <runs>/genome_annotation/<out>                    (keyed by species+strain)
    <runs>/genemark_store/<out>
    <runs>/genome_annotation_training/<out>
`<out>` is built like makeSampleTag in nextflow/modules/common/utils.nf.
A per-`out` folder is listed only if no kept row in the new accessions table still
uses it. When the replacement uses the SAME `out`, the folder holds results for the
OLD assembly. Nextflow -resume could reuse it, so it is flagged STALE_SAME_OUT.

This script never deletes or moves anything. It writes:
    --report   TSV, one row per gone ASMID (default: replaced_assemblies.tsv)
    --script   shell script that MOVES each path into a dated quarantine folder
               (default: none; give a path to write it). Nothing is deleted.
               Each move keeps the original relative path and appends a row to
               <quarantine>/QUARANTINE_LOG.tsv (date, ASMID, path, size in KB).
    --quarantine   quarantine folder (default: <runs>/stale_<YYYYMMDD>)
    --append-removed   append rows to data/curation/removed_asmids.csv
To delete for good, remove the quarantine folder by hand after you have compared
the old and new results.

Check `squeue -u $USER` and any live Nextflow run before you run the shell script.
Do not scan genome_annotation with ls or find: it has about 23,600 folders.
The script tests only the specific paths it builds.
"""

import argparse
import csv
import datetime
import glob
import os
import re
import shlex
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sample_sanitize import clean_species, clean_strain, backfill_strain  # noqa: E402

DEF_RUNS = os.path.abspath(os.path.join(HERE, '..', '..', 'Fungi_BFD_runs'))
DEF_NCBI = os.path.abspath(os.path.join(HERE, '..', '..', '..', '1KFG', '2026', 'NCBI_fungi'))
PER_OUT_DIRS = ('genome_annotation', 'genemark_store', 'genome_annotation_training')
PER_ASMID_DIRS = ('input_clean_genomes', os.path.join('input_clean_genomes', 'clean'))


def make_sample_tag(species, strain):
    """Port of makeSampleTag/cleanStrain (nextflow/modules/common/utils.nf)."""
    sp = (species or '').strip().replace("'", '').replace('"', '')
    st = (strain or '').split(';')[0].replace(':', ' ').strip()
    st = st.strip('*').replace('*', '-')
    joined = '_'.join(x for x in (sp, st) if x)
    return re.sub(r'[\s/#\[\]?{}]+', '_', joined)


def acc_base(asmid):
    """GCA_000000000 from GCA_000000000.1_name (None if the form is unexpected)."""
    m = re.match(r'^(GC[AF]_\d+)\.\d+', asmid)
    return m.group(1) if m else None


def read_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def load_taxonomy_species(path):
    """ASM_ACCESSION -> binomial SPECIES (first row wins)."""
    out = {}
    if not path or not os.path.exists(path):
        return out
    for row in read_csv(path):
        out.setdefault(row['ASM_ACCESSION'], (row.get('SPECIES') or '').strip())
    return out


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
    ap.add_argument('--samples', default='samples.csv', help='current samples.csv (before regeneration)')
    ap.add_argument('--accessions', default=f'{DEF_NCBI}/ncbi_accessions.csv')
    ap.add_argument('--taxonomy', default=f'{DEF_NCBI}/ncbi_accessions_taxonomy.csv')
    ap.add_argument('--exclude', default='data/curation/exclude_asmids.txt')
    ap.add_argument('--runs', default=DEF_RUNS, help='Fungi_BFD_runs directory')
    ap.add_argument('--report', default='replaced_assemblies.tsv')
    ap.add_argument('--script', default=None,
                    help='write a reviewable shell script that moves each path into the quarantine folder')
    ap.add_argument('--new-samples', default=None,
                    help='regenerated samples.csv. When given, "gone" means: in --samples but not here '
                         '(also catches rows dropped by curation or dedup, not only rows NCBI removed)')
    ap.add_argument('--quarantine', default=None,
                    help='quarantine folder (default: <runs>/stale_<YYYYMMDD>)')
    ap.add_argument('--append-removed', action='store_true',
                    help='append rows to data/curation/removed_asmids.csv')
    ap.add_argument('--removed-csv', default='data/curation/removed_asmids.csv')
    args = ap.parse_args()

    for p in (args.samples, args.accessions):
        if not os.path.exists(p):
            ap.error(f'input not found: {p}')

    samples = read_csv(args.samples)
    sample_by_id = {r['ASMID']: r for r in samples}
    acc = read_csv(args.accessions)
    new_ids = {r['ASM_FOLDER'] for r in acc}
    acc_by_id = {r['ASM_FOLDER']: r for r in acc}
    tax_species = load_taxonomy_species(args.taxonomy)
    excluded = load_ids(args.exclude)

    new_rows = None
    if args.new_samples:
        if not os.path.exists(args.new_samples):
            ap.error(f'input not found: {args.new_samples}')
        new_rows = {r['ASMID']: r for r in read_csv(args.new_samples)}
    keep_ids = set(new_rows) if new_rows is not None else new_ids

    # index assemblies that the old samples.csv does not have yet
    candidates = [a for a in keep_ids if a not in sample_by_id]
    by_base = {}
    for a in candidates:
        by_base.setdefault(acc_base(a), []).append(a)

    def new_out(asmid):
        """Predicted `out` for a new assembly (uses the same sanitizers as create_samples_file.py)."""
        if new_rows is not None and asmid in new_rows:
            return make_sample_tag(new_rows[asmid]['SPECIES'], new_rows[asmid]['STRAIN'])
        r = acc_by_id[asmid]
        binomial = clean_species(tax_species.get(asmid) or r['SPECIES'])
        strain = backfill_strain(clean_strain(r['STRAIN']), clean_species(r['SPECIES']), binomial)
        return make_sample_tag(binomial, strain)

    def sp_strain_key(asmid):
        if new_rows is not None and asmid in new_rows:
            return (clean_species(new_rows[asmid]['SPECIES_IN']).lower(),
                    (new_rows[asmid]['STRAIN'] or '').lower())
        r = acc_by_id[asmid]
        return (clean_species(r['SPECIES']).lower(), clean_strain(r['STRAIN']).lower())

    by_sp_strain = {}
    for a in candidates:
        k = sp_strain_key(a)
        if k[1]:
            by_sp_strain.setdefault(k, []).append(a)

    # `out` values still used by rows that stay
    kept_out = {}
    for r in samples:
        if r['ASMID'] in keep_ids:
            kept_out.setdefault(make_sample_tag(r['SPECIES'], r['STRAIN']), []).append(r['ASMID'])

    gone = sorted(a for a in sample_by_id if a not in keep_ids)
    ref_name = args.new_samples or args.accessions
    print(f'[samples] {len(sample_by_id)} ASMIDs; {len(gone)} not in {ref_name}', file=sys.stderr)

    rows, rm_paths = [], []
    for a in gone:
        s = sample_by_id[a]
        old_out = make_sample_tag(s['SPECIES'], s['STRAIN'])
        rep, cls = '', ('removed_upstream' if a not in new_ids else 'dropped_by_curation_or_dedup')
        same_base = by_base.get(acc_base(a), [])
        if same_base:
            rep, cls = sorted(same_base)[-1], 'version_update'
        else:
            k = (clean_species(s['SPECIES_IN']).lower(), (s['STRAIN'] or '').lower())
            hits = by_sp_strain.get(k, []) if k[1] else []
            if hits:
                rep, cls = sorted(hits)[-1], 'replaced_by_other'
        if not rep and new_rows is not None:
            # no accession or species+strain match: still look for a new row that
            # would write into the same results folder
            same_out = [c for c in candidates if new_out(c) == old_out]
            if same_out:
                rep = sorted(same_out)[-1]
        rep_out = new_out(rep) if rep else ''

        flags = []
        if a in excluded:
            flags.append('IN_EXCLUDE_LIST')
        paths = []
        for d in PER_ASMID_DIRS:
            paths += sorted(glob.glob(os.path.join(args.runs, d, a + '.*')))
        shared = kept_out.get(old_out)
        if shared:
            flags.append('OUT_SHARED_WITH_' + '+'.join(shared))   # keep the folder
        else:
            if rep and rep_out == old_out:
                flags.append('STALE_SAME_OUT')
            for d in PER_OUT_DIRS:
                p = os.path.join(args.runs, d, old_out)
                if os.path.isdir(p):
                    paths.append(p)
        if rep and rep_out != old_out:
            flags.append('OUT_CHANGES')
        has_predict = os.path.isdir(os.path.join(args.runs, 'genome_annotation', old_out, 'predict_results'))
        rows.append({
            'ASMID': a, 'CLASS': cls, 'REPLACEMENT_ASMID': rep,
            'SPECIES': s['SPECIES'], 'STRAIN': s['STRAIN'],
            'OLD_OUT': old_out, 'NEW_OUT': rep_out,
            'HAS_PREDICT_RESULTS': int(has_predict),
            'N_PATHS': len(paths), 'FLAGS': ';'.join(flags),
        })
        rm_paths += [(a, p) for p in paths]

    fields = ['ASMID', 'CLASS', 'REPLACEMENT_ASMID', 'SPECIES', 'STRAIN', 'OLD_OUT', 'NEW_OUT',
              'HAS_PREDICT_RESULTS', 'N_PATHS', 'FLAGS']
    with open(args.report, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)

    counts = {}
    for r in rows:
        counts[r['CLASS']] = counts.get(r['CLASS'], 0) + 1
    print(f'[report] {len(rows)} rows -> {args.report}  {counts}', file=sys.stderr)
    print(f'[paths] {len(rm_paths)} paths would be cleaned', file=sys.stderr)

    if args.script:
        today = datetime.date.today()
        quarantine = args.quarantine or os.path.join(args.runs, f'stale_{today:%Y%m%d}')
        q = shlex.quote(quarantine)
        with open(args.script, 'w') as f:
            f.write('#!/bin/bash\n# Generated by find_replaced_assemblies.py on '
                    f'{today}. READ before running. Moves files; deletes nothing.\n'
                    '# Check: squeue -u $USER ; no live Nextflow run on these genomes.\n'
                    f'set -eu\nQ={q}\nmkdir -p "$Q"\n'
                    'if [ ! -e "$Q/QUARANTINE_LOG.tsv" ]; then '
                    'printf "DATE\\tASMID\\tPATH\\tSIZE_KB\\tQUARANTINE_PATH\\n" > "$Q/QUARANTINE_LOG.tsv"; fi\n')
            for a, p in rm_paths:
                rel = os.path.relpath(p, args.runs)
                dest = os.path.join(quarantine, rel)
                f.write(f'# {a}\n'
                        f'if [ -e {shlex.quote(p)} ]; then\n'
                        f'  kb=$(du -sk -- {shlex.quote(p)} | cut -f1)\n'
                        f'  mkdir -p {shlex.quote(os.path.dirname(dest))}\n'
                        f'  mv -n -- {shlex.quote(p)} {shlex.quote(dest)}\n'
                        f'  printf "%s\\t%s\\t%s\\t%s\\t%s\\n" {today} {shlex.quote(a)} '
                        f'{shlex.quote(p)} "$kb" {shlex.quote(dest)} >> "$Q/QUARANTINE_LOG.tsv"\n'
                        f'fi\n')
        os.chmod(args.script, 0o750)
        print(f'[script] {args.script} (not executed); quarantine = {quarantine}', file=sys.stderr)

    if args.append_removed:
        have = set()
        if os.path.exists(args.removed_csv):
            have = {r['ASMID'] for r in read_csv(args.removed_csv)}
        new_rows = [r for r in rows if r['ASMID'] not in have]
        with open(args.removed_csv, 'a', newline='') as f:
            w = csv.writer(f, lineterminator='\n')
            for r in new_rows:
                w.writerow([r['ASMID'], r['CLASS'], datetime.date.today().isoformat(),
                            r['SPECIES'], r['STRAIN'], r['REPLACEMENT_ASMID'],
                            r['FLAGS']])
        print(f'[removed] appended {len(new_rows)} rows to {args.removed_csv}', file=sys.stderr)


if __name__ == '__main__':
    main()
