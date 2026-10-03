#!/usr/bin/env python3
"""Reset species whose RNA-seq cache is empty although SRA lists usable single-end runs.

Before Fungi_BFD fa8c7e9 (2026-07-16, single-end routing), species with only SINGLE-layout
SRA runs got 0-byte reads files, and RNASEQ_PREPARE then wrote a 0-byte Trinity FASTA.
SRA_FETCH/SRA_FETCH_SE and RNASEQ_PREPARE cache through storeDir, so those empty files
stop every later run from fetching the reads again. Fungi_BFD_runs DECISIONS D137.

A species is selected when all hold:
  - <rnaseq_reads>/sra_query/<species>.sra_query.csv lists >= 1 SINGLE run whose accession
    is not tagged 'skip' in rnaseq_blacklist.csv;
  - it is not in rnaseq_force_no_rnaseq.csv;
  - none of <species>_norm_{R1,R2,SE}.fastq.gz has data (each is absent or 0 bytes).

For each selected species, --apply removes only 0-byte files:
  <rnaseq_reads>/<species>_norm_{R1,R2,SE}.fastq.gz and <rnaseq_data>/<species>.trinity-GG.fasta.
Non-empty files are never touched. The next run then fetches the SE reads (needs
enable_single_end = true) and rebuilds Trinity.

Also writes --samples-out: the rows of --samples for the selected species (minus any ASMID
in --exclude-samples), for a follow-up annotation run.

Dry run by default. Usage:
  reset_empty_se_rnaseq.py --runs-dir /bigdata/.../Fungi_BFD_runs \
      [--samples Fungi_BFD/samples.csv] [--exclude-samples do_annotation_wave1/samples_wave1.csv] \
      [--samples-out FILE] [--apply]
"""
import argparse, csv, datetime, glob, os

ap = argparse.ArgumentParser()
ap.add_argument('--runs-dir', required=True, help='folder with rnaseq_reads/, rnaseq_data/ and the CSV lists')
ap.add_argument('--samples', help='samples CSV to draw follow-up rows from')
ap.add_argument('--exclude-samples', nargs='*', default=[], help='samples CSVs whose ASMIDs are already covered')
ap.add_argument('--samples-out', help='default <runs-dir>/samples_reset_empty_se.csv')
ap.add_argument('--apply', action='store_true')
a = ap.parse_args()
B = os.path.abspath(a.runs_dir)
R, D = f'{B}/rnaseq_reads', f'{B}/rnaseq_data'

skip = set()
for row in csv.reader(open(f'{B}/rnaseq_blacklist.csv')):
    if len(row) >= 4 and not row[0].startswith('#') and row[3].strip() == 'skip':
        skip.add(row[0].strip())
force = set()
for row in csv.reader(open(f'{B}/rnaseq_force_no_rnaseq.csv')):
    if row and not row[0].startswith('#'):
        force.add(row[0].strip())


def size(p):
    return os.path.getsize(p) if os.path.exists(p) else None


selected, to_remove = [], []
for q in sorted(glob.glob(f'{R}/sra_query/*.sra_query.csv')):
    sp = os.path.basename(q)[:-len('.sra_query.csv')]
    if sp in force:
        continue
    runs = [r['sra_accession'].strip() for r in csv.DictReader(open(q))
            if (r.get('layout') or '').strip() == 'SINGLE' and r['sra_accession'].strip() not in skip]
    if not runs:
        continue
    reads = [f'{R}/{sp}_norm_{x}.fastq.gz' for x in ('R1', 'R2', 'SE')]
    if any((size(p) or 0) > 0 for p in reads):
        continue
    empties = [p for p in reads + [f'{D}/{sp}.trinity-GG.fasta'] if size(p) == 0]
    mt = [os.path.getmtime(p) for p in reads if os.path.exists(p)]
    selected.append((sp, runs, datetime.date.fromtimestamp(max(mt)).isoformat() if mt else 'absent', empties))
    to_remove += empties

print(f'species selected: {len(selected)}; 0-byte files to remove: {len(to_remove)}')
for sp, runs, when, empties in selected:
    print(f'  {sp}\treads emptied {when}\tSE runs {";".join(runs)}\t{len(empties)} empty files')

if a.samples:
    covered = set()
    for f in a.exclude_samples:
        covered |= {r['ASMID'].strip() for r in csv.DictReader(open(f))}
    names = {sp for sp, *_ in selected}
    with open(a.samples) as fh:
        rd = csv.DictReader(fh)
        rows = [r for r in rd if '_'.join(r['SPECIES'].split()) in names and r['ASMID'].strip() not in covered]
        fields = rd.fieldnames
    out = a.samples_out or f'{B}/samples_reset_empty_se.csv'
    print(f'follow-up genomes (not in excluded CSVs): {len(rows)} -> {out if a.apply else "(dry run, not written)"}')
    if a.apply:
        with open(out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader(); w.writerows(rows)

if a.apply:
    log = f'{B}/reset_empty_se_rnaseq.{datetime.datetime.now():%Y%m%d_%H%M%S}.log'
    with open(log, 'w') as fh:
        for p in to_remove:
            if size(p) == 0:            # re-check right before removing
                os.remove(p)
                fh.write(p + '\n')
    print(f'removed {len(to_remove)} files; list in {log}')
else:
    print('dry run; rerun with --apply to remove the files')
