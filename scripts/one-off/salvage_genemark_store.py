#!/usr/bin/env python3
"""Copy finished GENEMARK_RUN results from a run folder's Nextflow work dirs into
<run>/genemark_store/<out>/, where FUNANNOTATE_PREDICTION reuses them instead of running
GeneMark again (nextflow/modules/funannotate/utils.nf storedGenemarkFor; Fungi_BFD_runs
DECISIONS D135).

Source tasks: GENEMARK_RUN / GENEMARK_RUN_SIB rows with status COMPLETED or CACHED in
<run>/logs/nextflow/*trace*.txt. Work dirs come from the "submitted process" lines of
<run>/.nextflow.log*. Per out, the newest GTF wins. Files are copied with their mtimes
(the reuse check needs the GTF to be newer than the masked genome), and only when the
store has no newer GTF. A task whose work dir lacks <out>.genemark.gtf or
<out>.other.gff3 is skipped.

Dry run by default; --apply copies.
Usage: salvage_genemark_store.py --run-dir <run folder> [--store DIR] [--apply]
"""
import argparse, glob, os, re, shutil, sys

ap = argparse.ArgumentParser()
ap.add_argument('--run-dir', required=True)
ap.add_argument('--store', help='default <run-dir>/genemark_store')
ap.add_argument('--apply', action='store_true')
a = ap.parse_args()
run = os.path.abspath(a.run_dir)
store = os.path.abspath(a.store or os.path.join(run, 'genemark_store'))

wd = {}
pat = re.compile(r'GENEMARK_RUN(?:_SIB)? \((.+?)\) > jobId: \d+; workDir: (\S+/(\w\w)/(\w{6})\w*)')
for log in glob.glob(f'{run}/.nextflow.log*'):
    with open(log, errors='replace') as fh:
        for line in fh:
            if 'submitted process' in line and 'GENEMARK_RUN' in line:
                m = pat.search(line)
                if m:
                    wd[f'{m.group(3)}/{m.group(4)}'] = (m.group(1), m.group(2))
print(f'work dirs from .nextflow.log*: {len(wd)}')

done = set()
for tf in glob.glob(f'{run}/logs/nextflow/*trace*.txt'):
    with open(tf) as fh:
        hdr = fh.readline().rstrip('\n').split('\t')
        for line in fh:
            r = dict(zip(hdr, line.rstrip('\n').split('\t')))
            if re.search(r':GENEMARK_RUN(_SIB)? \(', r.get('name', '')) and r.get('status') in ('COMPLETED', 'CACHED'):
                done.add(r['hash'])
print(f'GENEMARK_RUN COMPLETED/CACHED hashes in traces: {len(done)}')

best, missing_dir, missing_files = {}, 0, 0
for h in done:
    if h not in wd:
        missing_dir += 1
        continue
    out, d = wd[h]
    gtf, other = f'{d}/{out}.genemark.gtf', f'{d}/{out}.other.gff3'
    if not (os.path.exists(gtf) and os.path.exists(other)):
        missing_files += 1
        continue
    t = os.path.getmtime(gtf)
    if out not in best or t > best[out][0]:
        best[out] = (t, d)
print(f'hash not in logs: {missing_dir}; work dir without GTF/other.gff3: {missing_files}')
print(f'genomes with a salvageable result: {len(best)}')

copied = skipped_newer = empty = with_mod = 0
for out, (t, d) in sorted(best.items()):
    dest = os.path.join(store, out)
    dgtf = os.path.join(dest, f'{out}.genemark.gtf')
    if os.path.exists(dgtf) and os.path.getmtime(dgtf) >= t:
        skipped_newer += 1
        continue
    if os.path.getsize(f'{d}/{out}.genemark.gtf') == 0:
        empty += 1
    files = [f'{out}.genemark.gtf', f'{out}.other.gff3']
    if os.path.exists(f'{d}/{out}.genemark.mod'):
        files.append(f'{out}.genemark.mod')
        with_mod += 1
    copied += 1
    if a.apply:
        os.makedirs(dest, exist_ok=True)
        for f in files:
            shutil.copy2(f'{d}/{f}', os.path.join(dest, f))
verb = 'copied' if a.apply else 'would copy'
print(f'{verb}: {copied} (with .mod {with_mod}; empty GTF = GeneMark "too small" skip {empty}); '
      f'store already newer: {skipped_newer}')
if not a.apply:
    print(f'dry run; rerun with --apply to write {store}')
