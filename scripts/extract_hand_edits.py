#!/usr/bin/env python3
"""Find hand edits in a samples.csv and record them as curated overrides.

A "hand edit" is a field in the reference samples.csv that the generator does not
reproduce from the SAME inputs. So run this with the accessions and taxonomy
snapshot that were current when the reference file was built (not today's files).
Differences caused by newer NCBI data then do not show up as hand edits.

Steps:
  1. Run create_samples_file.py with the given snapshot into a temporary file.
  2. Compare it with the reference, field by field, for ASMIDs present in both.
  3. Propose one data/curation/overrides.csv row (ASMID,FIELD,VALUE) per difference.
  4. With --write: append the rows and log them in data/curation/overrides_log.csv
     (date, ASMID, field, generated value, kept value, source, note).

It is safe to run twice. Once an override exists, the generator reproduces the
value, so the second run finds nothing new. To undo a hand edit, delete its row
from overrides.csv.

Fields SUBPHYLUM and SUBCLASS are ignored by default because older taxonomy
snapshots did not carry them (see lib/subphylum_fill.csv in NCBI_fungi).
`out` (the results folder) is built from SPECIES and STRAIN, so an override of
either one pins the existing results folder name.
"""
import argparse
import csv
import datetime
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--reference', required=True, help='samples.csv that holds the hand edits')
    ap.add_argument('--accessions', required=True, help='accessions CSV from when the reference was built')
    ap.add_argument('--taxonomy', required=True, help='taxonomy CSV from when the reference was built')
    ap.add_argument('--overrides', default='data/curation/overrides.csv')
    ap.add_argument('--log', default='data/curation/overrides_log.csv')
    ap.add_argument('--ignore-fields', default='SUBPHYLUM,SUBCLASS')
    ap.add_argument('--source', default='', help='note for the log, e.g. the reference file and date')
    ap.add_argument('--write', action='store_true')
    args = ap.parse_args()

    ignore = {x for x in args.ignore_fields.split(',') if x}
    with tempfile.TemporaryDirectory() as td:
        gen_path = os.path.join(td, 'generated.csv')
        cmd = [sys.executable, '-I', os.path.join(HERE, 'create_samples_file.py'),
               '--accessions', args.accessions, '--taxonomy', args.taxonomy,
               '--outfile', gen_path,
               # the point of this script is to regenerate from an OLD snapshot,
               # which may predate the HAPLOTYPE column
               '--allow-missing-haplotype']
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit('generator failed:\n' + r.stderr[-2000:])
        generated = {x['ASMID']: x for x in read(gen_path)}
    ref = {x['ASMID']: x for x in read(args.reference)}

    have = {}
    if os.path.exists(args.overrides):
        for row in read(args.overrides):
            have[(row['ASMID'], row['FIELD'])] = row['VALUE']

    both = [a for a in ref if a in generated]
    print(f"[hand-edits] reference {len(ref)} rows; generated {len(generated)}; compared {len(both)}; "
          f"only in reference {len(set(ref) - set(generated))}; only in generated {len(set(generated) - set(ref))}",
          file=sys.stderr)
    fields = [c for c in next(iter(ref.values())).keys() if c != 'ASMID' and c not in ignore]
    found = []
    for a in sorted(both):
        for c in fields:
            if ref[a][c] != generated[a][c] and have.get((a, c)) != ref[a][c]:
                found.append((a, c, generated[a][c], ref[a][c]))
    for a, c, g, v in found:
        print(f"[hand-edits] {a} {c}: generated={g!r} reference={v!r}", file=sys.stderr)
    print(f"[hand-edits] {len(found)} field differences not yet in {args.overrides}", file=sys.stderr)
    locus = [x for x in found if x[1] == 'LOCUSTAG']
    if locus:
        print(f"[hand-edits] WARNING {len(locus)} LOCUSTAG differences: a changed locus tag breaks "
              "existing annotations. Check before accepting.", file=sys.stderr)
    if not args.write or not found:
        if found:
            print('[hand-edits] dry run; use --write to record them', file=sys.stderr)
        return

    today = datetime.date.today().isoformat()
    with open(args.overrides, 'a', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerows([(a, c, v) for a, c, g, v in found])
    new_log = not os.path.exists(args.log)
    with open(args.log, 'a', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        if new_log:
            w.writerow(['DATE', 'ASMID', 'FIELD', 'GENERATED_VALUE', 'KEPT_VALUE', 'SOURCE'])
        w.writerows([(today, a, c, g, v, args.source) for a, c, g, v in found])
    print(f"[hand-edits] appended {len(found)} rows to {args.overrides} and {args.log}", file=sys.stderr)


if __name__ == '__main__':
    main()
