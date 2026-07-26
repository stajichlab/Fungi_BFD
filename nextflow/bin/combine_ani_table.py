#!/usr/bin/env python3
"""
combine_ani_table.py — Merge every per-group ANI pair table + names lookup produced
by compare_ANI.nf into one labeled, queryable table (CSV + SQLite).

Inputs are two manifest files (one path per line, TAB-delimited
"<absolute_path>\\t<mtime_ms>\\t<size_bytes>", produced by compare_ANI.nf's
toManifest()) listing every currently-published:
  {group}.full.ani.tsv  or  {group}.ani.tsv   (query<TAB>reference<TAB>ANI[...])
  {group}_genome_names.tsv                    (filename, asmid, genus, species, strain)

Manifests (not a live directory scan/live channel) are the point: they're built from
a filesystem glob over every group's published output, independent of which specific
group(s) were actively recomputed in the current pipeline invocation -- this is what
lets a -resume run with only one changed group still produce a complete, all-groups
combined table instead of silently losing every other group. See
.living/learnings.md (2026-07-23) for the bug this replaced. Only the path (field 1)
is used here; mtime/size exist purely so the manifest's own content -- and therefore
this process's Nextflow input hash -- changes whenever any input file changes.

Group name is recovered from the filename by stripping the known suffixes, so this
script works unmodified across all --ani_method backends (skani/mash/sourmash write
'{group}.full.ani.tsv'; fastani's MERGE_ANI writes '{group}.ani.tsv').

Output columns (both CSV and SQLite table 'ani_pairs'):
  compare_level, taxon_group,
  query_filename, query_asmid, query_genus, query_species, query_strain,
  ref_filename,   ref_asmid,   ref_genus,   ref_species,   ref_strain,
  ani

SQLite indexes are created on (query_species), (ref_species), (query_asmid),
(ref_asmid), (taxon_group) to make ad-hoc lookups fast without a full scan.
"""

import argparse
import csv
import os
import sqlite3
import sys


ANI_SUFFIXES   = ['.full.ani.tsv', '.ani.tsv']
NAMES_SUFFIX   = '_genome_names.tsv'
NAMES_PREFIX   = 'names_'


def group_from_ani_path(path):
    base = os.path.basename(path)
    for suf in ANI_SUFFIXES:
        if base.endswith(suf):
            return base[: -len(suf)]
    return None


def group_from_names_path(path):
    """Names files arrive under two conventions depending on caller:
    REPORT_ANI's copy is '{group}_genome_names.tsv'; the raw names_map
    channel COMBINE_ANI_TABLE consumes is 'names_{group}.tsv'."""
    base = os.path.basename(path)
    if base.endswith(NAMES_SUFFIX):
        return base[: -len(NAMES_SUFFIX)]
    if base.startswith(NAMES_PREFIX) and base.endswith('.tsv'):
        return base[len(NAMES_PREFIX): -len('.tsv')]
    return None


def parse_names_tsv(path):
    """Return dict[filename] = {'asmid','genus','species','strain'}."""
    names = {}
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        cols = {name.strip().lower(): i for i, name in enumerate(header)}

        def get(parts, key):
            i = cols.get(key)
            return parts[i].strip() if i is not None and i < len(parts) else ''

        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if not parts or not parts[0]:
                continue
            fn = parts[0]
            names[fn] = {
                'asmid':   get(parts, 'asmid'),
                'genus':   get(parts, 'genus'),
                'species': get(parts, 'species'),
                'strain':  get(parts, 'strain'),
            }
    return names


def iter_ani_rows(path):
    """Yield (query_filename, ref_filename, ani_float), skipping headers/self-hits."""
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                continue
            q = os.path.basename(parts[0])
            r = os.path.basename(parts[1])
            try:
                ani = float(parts[2])
            except ValueError:
                continue  # header line (e.g. skani 'ANI') or malformed row
            if q == r:
                continue
            yield q, r, ani


def blank_names():
    return {'asmid': '', 'genus': '', 'species': '', 'strain': ''}


def read_manifest(path):
    """Read a toManifest()-produced file: one "<path>\\t<mtime>\\t<size>" record
    per line. Returns a list of (path, mtime) tuples. mtime is kept (not just
    discarded) because compare_ANI.nf's workflow unions two sources into this
    manifest -- this run's live channel items and a disk glob of everything
    already published -- and the same group can legitimately appear via both
    (see dedupe_by_group below); mtime is what lets us keep the freshest one."""
    entries = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            mtime = float(parts[1]) if len(parts) > 1 else 0.0
            entries.append((parts[0], mtime))
    return entries


def dedupe_by_group(entries, group_fn):
    """Given [(path, mtime), ...], keep only the freshest (max mtime) entry per
    group_fn(path). Groups with a null group_fn result are dropped entirely
    (same as the original glob-based filtering)."""
    best = {}
    for path, mtime in entries:
        gn = group_fn(path)
        if gn is None:
            continue
        if gn not in best or mtime > best[gn][1]:
            best[gn] = (path, mtime)
    return sorted(p for p, _m in best.values())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--ani-manifest',   required=True,
                     help='Manifest file listing *.ani.tsv / *.full.ani.tsv paths (see toManifest() in compare_ANI.nf)')
    ap.add_argument('--names-manifest', required=True,
                     help='Manifest file listing *_genome_names.tsv paths')
    ap.add_argument('--compare-level', required=True, help='Taxonomic rank used for grouping (e.g. GENUS)')
    ap.add_argument('--csv-output',    required=True)
    ap.add_argument('--db-output',     required=True)
    args = ap.parse_args()

    ani_paths = dedupe_by_group(read_manifest(args.ani_manifest), group_from_ani_path)
    names_paths = dedupe_by_group(read_manifest(args.names_manifest), group_from_names_path)

    names_by_group = {}
    for p in names_paths:
        gn = group_from_names_path(p)
        if gn is None:
            continue
        names_by_group[gn] = parse_names_tsv(p)

    if not ani_paths:
        print("Warning: no *.ani.tsv files found — writing empty outputs.", file=sys.stderr)

    header = [
        'compare_level', 'taxon_group',
        'query_filename', 'query_asmid', 'query_genus', 'query_species', 'query_strain',
        'ref_filename',   'ref_asmid',   'ref_genus',   'ref_species',   'ref_strain',
        'ani',
    ]

    col_defs = ', '.join('"%s"%s' % (c, ' REAL' if c == 'ani' else ' TEXT') for c in header)
    db = sqlite3.connect(args.db_output)
    db.execute(f"CREATE TABLE ani_pairs ({col_defs})")

    n_rows = 0
    n_groups = 0
    with open(args.csv_output, 'w', newline='') as csv_fh:
        writer = csv.writer(csv_fh)
        writer.writerow(header)

        for ani_path in ani_paths:
            group = group_from_ani_path(ani_path)
            names = names_by_group.get(group, {})
            if group not in names_by_group:
                print(f"Warning: no names file for group '{group}' "
                      f"(expected {group}{NAMES_SUFFIX}) — labels will be blank",
                      file=sys.stderr)
            n_groups += 1
            rows = []
            for q, r, ani in iter_ani_rows(ani_path):
                qi = names.get(q, blank_names())
                ri = names.get(r, blank_names())
                row = [
                    args.compare_level, group,
                    q, qi['asmid'], qi['genus'], qi['species'], qi['strain'],
                    r, ri['asmid'], ri['genus'], ri['species'], ri['strain'],
                    ani,
                ]
                writer.writerow(row)
                rows.append(row)
                n_rows += 1
            if rows:
                db.executemany(
                    f"INSERT INTO ani_pairs VALUES ({', '.join('?' for _ in header)})", rows
                )

    for col in ('query_species', 'ref_species', 'query_asmid', 'ref_asmid', 'taxon_group'):
        db.execute(f'CREATE INDEX idx_ani_pairs_{col} ON ani_pairs("{col}")')
    db.commit()
    db.close()

    print(f"Combined {n_rows} pairs across {n_groups} groups -> {args.csv_output}, {args.db_output}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
