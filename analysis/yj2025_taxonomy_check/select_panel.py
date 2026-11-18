#!/usr/bin/env python3
"""Pick the reference proteomes for the YJ-2025a/b placement check.

Reads BFD species.parquet + asm_stats.parquet + busco_genome.parquet (read only) and writes a
TSV of reference genomes with the path of each annotated proteome:

  * one genome per GENUS (best BUSCO complete, then N50) for every Pucciniomycotina class
    (Microbotryomycetes is the class the V18 tree points to; the other classes are the test of it),
  * one genome per SPECIES for Rhodotorula (the genus of the V18 sister, R. pacifica),
  * a few controls: Saccharomycotina (the NCBI label), Agaricomycotina, Ustilaginomycotina.

Queries (the two YJ-2025 genomes and R. pacifica as a calibration query) are listed with role=query.

Usage: select_panel.py --runs /bigdata/.../Fungi_BFD_runs --out panel.tsv
"""
import argparse
import os
import sys

import duckdb

PUCC_CLASSES = ['Microbotryomycetes', 'Cystobasidiomycetes', 'Agaricostilbomycetes', 'Pucciniomycetes',
                'Mixiomycetes', 'Tritirachiomycetes', 'Atractiellomycetes', 'Classiculomycetes',
                'Spiculogloeomycetes', 'Cryptomycocolacomycetes']
CONTROL_CLASSES = {  # class -> how many genera to take (largest BUSCO first)
    'Saccharomycetes': 2, 'Pichiomycetes': 2, 'Tremellomycetes': 2, 'Agaricomycetes': 2,
    'Ustilaginomycetes': 2, 'Exobasidiomycetes': 1}
QUERIES = ['FF90B97F', 'F6EB82EF', 'FBBC614A']  # YJ-2025a, YJ-2025b, Rhodotorula pacifica JCM10908


def pep_path(runs, species, strain, locustag):
    """input/pep/<SPECIES>_<STRAIN>.proteins.fa.gz with spaces -> '_'; verified by the LOCUSTAG of the first id."""
    import gzip
    import glob
    base = os.path.join(runs, 'input', 'pep')
    cands = []
    if strain:
        cands.append(os.path.join(base, f'{species} {strain}'.replace(' ', '_') + '.proteins.fa.gz'))
    cands += sorted(glob.glob(os.path.join(base, species.replace(' ', '_') + '*.proteins.fa.gz')))
    for c in cands:
        if not os.path.exists(c):
            continue
        with gzip.open(c, 'rt') as fh:
            first = fh.readline()
        if first.startswith('>' + locustag + '_'):
            return c
    return ''


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--runs', default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    T = os.path.join(a.runs, 'tables')
    con = duckdb.connect()
    con.execute(f"""
      CREATE TABLE g AS
      SELECT s.LOCUSTAG, s.PHYLUM, s.CLASS, s."ORDER" AS ORDER_, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
             b.complete_pct AS busco, a.N50_bp
      FROM '{T}/species.parquet' s
      JOIN '{T}/asm_stats.parquet' a ON a.ASMID = s.ASMID
      LEFT JOIN '{T}/busco_genome.parquet' b ON b.ASMID = s.ASMID
      WHERE a.annotation_status = 'annotated'""")
    rank = "row_number() OVER (PARTITION BY {k} ORDER BY busco DESC NULLS LAST, N50_bp DESC NULLS LAST, LOCUSTAG)"
    rows = con.execute(f"""
      WITH pucc AS (SELECT *, 'pucciniomycotina' AS role, {rank.format(k='GENUS')} AS r FROM g
                    WHERE CLASS IN ({','.join("'%s'" % c for c in PUCC_CLASSES)}) AND GENUS IS NOT NULL AND GENUS <> ''),
           rhodo AS (SELECT *, 'rhodotorula' AS role, {rank.format(k='SPECIES')} AS r FROM g WHERE GENUS = 'Rhodotorula'),
           ctl_g AS (SELECT *, {rank.format(k='GENUS')} AS rg FROM g WHERE CLASS IN ({','.join("'%s'" % c for c in CONTROL_CLASSES)})),
           ctl AS (SELECT * EXCLUDE (rg), 'control' AS role,
                          row_number() OVER (PARTITION BY CLASS ORDER BY busco DESC NULLS LAST, LOCUSTAG) AS r
                   FROM ctl_g WHERE rg = 1),
           q AS (SELECT *, 'query' AS role, 1 AS r FROM g WHERE LOCUSTAG IN ({','.join("'%s'" % x for x in QUERIES)}))
      SELECT LOCUSTAG, role, PHYLUM, CLASS, ORDER_, FAMILY, GENUS, SPECIES, STRAIN, busco FROM q
      UNION ALL SELECT LOCUSTAG, role, PHYLUM, CLASS, ORDER_, FAMILY, GENUS, SPECIES, STRAIN, busco FROM pucc WHERE r = 1
      UNION ALL SELECT LOCUSTAG, role, PHYLUM, CLASS, ORDER_, FAMILY, GENUS, SPECIES, STRAIN, busco FROM rhodo WHERE r = 1
      UNION ALL SELECT LOCUSTAG, role, PHYLUM, CLASS, ORDER_, FAMILY, GENUS, SPECIES, STRAIN, busco FROM ctl
                WHERE r <= {max(CONTROL_CLASSES.values())}
      """).fetchall()
    # keep the per-class control cap and drop duplicates (a Rhodotorula genus rep is also in pucc)
    seen, out, ctl_n = set(), [], {}
    for r in rows:
        loc, role, cls = r[0], r[1], r[3]
        if role == 'control':
            ctl_n[cls] = ctl_n.get(cls, 0) + 1
            if ctl_n[cls] > CONTROL_CLASSES[cls]:
                continue
        if loc in seen:
            continue
        seen.add(loc)
        out.append(r)
    missing = [q for q in QUERIES if q not in seen]
    if missing:
        sys.exit(f'query genomes not found or not annotated: {missing}')
    n_nopep = 0
    with open(a.out, 'w') as o:
        o.write('LOCUSTAG\trole\tPHYLUM\tCLASS\tORDER\tFAMILY\tGENUS\tSPECIES\tSTRAIN\tbusco_complete\tpep\n')
        for r in out:
            p = pep_path(a.runs, r[7] or '', r[8] or '', r[0])
            if not p:
                n_nopep += 1
                print(f'no proteome file for {r[0]} {r[7]} {r[8]}', file=sys.stderr)
                if r[1] == 'query':
                    sys.exit('a query proteome is missing')
                continue
            o.write('\t'.join('' if v is None else str(v) for v in r) + '\t' + p + '\n')
    print(f'{len(out)} genomes selected, {n_nopep} without a proteome file, written {a.out}', file=sys.stderr)


if __name__ == '__main__':
    main()
