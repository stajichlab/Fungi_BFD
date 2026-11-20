#!/usr/bin/env python3
"""Check a freshly built BFD database before it is published.

Ported from BFD_dashboard/workflows/build_bfd_db.sh (2026-11-19), which ran
these checks only for the dashboard copy. When the SQL layer moved into the
pipeline the checks did not come with it, which left two silent-wrong-build
paths open: a REF_DIR holding a truncated or empty CSV makes read_csv return
zero rows with no error, so every genome gets NULL traits; and a duplicate
LOCUSTAG in species silently multiplies v_genome_labels. Both produce a build
that succeeds and is wrong, and both callers publish it -- BUILD_DUCKDB by
publishDir, rebuild_pipeline_duckdb.sh after only a size check.

Everything runs through the duckdb CLI rather than the Python module, so this
uses the same binary that built the database and cannot hit a storage-version
mismatch (the CLI on the cluster is 1.1.3; the Python module is 1.4.4).

Exit 0 and print "verify ok", or exit 1 listing every problem found -- all
checks run even after one fails, so a single build reports everything.

Usage:
    verify_BFD_duckDB.py --db DB --tables-dir DIR [--prev-db DB] [--expect-views N]
"""
import argparse
import shutil
import subprocess
import sys
from pathlib import Path

# Views the build is not allowed to omit whatever else is missing. Everything
# else in BFD is optional by design: see the mkview "SKIP:" path in
# build_BFD_duckDB.sh and filter_sql_for_db.py.
CORE_VIEWS = ["species", "asm_stats", "busco_genome", "gene_info", "gene_proteins",
              "v_protein_annotation", "v_genome_labels", "v_genome_traits"]

REF_TABLES = ["ref_fungaltraits", "ref_growth_form_override",
              "ref_class_subphylum", "ref_phylum_subphylum"]

# Views intentionally removed; not a regression if the previous database had them.
RETIRED = {"v_gene_introns_phase"}   # dropped 2026-11-17 with the gene_introns phase fix


def duck(db, sql, read_only=True):
    """Run one statement through the CLI. Returns (ok, output)."""
    cmd = ["duckdb", "-csv", "-noheader"]
    if read_only:
        cmd.append("-readonly")
    cmd += [str(db), "-c", sql]
    p = subprocess.run(cmd, capture_output=True, text=True)
    return p.returncode == 0, (p.stdout if p.returncode == 0 else p.stderr).strip()


def scalar(db, sql):
    ok, out = duck(db, sql)
    if not ok:
        return None, out
    first = out.splitlines()[0] if out else ""
    try:
        return int(first), None
    except ValueError:
        return None, f"not a number: {first!r}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True)
    ap.add_argument("--tables-dir", required=True)
    ap.add_argument("--prev-db", default=None,
                    help="the database this one will replace; its view set must not shrink")
    ap.add_argument("--expect-views", type=int, default=None)
    args = ap.parse_args()

    if shutil.which("duckdb") is None:
        sys.exit("VERIFY FAILED: duckdb is not on PATH (module load duckdb)")
    db, tables = Path(args.db), Path(args.tables_dir)
    if not db.exists():
        sys.exit(f"VERIFY FAILED: {db} does not exist")

    problems = []

    ok, out = duck(db, "SELECT view_name FROM duckdb_views() WHERE NOT internal ORDER BY view_name")
    if not ok:
        sys.exit(f"VERIFY FAILED: cannot read the view list: {out}")
    views = [v for v in out.splitlines() if v]
    print(f"new db: {len(views)} views")
    if args.expect_views is not None and len(views) != args.expect_views:
        problems.append(f"expected {args.expect_views} views, found {len(views)}")

    ok, out = duck(db, "SELECT table_name FROM duckdb_tables() ORDER BY table_name")
    tbls = [t for t in out.splitlines() if t] if ok else []
    print(f"          {len(tbls)} tables: {', '.join(tbls) or 'none'}")

    # A: every view binds. CREATE binds eagerly, but a view created before its
    # source was replaced can still be broken at this point.
    for v in views:
        good, err = duck(db, f'SELECT * FROM "{v}" LIMIT 0')
        if not good:
            problems.append(f"view {v} does not bind: {err[:120]}")

    # B: the non-negotiable set
    present = set(views) | set(tbls)
    for name in CORE_VIEWS:
        if name not in present:
            problems.append(f"core view {name} missing")

    # C: view row count equals its parquet's, for the three schema-driven tables
    for view, pq in [("species", "species.parquet"), ("asm_stats", "asm_stats.parquet"),
                     ("busco_genome", "busco_genome.parquet")]:
        if view not in present:
            continue            # already reported by B if it is core
        pq_path = tables / pq
        if not pq_path.exists():
            problems.append(f"{view} exists but {pq_path} does not")
            continue
        nv, e1 = scalar(db, f"SELECT count(*) FROM {view}")
        npq, e2 = scalar(db, f"SELECT count(*) FROM read_parquet('{pq_path}')")
        if nv is None or npq is None:
            problems.append(f"{view}: count failed ({e1 or e2})")
            continue
        print(f"  {view:14s} view={nv:,d} parquet={npq:,d}")
        if nv != npq:
            problems.append(f"{view}: view {nv} != parquet {npq}")

    # D: the trait tables exist AND hold rows. read_csv on a truncated or empty
    # file returns zero rows without raising, so "exists" is not enough.
    for name in REF_TABLES:
        n, err = scalar(db, f"SELECT count(*) FROM {name}")
        if n is None:
            problems.append(f"table {name} missing or unreadable: {(err or '')[:120]}")
            continue
        print(f"  {name:26s} {n:,d} rows")
        if n == 0:
            problems.append(f"table {name} is empty (check REF_DIR)")

    # E: the labels join must not multiply or lose genomes
    if {"v_genome_labels", "species"} <= present:
        n_sp, _ = scalar(db, "SELECT count(*) FROM species")
        n_lab, _ = scalar(db, "SELECT count(*) FROM v_genome_labels")
        n_key, _ = scalar(db, "SELECT count(DISTINCT LOCUSTAG) FROM v_genome_labels")
        n_nosub, _ = scalar(db, "SELECT count(*) FROM v_genome_labels WHERE SUBPHYLUM IS NULL")
        if None in (n_sp, n_lab, n_key):
            problems.append("v_genome_labels check failed to run")
        else:
            print(f"  v_genome_labels rows={n_lab:,d} distinct LOCUSTAG={n_key:,d} "
                  f"species={n_sp:,d}; SUBPHYLUM NULL for {n_nosub if n_nosub is not None else '?'}")
            if not (n_lab == n_key == n_sp):
                problems.append(f"v_genome_labels rows {n_lab}, distinct LOCUSTAG {n_key}, "
                                f"species rows {n_sp} must be equal")

    # F: the view set must not shrink against the database being replaced
    if args.prev_db and Path(args.prev_db).exists():
        ok, out = duck(args.prev_db,
                       "SELECT view_name FROM duckdb_views() WHERE NOT internal ORDER BY view_name")
        if not ok:
            print(f"previous db: unreadable, skipping the shrink check ({out[:100]})")
        else:
            old = {v for v in out.splitlines() if v}
            lost = sorted(old - set(views) - RETIRED)
            print(f"previous db: {len(old)} views; views lost: {lost or 'none'}")
            if lost:
                problems.append(f"views present before but not now: {lost}")
    else:
        print("previous db: not given, shrink check skipped")

    if problems:
        sys.exit("VERIFY FAILED:\n  " + "\n  ".join(problems))
    print("verify ok")


if __name__ == "__main__":
    main()
