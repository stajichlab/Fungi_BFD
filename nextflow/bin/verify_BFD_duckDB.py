#!/usr/bin/env python3
"""Check a freshly built BFD database before it is published.

Ported from BFD_dashboard/workflows/build_bfd_db.sh (2026-11-19), which ran
these checks only for the dashboard copy. A build that succeeds and is wrong is
published by both callers -- BUILD_DUCKDB by publishDir,
rebuild_pipeline_duckdb.sh after only a size check. The labels and traits
checks (ref_* tables, v_genome_labels row counts) moved back to BFD_dashboard
on 2026-11-20 with the labels layer.

Everything runs through the duckdb CLI rather than the Python module, so this
uses the same binary that built the database and cannot hit a storage-version
mismatch (the CLI on the cluster is 1.1.3; the Python module is 1.4.4).

Exit 0 and print "verify ok", or exit 1 listing every problem found -- all
checks run even after one fails, so a single build reports everything.

Usage:
    verify_BFD_duckDB.py --db DB --tables-dir DIR [--prev-db DB] [--expect-views N]
"""
import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

# Views the build is not allowed to omit whatever else is missing. Everything
# else in BFD is optional by design: see the mkview "SKIP:" path in
# build_BFD_duckDB.sh and filter_sql_for_db.py.
CORE_VIEWS = ["species", "asm_stats", "busco_genome", "gene_info", "gene_proteins",
              "v_protein_annotation"]

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

    # G: the predgpi join in v_protein_annotation must filter on `feature`.
    #
    # PredGPI writes a `Chain` row for every protein it reads and a
    # `GPI-anchor` row only for the ones it predicts an anchor on. Without the
    # filter every protein matches the join and is_gpi_anchored is a constant
    # TRUE -- 223,938,474 proteins instead of the 3,299,703 that are anchored
    # in the 2026-11 build. Nothing errors; the column is simply wrong.
    #
    # This checks the database that was produced, not the script that produced
    # it. rebuild_pipeline_duckdb.sh used to grep the builder for the literal
    # clause, which broke twice over: it passed for a database built by some
    # other checkout, and it would fail on a correct build if the SQL were
    # merely reformatted or moved into one of the sql/*.sql files. DuckDB
    # stores the view's normalised SQL, so the clause is readable here however
    # it was written. nextflow/tests/test_gpi_anchor_flag.py covers the
    # behaviour on fixture data; this covers the shipped artefact.
    if "v_protein_annotation" in present:
        ok, vsql = duck(db, "SELECT sql FROM duckdb_views() "
                            "WHERE view_name = 'v_protein_annotation'")
        if not ok or not vsql:
            problems.append("could not read the SQL of v_protein_annotation")
        elif re.search(r"\bpredgpi\b", vsql, re.I):
            # look at what follows each mention of the table, not the whole view
            windows = [vsql[m.end():m.end() + 200]
                       for m in re.finditer(r"\bpredgpi\b", vsql, re.I)]
            if not any(re.search(r"\bfeature\b", w, re.I) for w in windows):
                problems.append(
                    "v_protein_annotation joins predgpi with no filter on `feature`: "
                    "is_gpi_anchored will be TRUE for every protein (Fungi_BFD PR #75)")
            else:
                print("  predgpi join           filtered on feature")

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
