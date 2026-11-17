#!/usr/bin/env python3
"""test_gpi_anchor_flag.py -- the predgpi join in build_BFD_duckDB.sh.

Usage: python3 nextflow/tests/test_gpi_anchor_flag.py
Requires: `duckdb` CLI on PATH (module load duckdb on HPCC).

Regression under test
---------------------
The join read

    LEFT JOIN (SELECT DISTINCT protein_id, TRUE AS is_gpi FROM predgpi) ...

with no filter on `feature`. PredGPI emits a `Chain` row for EVERY protein it
reads and a `GPI-anchor` row only for the ones it predicts an anchor on, so
every protein matched and `is_gpi_anchored` came back TRUE for all of them:
223,938,474 proteins instead of the 3,299,703 (1.5 %) that are anchored, in the
2026-11 build. Nothing errored; the column was simply a constant.

The clause is extracted from the shipped script rather than copied here, so
this test fails if the fix is reverted or edited away.
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
BUILD_SH = TESTS_DIR.parent / "bin" / "build_BFD_duckDB.sh"

# Every protein gets a Chain row; only p2 and p4 are predicted GPI-anchored.
PREDGPI_ROWS = """\
species_prefix\tprotein_id\tsource\tfeature\tscore
F0000001\tp1\tPredGPI\tChain\t1.0
F0000001\tp2\tPredGPI\tChain\t1.0
F0000001\tp2\tPredGPI\tGPI-anchor\t0.99
F0000001\tp3\tPredGPI\tChain\t1.0
F0000001\tp4\tPredGPI\tChain\t1.0
F0000001\tp4\tPredGPI\tGPI-anchor\t0.87
F0000001\tp4\tPredGPI\tGPI-anchor\t0.55
"""

PROTEIN_ROWS = """\
protein_id
p1
p2
p3
p4
"""


def fixture_parquet(path, rows_tsv):
    tsv = path.with_suffix(".tsv")
    tsv.write_text(rows_tsv)
    subprocess.run(
        ["duckdb", "-c",
         f"COPY (SELECT * FROM read_csv_auto('{tsv}', delim='\\t', sample_size=-1)) "
         f"TO '{path}' (FORMAT PARQUET);"],
        check=True, capture_output=True)


def extract_predgpi_join():
    """Pull the predgpi LEFT JOIN out of the shipped build script."""
    text = BUILD_SH.read_text()
    m = re.search(r"^LEFT JOIN \(SELECT DISTINCT protein_id, TRUE AS is_gpi.*?\) gpi USING\(protein_id\)",
                  text, re.MULTILINE | re.DOTALL)
    if not m:
        sys.exit("FAIL: could not find the predgpi LEFT JOIN in {}".format(BUILD_SH))
    return m.group(0)


def main():
    clause = extract_predgpi_join()

    failures = []

    # 1. Static guard: the clause must filter on the feature column at all.
    if "feature" not in clause:
        failures.append(
            "clause does not mention `feature`; every protein will match:\n  " + clause)

    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        predgpi = d / "predgpi.parquet"
        proteins = d / "proteins.parquet"
        fixture_parquet(predgpi, PREDGPI_ROWS)
        fixture_parquet(proteins, PROTEIN_ROWS)

        sql = (
            f"CREATE VIEW predgpi AS SELECT * FROM '{predgpi}';\n"
            f"SELECT gp.protein_id, COALESCE(gpi.is_gpi, FALSE) AS is_gpi_anchored\n"
            f"FROM (SELECT * FROM '{proteins}') gp\n"
            f"{clause}\n"
            f"ORDER BY 1;"
        )
        out = subprocess.run(["duckdb", "-noheader", "-list", "-c", sql],
                             capture_output=True, text=True)
        if out.returncode != 0:
            sys.exit("FAIL: query errored:\n{}\nSQL was:\n{}".format(out.stderr, sql))

        got = dict(line.split("|") for line in out.stdout.strip().splitlines())
        want = {"p1": "false", "p2": "true", "p3": "false", "p4": "true"}

        for pid in sorted(want):
            if got.get(pid) != want[pid]:
                failures.append(
                    "{}: expected is_gpi_anchored={}, got {}".format(
                        pid, want[pid], got.get(pid)))

        # 2. The headline property: not everything is anchored.
        n_true = sum(1 for v in got.values() if v == "true")
        if n_true == len(got):
            failures.append(
                "every protein came back anchored ({}/{}) -- this is the original bug"
                .format(n_true, len(got)))

        # 3. A protein with two GPI-anchor rows must not be duplicated.
        if len(got) != 4:
            failures.append(
                "expected 4 rows, got {} -- the DISTINCT is not deduplicating".format(len(got)))

    if failures:
        print("FAIL")
        for f in failures:
            print("  -", f)
        sys.exit(1)
    print("ok  4 proteins, 2 anchored, no duplication; clause taken from build_BFD_duckDB.sh")


if __name__ == "__main__":
    main()
