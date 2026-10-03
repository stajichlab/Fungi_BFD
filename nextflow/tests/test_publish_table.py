#!/usr/bin/env python3
"""test_publish_table.py -- unit tests for nextflow/bin/publish_table.py.

Usage: python3 nextflow/tests/test_publish_table.py
Requires: `duckdb` CLI on PATH (module load duckdb on HPCC).
"""
import csv
import io
import json
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
BIN_DIR = TESTS_DIR.parent / "bin"
REPO_ROOT = TESTS_DIR.parent.parent
SCHEMA_PATH = str(REPO_ROOT / "sql" / "table_schema.json")


def make_fixture_parquet(path, rows_tsv):
    tsv = path.with_suffix(".tsv")
    tsv.write_text(rows_tsv)
    subprocess.run(
        ["duckdb", "-c",
         f"COPY (SELECT * FROM read_csv_auto('{tsv}', delim='\\t', sample_size=-1)) "
         f"TO '{path}' (FORMAT PARQUET);"],
        check=True, capture_output=True,
    )


def run_publish(local, tables_dir, run_id, extra_args=()):
    return subprocess.run(
        [sys.executable, str(BIN_DIR / "publish_table.py"),
         "--table", "busco_genome",
         "--local-parquet", str(local),
         "--tables-dir", str(tables_dir),
         "--built-by", "TEST",
         "--merge-run-id", run_id,
         "--schema", SCHEMA_PATH,
         *extra_args],
        capture_output=True, text=True,
    )


def test_publish_writes_parquet_and_manifest_atomically():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        local = td / "busco_genome.parquet"
        make_fixture_parquet(local, "ASMID\tcomplete_pct\nGCA_1\t99.0\n")
        tables_dir = td / "tables"
        tables_dir.mkdir()

        result = run_publish(local, tables_dir, "test-run-1")
        assert result.returncode == 0, result.stderr

        published = tables_dir / "busco_genome.parquet"
        manifest_file = tables_dir / "_manifest" / "busco_genome.json"
        assert published.is_file()
        assert manifest_file.is_file()
        manifest = json.loads(manifest_file.read_text())
        assert manifest["table"] == "busco_genome"
        assert manifest["row_count"] == 1
        assert manifest["built_by"] == "TEST"
        assert manifest["merge_run_id"] == "test-run-1"
        assert manifest["schema_version"] == 1
        assert "ASMID" in manifest["columns"]  # derived from the real Parquet file via DESCRIBE
        assert len(manifest["parquet_sha256"]) == 64

        leftovers = list(tables_dir.glob(".tmp.*")) + list((tables_dir / "_manifest").glob(".tmp.*"))
        assert leftovers == [], leftovers


def test_publish_rejects_empty_parquet():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        local = td / "busco_genome.parquet"
        make_fixture_parquet(local, "ASMID\tcomplete_pct\n")  # header only, 0 rows
        tables_dir = td / "tables"
        tables_dir.mkdir()

        result = run_publish(local, tables_dir, "test-run-2")
        assert result.returncode != 0
        assert "0 rows" in result.stderr or "0 rows" in result.stdout


def test_publish_shrink_check_is_opt_in_and_overridable():
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        tables_dir = td / "tables"
        tables_dir.mkdir()

        big = td / "big.parquet"
        make_fixture_parquet(big, "ASMID\tcomplete_pct\nA1\t1.0\nA2\t2.0\n")
        small = td / "small.parquet"
        make_fixture_parquet(small, "ASMID\tcomplete_pct\nA1\t1.0\n")

        assert run_publish(big, tables_dir, "run-a", ["--enforce-no-shrink"]).returncode == 0

        # Without --enforce-no-shrink (the default, matching merge_all=false runs),
        # a smaller table is accepted -- this is normal behavior for a
        # current-run-only merge, not a data-integrity violation.
        result_default = run_publish(small, tables_dir, "run-b")
        assert result_default.returncode == 0, result_default.stderr

        # Re-publish the big one, then confirm --enforce-no-shrink actually
        # rejects a real shrink when asked to.
        assert run_publish(big, tables_dir, "run-c", ["--enforce-no-shrink"]).returncode == 0
        result_shrink = run_publish(small, tables_dir, "run-d", ["--enforce-no-shrink"])
        assert result_shrink.returncode != 0

        result_override = run_publish(
            small, tables_dir, "run-e", ["--enforce-no-shrink", "--allow-shrink"]
        )
        assert result_override.returncode == 0, result_override.stderr


if __name__ == "__main__":
    test_publish_writes_parquet_and_manifest_atomically()
    test_publish_rejects_empty_parquet()
    test_publish_shrink_check_is_opt_in_and_overridable()
    print("All tests passed.")
