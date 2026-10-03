#!/usr/bin/env python3
"""test_bfd_common.py -- unit tests for nextflow/bin/bfd_common.py.

Usage: python3 nextflow/tests/test_bfd_common.py
"""
import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent / "bin"))

from bfd_common import assert_unique_key, load_schema  # noqa: E402


def test_load_schema_reads_real_contract():
    repo_root = TESTS_DIR.parent.parent
    schema = load_schema(str(repo_root / "sql" / "table_schema.json"))
    assert set(schema) == {"species", "asm_stats", "busco_genome"}, schema.keys()
    assert schema["asm_stats"]["key"] == "ASMID"
    assert schema["asm_stats"]["join"] == {"table": "species", "on": "ASMID", "how": "left"}


def test_assert_unique_key_passes_on_unique_rows():
    rows = [{"ASMID": "A1", "x": 1}, {"ASMID": "A2", "x": 2}]
    assert_unique_key(rows, "ASMID", "test_table")  # must not raise/exit


def _run_snippet(code):
    bin_dir = str(TESTS_DIR.parent / "bin")
    return subprocess.run(
        [sys.executable, "-c", f"import sys; sys.path.insert(0, {bin_dir!r}); {code}"],
        capture_output=True, text=True,
    )


def test_assert_unique_key_fails_loud_on_duplicate():
    result = _run_snippet(
        "from bfd_common import assert_unique_key; "
        "assert_unique_key([{'ASMID':'A1'},{'ASMID':'A1'}], 'ASMID', 'test_table')"
    )
    assert result.returncode != 0, "expected non-zero exit on duplicate key"
    assert "duplicate" in result.stderr.lower(), result.stderr


def test_assert_unique_key_fails_loud_on_missing_key_field():
    result = _run_snippet(
        "from bfd_common import assert_unique_key; "
        "assert_unique_key([{'OTHER':'x'}], 'ASMID', 'test_table')"
    )
    assert result.returncode != 0, "expected non-zero exit, not an uncaught KeyError"
    assert "missing key field" in result.stderr.lower(), result.stderr
    assert "Traceback" not in result.stderr, "should be a clean sys.exit, not a raw traceback"


if __name__ == "__main__":
    test_load_schema_reads_real_contract()
    test_assert_unique_key_passes_on_unique_rows()
    test_assert_unique_key_fails_loud_on_duplicate()
    test_assert_unique_key_fails_loud_on_missing_key_field()
    print("All tests passed.")
