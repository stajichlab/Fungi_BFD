#!/usr/bin/env python3
"""test_bfd_common.py -- unit tests for nextflow/bin/bfd_common.py.

Usage: python3 nextflow/tests/test_bfd_common.py
"""
import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent / "bin"))

from bfd_common import assert_unique_key, gff_stem, is_gff_name, load_schema  # noqa: E402


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


def test_gff_stem_ignores_gz_suffix():
    # X.gff3.gz once became "X..gz" (str.replace), which broke the FASTA lookup
    # and failed every CALC_GENE_STATS task (2026-10-22).
    assert gff_stem("Aspergillus_fumigatus_Af293.gff3") == "Aspergillus_fumigatus_Af293"
    assert gff_stem("Aspergillus_fumigatus_Af293.gff3.gz") == "Aspergillus_fumigatus_Af293"
    assert gff_stem("/some/dir/Aspergillus_fumigatus_Af293.gff3.gz") == "Aspergillus_fumigatus_Af293"
    assert gff_stem("X.gff3.gz", ".gff3") == "X"
    assert gff_stem("X.gff.gz", "gff") == "X"


def test_gff_stem_keeps_dots_inside_the_name():
    assert gff_stem("Ascomycota_sp._Mars-9.1.gff3.gz") == "Ascomycota_sp._Mars-9.1"
    assert gff_stem("Candida_sp._CBS_1.2.gff3") == "Candida_sp._CBS_1.2"


def test_is_gff_name_accepts_plain_and_gz_only():
    assert is_gff_name("X.gff3")
    assert is_gff_name("X.gff3.gz")
    assert not is_gff_name("X.gff3.bak")
    assert not is_gff_name("X.gff3.gz.tmp")
    assert not is_gff_name("X.gff")
    assert is_gff_name("X.gff.gz", "gff")
    assert not is_gff_name("Xgff3")


if __name__ == "__main__":
    test_load_schema_reads_real_contract()
    test_assert_unique_key_passes_on_unique_rows()
    test_assert_unique_key_fails_loud_on_duplicate()
    test_assert_unique_key_fails_loud_on_missing_key_field()
    test_gff_stem_ignores_gz_suffix()
    test_gff_stem_keeps_dots_inside_the_name()
    test_is_gff_name_accepts_plain_and_gz_only()
    print("All tests passed.")
