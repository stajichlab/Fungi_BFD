"""Tests for nextflow/bin/filter_sql_for_db.py.

The filter exists because bfd_compat_views.sql, bfd_traits.sql and
bfd_extra_views.sql came from BFD_dashboard, where they only ever ran against a
complete production tables/ directory. Unguarded they make ~10 optional tables
mandatory: COMMENT ON VIEW and CREATE VIEW ... read_parquet() both raise on a
missing target, and the duckdb CLI then exits 1 with the file half applied.
"""
import subprocess
import sys
from pathlib import Path

import pytest

BIN = Path(__file__).resolve().parents[1] / "bin"
SQL_DIR = Path(__file__).resolve().parents[2] / "sql"
FILTER = BIN / "filter_sql_for_db.py"


def run(sql_text, existing, tmp_path):
    sql = tmp_path / "in.sql"
    sql.write_text(sql_text)
    names = tmp_path / "names.txt"
    names.write_text("\n".join(existing) + "\n")
    p = subprocess.run([sys.executable, str(FILTER), "--sql", str(sql),
                        "--existing-file", str(names)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return p.stdout, p.stderr


def test_comment_on_missing_view_is_dropped(tmp_path):
    out, err = run("COMMENT ON VIEW gone IS 'x';\n", ["here"], tmp_path)
    assert "COMMENT ON VIEW gone" not in out
    assert "SKIP" in err and "gone" in err


def test_comment_on_existing_view_is_kept(tmp_path):
    out, _ = run("COMMENT ON VIEW here IS 'x';\n", ["here"], tmp_path)
    assert "COMMENT ON VIEW here" in out


def test_view_over_missing_parquet_is_dropped_with_its_comment(tmp_path):
    missing = tmp_path / "nope.parquet"
    sql = (f"CREATE OR REPLACE VIEW v AS SELECT * FROM read_parquet('{missing}');\n"
           "COMMENT ON VIEW v IS 'doc';\n")
    out, err = run(sql, [], tmp_path)
    assert "CREATE OR REPLACE VIEW v" not in out
    # the COMMENT must go too: its target was never created
    assert "COMMENT ON VIEW v" not in out
    assert err.count("SKIP") == 2


def test_view_over_present_parquet_keeps_both(tmp_path):
    present = tmp_path / "yes.parquet"
    present.write_bytes(b"PAR1")          # existence is all the filter checks
    sql = (f"CREATE OR REPLACE VIEW v AS SELECT * FROM read_parquet('{present}');\n"
           "COMMENT ON VIEW v IS 'doc';\n")
    out, err = run(sql, [], tmp_path)
    assert "CREATE OR REPLACE VIEW v" in out
    assert "COMMENT ON VIEW v" in out
    assert "SKIP" not in err


def test_alias_over_missing_source_is_dropped(tmp_path):
    out, _ = run("CREATE OR REPLACE VIEW alias AS SELECT * FROM absent_src;\n", ["other"], tmp_path)
    assert "alias" not in out


def test_alias_over_present_source_is_kept(tmp_path):
    out, _ = run("CREATE OR REPLACE VIEW alias AS SELECT * FROM real_src;\n", ["real_src"], tmp_path)
    assert "CREATE OR REPLACE VIEW alias" in out


def test_semicolon_inside_a_string_literal_does_not_split(tmp_path):
    """The real COMMENT text contains semicolons and doubled quotes; a naive
    split on ';' would cut a statement in half and emit invalid SQL."""
    sql = "COMMENT ON VIEW here IS 'a; b; and ''quoted'' text';\n"
    out, _ = run(sql, ["here"], tmp_path)
    assert out.count(";") == sql.count(";")
    assert "''quoted''" in out


def test_unrecognised_statements_pass_through(tmp_path):
    """The filter guards known-optional data; it is not a validator. Anything
    it cannot classify must survive so DuckDB still reports real errors."""
    sql = "CREATE OR REPLACE MACRO in_genome(x) AS (x);\nSET threads=4;\n"
    out, _ = run(sql, [], tmp_path)
    assert "MACRO in_genome" in out
    assert "SET threads=4" in out


@pytest.mark.parametrize("name", ["bfd_compat_views.sql", "bfd_traits.sql", "bfd_extra_views.sql"])
def test_real_files_survive_intact_when_everything_exists(name, tmp_path):
    """With every name present and every path readable, the filter must be a
    no-op -- otherwise it would quietly thin a correct production build.

    The builder rewrites each file's line-1 path before filtering, so this does
    the same: the declared default points at the shared checkout, which is not
    where this working copy's data lives.
    """
    import re
    text = (SQL_DIR / name).read_text()
    first = text.splitlines()[0]
    m = re.fullmatch(r"-- (TABLES_DIR|REF_DIR)=(\S+)", first)
    assert m, first
    real = (str(Path(__file__).resolve().parents[2] / "lib" / "refs") if m.group(1) == "REF_DIR"
            else "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables")
    text = text.replace(m.group(2), real)

    names = set(re.findall(r"(?:COMMENT ON (?:VIEW|TABLE)|CREATE (?:OR REPLACE )?(?:VIEW|TABLE))\s+(\w+)", text))
    names |= set(re.findall(r"\b(?:FROM|JOIN)\s+([a-z_][\w]*)", text))
    sql = tmp_path / name
    sql.write_text(text)
    nf = tmp_path / "names.txt"
    nf.write_text("\n".join(sorted(names)) + "\n")
    p = subprocess.run([sys.executable, str(FILTER), "--sql", str(sql), "--existing-file", str(nf)],
                       capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    assert "SKIP" not in p.stderr, p.stderr
    assert p.stdout == text
