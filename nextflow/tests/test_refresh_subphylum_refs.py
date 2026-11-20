"""Tests for scripts/refresh_subphylum_refs.sh.

Two layers. The text tests pin the statement slice against the real
sql/bfd_traits.sql and always run. The end-to-end tests build a miniature
database and run the script against it; they need the duckdb CLI and skip
without it.
"""
import os
import re
import shutil
import subprocess
import textwrap
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "refresh_subphylum_refs.sh"
TRAITS_SQL = REPO / "sql" / "bfd_traits.sql"
REFS = REPO / "lib" / "refs"

SLICE_START = "/^CREATE OR REPLACE TABLE ref_class_subphylum/"
SLICE_END = "/^COMMENT ON TABLE ref_phylum_subphylum/"

# DuckDB >= 1.4 enables an interactive "agent mode" when AI_AGENT is set and
# stdout is a pipe, which would corrupt the values the script parses.
ENV = {k: v for k, v in os.environ.items() if k != "AI_AGENT"}
HAVE_DUCKDB = shutil.which("duckdb") is not None
needs_duckdb = pytest.mark.skipif(not HAVE_DUCKDB, reason="duckdb CLI not available")


def _slice(text):
    """The statements the script keeps, applied to arbitrary SQL text."""
    out, keeping = [], False
    for line in text.splitlines(keepends=True):
        if line.startswith("CREATE OR REPLACE TABLE ref_class_subphylum"):
            keeping = True
        if keeping:
            out.append(line)
        if line.startswith("COMMENT ON TABLE ref_phylum_subphylum"):
            keeping = False
    return "".join(out)


# --------------------------------------------------------------------------
# the slice, against the real SQL
# --------------------------------------------------------------------------

def test_script_is_executable():
    assert SCRIPT.exists()
    assert os.access(SCRIPT, os.X_OK)


def test_slice_takes_exactly_the_two_subphylum_tables():
    got = _slice(TRAITS_SQL.read_text(encoding="utf-8"))
    creates = re.findall(r"^CREATE OR REPLACE TABLE (\w+)", got, re.M)
    assert creates == ["ref_class_subphylum", "ref_phylum_subphylum"]


def test_slice_excludes_the_other_ref_tables():
    """The FungalTraits tables are large and unrelated; refreshing them here
    would make a 4-second operation a slow one and change traits unasked."""
    got = _slice(TRAITS_SQL.read_text(encoding="utf-8"))
    assert "ref_fungaltraits" not in got
    assert "ref_growth_form_override" not in got


def test_slice_excludes_view_definitions():
    """v_genome_traits must not be recreated by this path -- see the script's
    'WHAT THIS DOES NOT DO'."""
    got = _slice(TRAITS_SQL.read_text(encoding="utf-8"))
    assert "CREATE OR REPLACE VIEW" not in got


def test_slice_reads_both_tsvs():
    got = _slice(TRAITS_SQL.read_text(encoding="utf-8"))
    reads = re.findall(r"read_csv\('([^']+)'", got)
    assert sorted(Path(p).name for p in reads) == [
        "class_subphylum.tsv",
        "phylum_subphylum.tsv",
    ]


def test_slice_is_empty_if_the_table_is_renamed():
    """Guards the guard: if bfd_traits.sql renames the first table, the slice
    collapses and the script's CREATE-count check must be what catches it."""
    mangled = TRAITS_SQL.read_text(encoding="utf-8").replace(
        "CREATE OR REPLACE TABLE ref_class_subphylum",
        "CREATE OR REPLACE TABLE ref_class_to_subphylum",
    )
    assert _slice(mangled).strip() == ""


def test_ref_dir_header_is_substitutable():
    first = TRAITS_SQL.read_text(encoding="utf-8").splitlines()[0]
    assert first.startswith("-- REF_DIR=")
    default = first[len("-- REF_DIR="):]
    assert default
    # every read_csv in the slice must sit under the declared prefix, or the
    # substitution would leave a stale absolute path behind
    got = _slice(TRAITS_SQL.read_text(encoding="utf-8"))
    for p in re.findall(r"read_csv\('([^']+)'", got):
        assert p.startswith(default), p


# --------------------------------------------------------------------------
# end to end, against a miniature database
# --------------------------------------------------------------------------

CLASS_TSV_SPARSE = "CLASS\tSUBPHYLUM\tsource\tnote\nAgaricomycetes\tAgaricomycotina\ttaxonomy_py\t\n"
CLASS_TSV_FULL = CLASS_TSV_SPARSE + "Sordariomycetes\tPezizomycotina\ttaxonomy_py\t\n"
PHYLUM_TSV_SPARSE = "PHYLUM\tSUBPHYLUM\tsource\tnote\nGlomeromycota\tGlomeromycotina\tproposed\t\n"
PHYLUM_TSV_FULL = PHYLUM_TSV_SPARSE + "Mucoromycota\tMucoromycotina\tproposed\t\n"


def _duckdb(db, sql):
    return subprocess.run(
        ["duckdb", str(db), "-c", sql],
        capture_output=True, text=True, env=ENV, check=True,
    ).stdout


def _make_db(tmp_path, with_view=True, view_reads_phylum=True):
    """A miniature stand-in: species plus the two lookup tables, and a
    v_genome_traits shaped like the real one's SUBPHYLUM arm."""
    db = tmp_path / "mini.duckdb"
    phylum_arm = (
        "LEFT JOIN ref_phylum_subphylum ph ON s.PHYLUM = ph.PHYLUM"
        if view_reads_phylum else ""
    )
    coalesce = (
        "coalesce(s.SUBPHYLUM, c.SUBPHYLUM, ph.SUBPHYLUM)"
        if view_reads_phylum else "coalesce(s.SUBPHYLUM, c.SUBPHYLUM)"
    )
    sql = textwrap.dedent(f"""
        CREATE TABLE species (LOCUSTAG VARCHAR, PHYLUM VARCHAR, CLASS VARCHAR, SUBPHYLUM VARCHAR);
        INSERT INTO species VALUES
            ('A', 'Basidiomycota', 'Agaricomycetes',  NULL),   -- via class table
            ('B', 'Ascomycota',    'Sordariomycetes', NULL),   -- via class table, only once full
            ('C', 'Mucoromycota',  NULL,              NULL),   -- via phylum table, only once full
            ('D', 'Ascomycota',    'Eurotiomycetes',  'Pezizomycotina'),  -- already filled
            ('E', 'Chytridiomycota', 'Nowhere',       NULL);   -- stays NULL
        CREATE TABLE ref_class_subphylum  (CLASS VARCHAR, SUBPHYLUM VARCHAR, source VARCHAR, note VARCHAR);
        CREATE TABLE ref_phylum_subphylum (PHYLUM VARCHAR, SUBPHYLUM VARCHAR, source VARCHAR, note VARCHAR);
        INSERT INTO ref_class_subphylum  VALUES ('Agaricomycetes', 'Agaricomycotina', 'taxonomy_py', NULL);
        INSERT INTO ref_phylum_subphylum VALUES ('Glomeromycota', 'Glomeromycotina', 'proposed', NULL);
    """)
    if with_view:
        sql += textwrap.dedent(f"""
            CREATE VIEW v_genome_traits AS
            SELECT s.LOCUSTAG, {coalesce} AS SUBPHYLUM
            FROM species s
            LEFT JOIN ref_class_subphylum c ON s.CLASS = c.CLASS
            {phylum_arm};
        """)
    _duckdb(db, sql)
    return db


def _refdir(tmp_path, class_tsv, phylum_tsv):
    d = tmp_path / "refs"
    d.mkdir(exist_ok=True)
    (d / "class_subphylum.tsv").write_text(class_tsv, encoding="utf-8")
    (d / "phylum_subphylum.tsv").write_text(phylum_tsv, encoding="utf-8")
    return d


def _run(db, refdir, **extra):
    env = dict(ENV, PIPELINE=str(REPO), REF_DIR=str(refdir), **extra)
    return subprocess.run([str(SCRIPT), str(db)],
                          capture_output=True, text=True, env=env)


def _nulls(db):
    # -csv -noheader, not the default box output: the box prints the column
    # TYPE ("int64") above the value, and a bare \d+ search matches the 64.
    out = subprocess.run(
        ["duckdb", "-csv", "-noheader", str(db), "-c",
         "SELECT count(*) FROM v_genome_traits WHERE SUBPHYLUM IS NULL;"],
        capture_output=True, text=True, env=ENV, check=True,
    ).stdout
    return int(out.strip())


@needs_duckdb
def test_refresh_fills_subphylum_from_both_tables(tmp_path):
    db = _make_db(tmp_path)
    assert _nulls(db) == 3          # B, C and E start NULL
    r = _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    assert r.returncode == 0, r.stderr
    assert _nulls(db) == 1          # only E, which no table covers
    assert "2 genome(s) gained a SUBPHYLUM" in r.stdout


@needs_duckdb
def test_refresh_reports_before_and_after(tmp_path):
    db = _make_db(tmp_path)
    r = _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    assert "before:    ref_class=1 ref_phylum=1 subphylum_null=3" in r.stdout
    assert "after:     ref_class=2 ref_phylum=2 subphylum_null=1" in r.stdout


@needs_duckdb
def test_refresh_leaves_species_subphylum_alone(tmp_path):
    """species.SUBPHYLUM is the pipeline's column; this script must not edit it."""
    db = _make_db(tmp_path)
    before = _duckdb(db, "SELECT SUBPHYLUM FROM species WHERE LOCUSTAG = 'D';")
    _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    assert _duckdb(db, "SELECT SUBPHYLUM FROM species WHERE LOCUSTAG = 'D';") == before


@needs_duckdb
def test_shrinking_lookup_tables_is_refused(tmp_path):
    db = _make_db(tmp_path)
    _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    r = _run(db, _refdir(tmp_path, CLASS_TSV_SPARSE, PHYLUM_TSV_SPARSE))
    assert r.returncode == 1
    assert "lookup rows dropped" in r.stderr


@needs_duckdb
def test_shrink_can_be_overridden(tmp_path):
    db = _make_db(tmp_path)
    _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    r = _run(db, _refdir(tmp_path, CLASS_TSV_SPARSE, PHYLUM_TSV_SPARSE), ALLOW_SHRINK="1")
    assert r.returncode == 0, r.stderr
    assert _nulls(db) == 3


@needs_duckdb
def test_missing_view_is_refused(tmp_path):
    db = _make_db(tmp_path, with_view=False)
    r = _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    assert r.returncode == 1
    assert "no v_genome_traits view" in r.stderr


@needs_duckdb
def test_view_that_dropped_a_lookup_table_is_refused(tmp_path):
    """If the SQL moves on from the two-table coalesce, this shortcut is no
    longer equivalent to a build and must say so instead of half-working."""
    db = _make_db(tmp_path, view_reads_phylum=False)
    r = _run(db, _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    assert r.returncode == 1
    assert "no longer reads ref_phylum_subphylum" in r.stderr
    assert "full build" in r.stderr


@needs_duckdb
def test_missing_ref_file_is_refused(tmp_path):
    db = _make_db(tmp_path)
    d = _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL)
    (d / "phylum_subphylum.tsv").unlink()
    r = _run(db, d)
    assert r.returncode == 1
    assert "phylum_subphylum.tsv" in r.stderr


@needs_duckdb
def test_missing_database_is_refused(tmp_path):
    r = _run(tmp_path / "nope.duckdb", _refdir(tmp_path, CLASS_TSV_FULL, PHYLUM_TSV_FULL))
    assert r.returncode == 1
    assert "no database at" in r.stderr


def test_usage_without_an_argument():
    r = subprocess.run([str(SCRIPT)], capture_output=True, text=True,
                       env=dict(ENV, PIPELINE=str(REPO)))
    assert r.returncode == 2
    assert "usage:" in r.stderr
