"""The verifier's check that v_protein_annotation filters the predgpi join.

Replaces the literal `grep -q "feature = 'GPI-anchor'"` that
rebuild_pipeline_duckdb.sh ran against the builder source. That grep passed for
a database built by any other checkout, and would have failed on a correct
build if the clause were merely reformatted. This asserts on the database.
"""
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

VERIFY = Path(__file__).resolve().parents[1] / "bin" / "verify_BFD_duckDB.py"

pytestmark = pytest.mark.skipif(shutil.which("duckdb") is None,
                                reason="duckdb CLI not on PATH (module load duckdb)")


def build(tmp_path, predgpi_where):
    """A minimal database with just enough for the GPI check to have an opinion."""
    db = tmp_path / "t.duckdb"
    sql = f"""
    CREATE TABLE gene_proteins(protein_id VARCHAR, species_prefix VARCHAR, length INT);
    INSERT INTO gene_proteins VALUES ('p1','L1',100),('p2','L1',200);
    CREATE TABLE predgpi(protein_id VARCHAR, feature VARCHAR);
    INSERT INTO predgpi VALUES ('p1','Chain'),('p2','Chain'),('p2','GPI-anchor');
    CREATE VIEW v_protein_annotation AS
    SELECT gp.protein_id, COALESCE(gpi.is_gpi, FALSE) AS is_gpi_anchored
    FROM gene_proteins gp
    LEFT JOIN (SELECT DISTINCT protein_id, TRUE AS is_gpi
               FROM predgpi {predgpi_where}) gpi USING(protein_id);
    """
    p = subprocess.run(["duckdb", str(db), "-c", sql], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    return db


def run_verify(db, tmp_path):
    return subprocess.run([sys.executable, str(VERIFY), "--db", str(db),
                           "--tables-dir", str(tmp_path)],
                          capture_output=True, text=True)


def test_unfiltered_predgpi_join_is_rejected(tmp_path):
    """The actual 2026-11 regression: every protein matches, so the column is
    a constant TRUE and nothing errors."""
    db = build(tmp_path, "")
    r = run_verify(db, tmp_path)
    assert r.returncode != 0
    assert "no filter on `feature`" in r.stdout + r.stderr


def test_filtered_predgpi_join_passes_the_gpi_check(tmp_path):
    db = build(tmp_path, "WHERE feature = 'GPI-anchor'")
    r = run_verify(db, tmp_path)
    # other checks fail on this stub database; only the GPI one must not
    assert "no filter on `feature`" not in r.stdout + r.stderr
    assert "predgpi join" in r.stdout


def test_reformatted_filter_still_passes(tmp_path):
    """What the old grep got wrong: a correct build written differently.
    `grep -q "feature = 'GPI-anchor'"` fails on every one of these."""
    db = build(tmp_path, "WHERE feature IN ('GPI-anchor')")
    r = run_verify(db, tmp_path)
    assert "no filter on `feature`" not in r.stdout + r.stderr
