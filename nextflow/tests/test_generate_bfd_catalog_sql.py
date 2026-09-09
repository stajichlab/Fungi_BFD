#!/usr/bin/env python3
"""test_generate_bfd_catalog_sql.py -- generates the catalog SQL against a
tiny fixture Parquet dataset, feeds it into an in-memory duckdb, and checks
the resulting views return the expected joined/aliased/typed rows -- and
that a missing source file is skipped instead of crashing the whole build.

Usage: python3 nextflow/tests/test_generate_bfd_catalog_sql.py
Requires: `duckdb` CLI on PATH.
"""
import csv
import io
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
BIN_DIR = TESTS_DIR.parent / "bin"
REPO_ROOT = TESTS_DIR.parent.parent
SCHEMA_PATH = str(REPO_ROOT / "sql" / "table_schema.json")


def write_parquet(path, tsv_text):
    tsv = path.with_suffix(".tsv")
    tsv.write_text(tsv_text)
    subprocess.run(
        ["duckdb", "-c",
         f"COPY (SELECT * FROM read_csv_auto('{tsv}', delim='\\t', sample_size=-1)) "
         f"TO '{path}' (FORMAT PARQUET);"],
        check=True, capture_output=True,
    )


def query_csv_rows(db_path, sql):
    result = subprocess.run(["duckdb", "-csv", str(db_path), "-c", sql],
                             capture_output=True, text=True, check=True)
    reader = csv.DictReader(io.StringIO(result.stdout))
    return list(reader)


def test_catalog_views_join_cast_and_skip_missing_correctly():
    with tempfile.TemporaryDirectory() as td:
        tables_dir = Path(td) / "tables"
        tables_dir.mkdir()
        (tables_dir / "_manifest").mkdir()

        write_parquet(tables_dir / "species.parquet",
                      "LOCUSTAG\tASMID\tSPECIES\tSTRAIN\tPHYLUM\tSUBPHYLUM\tCLASS\tSUBCLASS\tORDER\tFAMILY\tGENUS\tNCBI_TAXONID\tBUSCO_LINEAGE\n"
                      "LOC1\tGCA_1\tfoo\tbar\tAscomycota\tNA\tNA\tNA\tNA\tNA\tFusarium\t123\tfungi_odb10\n")
        write_parquet(tables_dir / "asm_stats.parquet",
                      "ASMID\tcontig_count\ttotal_length_bp\tmin_contig_bp\tmax_contig_bp\tmedian_contig_bp\tmean_contig_bp\tL50\tN50_bp\tL90\tN90_bp\tgc_pct\tn_gap_count\ttotal_n_bases\tmasked_bases\tmasked_pct\tt2t_scaffolds\ttelomere_fwd\ttelomere_rev\n"
                      "GCA_1\t10\t1000000\t500\t200000\t50000\t100000\t3\t150000\t8\t50000\t48.5\t0\t0\t10000\t1.0\t2\t2\t2\n")
        # busco_genome.parquet deliberately NOT created -- run_busco_genome
        # defaults false in production, so its view must be skipped, not
        # crash the whole catalog build.
        (tables_dir / "_manifest" / "asm_stats.json").write_text(
            '{"table": "asm_stats", "row_count": 1, "built_at": "2026-09-09T00:00:00Z", '
            '"built_by": "TEST", "merge_run_id": "r1", "schema_version": 1, '
            '"columns": ["ASMID"], "parquet_sha256": "x"}'
        )

        gen = subprocess.run(
            [sys.executable, str(BIN_DIR / "generate_bfd_catalog_sql.py"),
             "--schema", SCHEMA_PATH, "--tables-dir", str(tables_dir)],
            capture_output=True, text=True, check=True,
        )
        assert "CREATE VIEW busco_genome" not in gen.stdout, gen.stdout
        assert "SKIP: busco_genome" in gen.stdout, gen.stdout

        db_path = Path(td) / "test.duckdb"
        subprocess.run(["duckdb", str(db_path)], input=gen.stdout, text=True, check=True)

        rows = query_csv_rows(db_path, "SELECT LOCUSTAG, ASMID, TOTAL_LENGTH, GC_PERCENT FROM asm_stats;")
        assert len(rows) == 1, rows
        assert rows[0]["LOCUSTAG"] == "LOC1"
        assert rows[0]["ASMID"] == "GCA_1"
        assert rows[0]["TOTAL_LENGTH"] == "1000000"
        assert rows[0]["GC_PERCENT"] == "48.5"

        manifest_rows = query_csv_rows(db_path, 'SELECT "table", row_count FROM _table_manifest;')
        assert manifest_rows == [{"table": "asm_stats", "row_count": "1"}], manifest_rows


if __name__ == "__main__":
    test_catalog_views_join_cast_and_skip_missing_correctly()
    print("All tests passed.")
