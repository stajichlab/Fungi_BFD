#!/usr/bin/env python3
"""test_merge_deeploc.py -- unit tests for nextflow/bin/merge_deeploc.py.

Usage: python3 nextflow/tests/test_merge_deeploc.py
Requires: `duckdb` CLI on PATH (module load duckdb on HPCC).

The fixture rows are real DeepLoc 2.1 Fast output for the package's own
test.fasta (container test, 2026-10-21), with funannotate-style ids.
"""
import gzip
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
SCRIPT = TESTS_DIR.parent / "bin" / "merge_deeploc.py"

HEADER = ("Protein_ID,Localizations,Signals,Membrane types,Cytoplasm,Nucleus,"
          "Extracellular,Cell membrane,Mitochondrion,Plastid,Endoplasmic reticulum,"
          "Lysosome/Vacuole,Golgi apparatus,Peroxisome,Peripheral,Transmembrane,"
          "Lipid anchor,Soluble")
ROWS_A = [
    "AAAA_000001-T1,Extracellular|Endoplasmic reticulum,Signal peptide,Soluble,"
    "0.26899999380111694,0.1404999941587448,0.707099974155426,0.2443999946117401,"
    "0.12189999967813492,0.010200000368058681,0.6635000109672546,0.4058000147342682,"
    "0.23389999568462372,0.004399999976158142,0.39430001378059387,0.06459999829530716,"
    "0.10559999942779541,0.8167999982833862",
    "AAAA_000002-T1,Peroxisome,Peroxisomal targeting signal,Soluble,"
    "0.40220001339912415,0.3111000061035156,0.13449999690055847,0.10530000180006027,"
    "0.3305000066757202,0.042399998754262924,0.13379999995231628,0.06419999897480011,"
    "0.01759999990463257,0.7369999885559082,0.5723000168800354,0.028999999165534973,"
    "0.0333000011742115,0.7875999808311462",
]
ROWS_B = [
    "BBBB_000007-T1,Mitochondrion,Mitochondrial transit peptide,Soluble,"
    "0.1,0.1,0.1,0.1,0.9,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.9",
]


def write_gz(path, header, rows):
    with gzip.open(path, "wt") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(r + "\n")


def run_merge(files, out):
    return subprocess.run([sys.executable, str(SCRIPT), "-o", str(out),
                           "--memory-limit", "1GB", "--threads", "1"]
                          + [str(f) for f in files],
                          capture_output=True, text=True, cwd=out.parent)


def query(parquet, sql):
    res = subprocess.run(["duckdb", "-json", "-c",
                          sql.replace("{T}", f"read_parquet('{parquet}')")],
                         capture_output=True, text=True, check=True)
    return json.loads(res.stdout or "[]")


def test_merge_columns_and_flags():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        a, b = d / "AAAA.deeploc_fast.csv.gz", d / "BBBB.deeploc_fast.csv.gz"
        write_gz(a, HEADER, ROWS_A)
        write_gz(b, HEADER, ROWS_B)
        out = d / "deeploc.parquet"
        r = run_merge([b, a], out)
        assert r.returncode == 0, r.stderr
        rows = query(out, "SELECT * FROM {T} ORDER BY protein_id")
        assert [x["protein_id"] for x in rows] == ["AAAA_000001-T1", "AAAA_000002-T1", "BBBB_000007-T1"]
        assert [x["species_prefix"] for x in rows] == ["AAAA", "AAAA", "BBBB"]
        assert all(x["model"] == "fast" for x in rows)
        assert [x["has_signal_peptide"] for x in rows] == [True, False, False]
        assert [x["has_mito_transit_peptide"] for x in rows] == [False, False, True]
        assert rows[0]["localizations"] == "Extracellular|Endoplasmic reticulum"
        assert abs(rows[0]["p_extracellular"] - 0.7071) < 1e-6
        assert abs(rows[1]["p_peroxisome"] - 0.737) < 1e-6
        assert abs(rows[2]["p_mitochondrion"] - 0.9) < 1e-6
        assert len(rows[0]) == 22


def test_rejects_mixed_models():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        a, b = d / "AAAA.deeploc_fast.csv.gz", d / "BBBB.deeploc_accurate.csv.gz"
        write_gz(a, HEADER, ROWS_A)
        write_gz(b, HEADER, ROWS_B)
        r = run_merge([a, b], d / "x.parquet")
        assert r.returncode != 0 and "mix DeepLoc models" in r.stderr + r.stdout


def test_rejects_changed_header():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        a = d / "AAAA.deeploc_fast.csv.gz"
        write_gz(a, HEADER.replace("Soluble", "Soluble,Extra"), ROWS_A)
        r = run_merge([a], d / "x.parquet")
        assert r.returncode != 0 and "unexpected header" in r.stderr + r.stdout


if __name__ == "__main__":
    if shutil.which("duckdb") is None:
        print("SKIP: duckdb CLI not on PATH (module load duckdb)")
        sys.exit(0)
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    for t in tests:
        t()
        print(f"  OK  {t.__name__}")
    print(f"All {len(tests)} tests passed.")
