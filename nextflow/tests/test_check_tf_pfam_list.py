"""Tests for bin/check_tf_pfam_list.py. Run: python3 -m pytest nextflow/tests/test_check_tf_pfam_list.py

The real-dat test is skipped when the Pfam dat file is absent.
"""
import subprocess
import sys
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).resolve().parent
NF = TESTS_DIR.parent
SCRIPT = NF / "bin" / "check_tf_pfam_list.py"
REAL_CSV = NF / "assets" / "fungal_tf_pfam_domains.csv"
REAL_DAT = Path("/srv/projects/db/pfam/current/Pfam-A.hmm.dat")

DAT = """# STOCKHOLM 1.0
#=GF ID   Zn_clus
#=GF AC   PF00172.25
#=GF DE   Fungal Zn(2)-Cys(6) binuclear cluster domain
//
# STOCKHOLM 1.0
#=GF ID   GnRH
#=GF AC   PF00446.23
#=GF DE   Gonadotropin-releasing hormone
//
"""


def run(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)],
                          capture_output=True, text=True)


def make(tmp_path, rows):
    dat = tmp_path / "Pfam-A.hmm.dat"
    dat.write_text(DAT)
    csv_path = tmp_path / "tf.csv"
    csv_path.write_text("pfam_acc,pfam_name,tf_family\n" + "".join(r + "\n" for r in rows))
    return csv_path, dat


def test_match(tmp_path):
    c, d = make(tmp_path, ["PF00172,Zn_clus,Zn cluster"])
    r = run("--csv", c, "--pfam-dat", d)
    assert r.returncode == 0, r.stderr
    assert "PF00172\tZn_clus\tZn_clus\tok" in r.stdout


def test_wrong_name_fails(tmp_path):
    c, d = make(tmp_path, ["PF00446,SBP,SBP zinc finger"])
    r = run("--csv", c, "--pfam-dat", d)
    assert r.returncode == 1
    assert "GnRH" in r.stdout and "FAIL name differs" in r.stdout


def test_missing_accession_fails(tmp_path):
    c, d = make(tmp_path, ["PF99999,Nothing,none"])
    r = run("--csv", c, "--pfam-dat", d)
    assert r.returncode == 1
    assert "absent from Pfam dat" in r.stdout


def test_explicit_missing_dat_is_error(tmp_path):
    c, _ = make(tmp_path, ["PF00172,Zn_clus,Zn cluster"])
    r = run("--csv", c, "--pfam-dat", tmp_path / "nope.dat")
    assert r.returncode == 2


@pytest.mark.skipif(not REAL_DAT.exists(), reason="Pfam dat not available")
def test_real_csv_against_real_dat():
    r = run("--csv", REAL_CSV, "--pfam-dat", REAL_DAT)
    assert r.returncode == 0, r.stdout + r.stderr
