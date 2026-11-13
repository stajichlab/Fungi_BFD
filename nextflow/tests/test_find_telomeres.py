"""Tests for bin/find_telomeres.py orientation, labels and de-duplication.
Run: python3 -m pytest nextflow/tests/test_find_telomeres.py  (needs biopython)
"""
import importlib.util
import io
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(BIN))
_spec = importlib.util.spec_from_file_location("find_telomeres", BIN / "find_telomeres.py")
ft = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ft)


def run(tmp_path, seq, patterns, both_ends):
    fa = tmp_path / "g.fa"
    fa.write_text(f">s1\n{seq}\n")
    out = tmp_path / "o.tsv"
    argv = [str(fa), "-o", str(out), "-p", *patterns, "-n", "3", "-l", "36",
            "-s", "5000"]
    if both_ends:
        argv.append("--both-ends")
    ft.main(argv)
    rows = [l.split("\t") for l in out.read_text().splitlines()[1:]]
    return rows


def test_gstrand_label_collapses_orientation():
    assert ft.gstrand_label("CCCTAA") == "TTAGGG"
    assert ft.gstrand_label("TTAGGG") == "TTAGGG"
    assert ft.gstrand_label("CCC[C]+TAA") == "TTA[G]+GGG"


def test_pattern_and_rc_deduplicate_with_both_ends():
    assert len(ft.build_pattern_set(["TTAGGG", "CCCTAA"], True, False)) == 4


def test_g_strand_pattern_found_at_5prime_with_both_ends(tmp_path):
    # 12-bp Flavi unit in C-strand orientation at the 5' end needs both-ends
    # when the pattern is supplied in G-strand form.
    seq = "TGTTGACCCTAA" * 5 + "ACGT" * 100
    assert run(tmp_path, seq, ["TTAGGGTCAACA"], False) == []
    rows = run(tmp_path, seq, ["TTAGGGTCAACA"], True)
    assert len(rows) == 1 and rows[0][3] == "TTAGGGTCAACA"


def test_overlapping_patterns_report_one_tract(tmp_path):
    seq = "CCCTAA" * 10 + "ACGT" * 100
    rows = run(tmp_path, seq, ["CCCTAA", "CCC[C]+TAA"], True)
    assert len(rows) == 1
