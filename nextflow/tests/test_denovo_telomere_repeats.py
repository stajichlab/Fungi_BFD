"""Tests for bin/denovo_telomere_repeats.py. Needs numpy."""
import importlib.util
import random
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
_spec = importlib.util.spec_from_file_location("denovo", BIN / "denovo_telomere_repeats.py")
d = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(d)
DET = dict(pmin=4, pmax=40, purity=0.85, tol=50, min_len=60)


def rnd(n, seed=3):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def test_canonical_strand_and_rotation():
    assert d.canonical("CCCTAA") == d.canonical("TTAGGG")
    u = "AAGGATGTCACGATCATTGGTGT"
    assert d.canonical(u) == d.canonical(d.rc(u))


def test_detects_short_and_long_units():
    r = d.detect("CCCTAA" * 20 + rnd(1000), **DET)
    assert r and d.canonical(r[3]) == d.canonical("TTAGGG")
    u = "AAGGATGTCACGATCATTGGTGT"
    r = d.detect(u * 6 + rnd(1000), **DET)
    assert r and d.canonical(r[3]) == d.canonical(u)


def test_random_sequence_gives_nothing():
    assert d.detect(rnd(1500), **DET) is None


def test_repeat_far_from_edge_is_ignored():
    assert d.detect(rnd(500) + "TTAGGG" * 30 + rnd(500), **DET) is None
