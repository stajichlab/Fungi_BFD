"""Experiment A (DECISIONS D81): within-genome titration of the number of complete
PASA training models. Shared helpers for make_inputs / subset.

Code under test is the frozen snapshot predict_arms/code_new4 (complete-first
getBestModel, R3/R5, review fixes, R6 option (b) on by default).
"""
import gzip, logging, os, random, shutil, sys

A = "/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/predict_arms"
X = os.path.dirname(A)
R = os.path.dirname(X)
CODE = f"{A}/code_new4"
sys.path.insert(0, CODE)
logging.basicConfig(level=logging.WARNING)
from funannotate import library as lib  # noqa: E402
from funannotate import train  # noqa: E402

assert lib.__file__.startswith(CODE), lib.__file__
lib.log = logging.getLogger("titration")

T = f"{A}/titration"
GENOMES = {
    "Neurospora_crassa_OR74A": "GCF_000182925.2_NC12",
    "Aspergillus_nidulans_FGSC_A4": "GCF_000011425.1_ASM1142v1",
    "Botrytis_cinerea_B05.10": "GCF_000143535.2_ASM14353v4",
    "Cryptococcus_neoformans_H99": "GCF_000149245.1_CNA3",
}
NS = [50, 100, 200, 300, 500, 750, 1000, 2000]


def draws_for(n):
    return 5 if n <= 500 else 3


def seed_for(genome, n, draw):
    return (sum(map(ord, genome)) * 1000003 + n * 101 + draw) % (2 ** 31)
