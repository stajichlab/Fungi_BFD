"""Tests for bin/merge_filter.py, the Python twin of loadMergeAllowedKeys() /
mergeFileAllowed() in modules/common/utils.nf. Run: python3 -m pytest nextflow/tests/test_merge_filter.py
"""
import importlib.util
import sys
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(BIN))
_spec = importlib.util.spec_from_file_location("merge_filter", BIN / "merge_filter.py")
mf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mf)

SAMPLES = """ASMID,SPECIES,STRAIN,LOCUSTAG
GCA_000001.1_ASM1v1,Fusarium graminearum,PH-1,F0000001
GCA_000002.2_ASM2v2,Saccharomyces cerevisiae,CEN.PK113-7D,F0000002
GCA_000003.1_ASM3v1,Saccharomyces cerevisiae,CEN.PK113-7D,F0000003
GCA_000004.1_ASM4v1,Blumeria graminis,JIW2,F0000004
GCA_000005.1_ASM5v1,Agmasoma carcini nom. nud.,v1,F0000005
"""
SUPPRESS = """# comment
GCA_000003.1_ASM3v1,duplicate of GCA_000002.2 (same strain)
GCA_000004.1_ASM4v1,too small
"""


def keys(tmp_path):
    s = tmp_path / "samples.csv"; s.write_text(SAMPLES)
    p = tmp_path / "suppress.txt"; p.write_text(SUPPRESS)
    return mf.allowed_keys(str(s), str(p))


def test_locustag_named_files(tmp_path):
    k = keys(tmp_path)
    assert mf.file_allowed(k, "F0000001.signalp.gff3.gz")
    assert mf.file_allowed(k, "F0000001.gene_exons.csv.gz")
    assert not mf.file_allowed(k, "F0000004.signalp.gff3.gz")      # suppressed
    assert not mf.file_allowed(k, "F9999999.signalp.gff3.gz")      # not in samples.csv


def test_old_name_copy_is_rejected(tmp_path):
    # LOCUSTAG followed by '_' (an older file-name format) must not match.
    k = keys(tmp_path)
    assert not mf.file_allowed(k, "F0000001_CBS_101.48_substr._RVII-324_met-.signalp.gff3.gz")
    assert not mf.file_allowed(k, "F0000001_RSA_919-.aiupred.txt.gz")


def test_targetp_summary_name(tmp_path):
    k = keys(tmp_path)
    assert mf.file_allowed(k, "F0000001_summary.targetp2.gz")
    assert not mf.file_allowed(k, "F0000004_summary.targetp2.gz")
    assert not mf.file_allowed(k, "F0000001_CBS_101.48_summary.targetp2.gz")


def test_asmid_named_files(tmp_path):
    k = keys(tmp_path)
    assert mf.file_allowed(k, "GCA_000001.1_ASM1v1.stats.txt")
    assert mf.file_allowed(k, "GCA_000001.1_ASM1v1.BUSCO_summary.fungi_odb12.txt")
    assert not mf.file_allowed(k, "GCA_000004.1_ASM4v1.telomeres.tsv.gz")   # suppressed
    # a prefix of an ASMID that ends at its version dot is not a key
    assert not mf.file_allowed(k, "GCA_000009.1_ASM9v1.stats.txt")


def test_out_named_files_and_shared_out(tmp_path):
    k = keys(tmp_path)
    assert mf.file_allowed(k, "Fusarium_graminearum_PH-1.overview.tsv.gz")
    # S. cerevisiae CEN.PK113-7D: one assembly suppressed, one kept, same out tag
    assert mf.file_allowed(k, "Saccharomyces_cerevisiae_CEN.PK113-7D.overview.tsv.gz")
    # Blumeria JIW2: its only assembly is suppressed
    assert not mf.file_allowed(k, "Blumeria_graminis_JIW2.rfam.tblout.gz")


def test_out_tag_with_dots(tmp_path):
    k = keys(tmp_path)
    tag = mf.make_sample_tag("Agmasoma carcini nom. nud.", "v1")
    assert tag in k
    assert mf.file_allowed(k, f"{tag}.substrates.tsv.gz")


def test_suppress_comments_and_blank(tmp_path):
    p = tmp_path / "s.txt"; p.write_text("# x\n\nGCA_1.1_A,why\n")
    assert mf.load_suppress(str(p)) == {"GCA_1.1_A"}
    assert mf.load_suppress(str(tmp_path / "missing.txt")) == set()
