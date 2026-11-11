"""Tests for the combined telomere call in bin/summarize_telomeres.py."""
import gzip
import importlib.util
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
_spec = importlib.util.spec_from_file_location("summarize_telomeres", BIN / "summarize_telomeres.py")
st = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(st)

HDR = ("ASMID\tunit\tperiod\tn_ends\tn_ends_examined\tfrac_ends\tmedian_tract_len\t"
       "bg_windows\tbg_hits\tenrichment\tstrand_biased_ends\tstatus\n")


def write_denovo(dirpath, asmid, cand_unit, n_ends, status, biased=0):
    rows = HDR
    if cand_unit:
        rows += f"{asmid}\t{cand_unit}\t{len(cand_unit)}\t{n_ends}\t20\t0.5\t100\t300\t0\tbg0\t{biased}\tcandidate\n"
    rows += f"{asmid}\t.\t0\t0\t20\t0\t0\t300\t0\tNA\t{biased}\t{status}\n"
    (dirpath / f"{asmid}.telomere_denovo.tsv").write_text(rows)


def tract_path(dirpath, asmid):
    return str(dirpath / f"{asmid}.telomeres.tsv.gz")


def test_load_denovo_missing_returns_empty(tmp_path):
    assert st.load_denovo(tract_path(tmp_path, "A1")) == {}
    assert st.telomere_call(0, {}) == ""
    assert st.telomere_call(3, {}) == "motif_found"


def test_calls(tmp_path):
    write_denovo(tmp_path, "L", "AGGGTTAGTT", 12, "repeat_found")
    d = st.load_denovo(tract_path(tmp_path, "L"))
    assert d["telomere_denovo_unit_len"] == 10 and d["telomere_denovo_n_ends"] == 12
    assert st.telomere_call(0, d) == "denovo_found"
    assert st.telomere_call(5, d) == "motif_found"

    write_denovo(tmp_path, "S", "ATT", 30, "repeat_found")
    assert st.telomere_call(0, st.load_denovo(tract_path(tmp_path, "S"))) == "short_unit_only"

    write_denovo(tmp_path, "B", None, 0, "strand_biased_only", biased=4)
    d = st.load_denovo(tract_path(tmp_path, "B"))
    assert d["telomere_denovo_strand_biased_ends"] == 4
    assert st.telomere_call(0, d) == "strand_biased_only"

    write_denovo(tmp_path, "N", None, 0, "no_terminal_repeat")
    assert st.telomere_call(0, st.load_denovo(tract_path(tmp_path, "N"))) == "no_terminal_repeat"
