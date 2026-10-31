"""Representative selection must prefer genomes that pass predict's pre-flight."""
import importlib.util
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bin"))

_spec = importlib.util.spec_from_file_location(
    "pick_representative_strain",
    Path(__file__).resolve().parents[1] / "bin" / "pick_representative_strain.py")
prs = importlib.util.module_from_spec(_spec)
sys.modules["pick_representative_strain"] = prs
_spec.loader.exec_module(prs)

PF = dict(prs.DEFAULT_PREFLIGHT)


def _stats(total, contigs, n50, busco, max_ctg=50_000):
    return {"complete_pct": busco, "n50_bp": n50, "total_bp": total,
            "contigs": contigs, "max_contig_bp": max_ctg}


def test_small_fragmented_fails_preflight():
    assert not prs.passes_preflight(_stats(4_271_025, 746, 7_226, 90.0), PF)


def test_small_but_contiguous_passes():
    assert prs.passes_preflight(_stats(5_610_204, 158, 97_321, 80.0), PF)


def test_large_fragmented_passes():
    assert prs.passes_preflight(_stats(30_000_000, 5000, 3_000, 90.0), PF)


def test_missing_stats_fail_open():
    assert prs.passes_preflight(None, PF)
    assert prs.passes_preflight({"complete_pct": 1.0, "n50_bp": 0, "total_bp": None}, PF)


def test_too_small_and_no_training_contig_fail():
    assert not prs.passes_preflight(_stats(50_000, 3, 20_000, 1.0), PF)
    assert not prs.passes_preflight(_stats(5_000_000, 900, 9_000, 1.0, max_ctg=8_000), PF)


def test_zero_threshold_disables_its_gate():
    no_train_gate = dict(PF, min_training_contigs=0)
    assert prs.passes_preflight(_stats(20_000_000, 900, 20_000, 1.0, max_ctg=8_000),
                                no_train_gate)
    no_small_gate = dict(PF, min_bp=0)
    assert prs.passes_preflight(_stats(4_271_025, 746, 7_226, 90.0), no_small_gate)
    no_floor = dict(PF, abs_min_bp=0)
    assert prs.passes_preflight(_stats(50_000, 3, 20_000, 1.0), no_floor)


def test_fallback_defaults_match_nextflow_config():
    # The pipeline passes params.predict_*; the script defaults are only for
    # hand runs, but must not drift from the values predict uses.
    import re
    cfg = (Path(__file__).resolve().parents[1] / "nextflow.config").read_text()
    def val(name):
        return int(re.search(rf"^\s*{name}\s*=\s*(\d+)", cfg, re.M).group(1))
    assert PF == {"min_bp": val("predict_min_asm_bp"),
                  "max_n50": val("predict_frag_max_n50"),
                  "max_contigs": val("predict_frag_max_contigs"),
                  "min_contig_len": val("predict_min_training_contig_len"),
                  "min_training_contigs": val("predict_min_training_contigs"),
                  "abs_min_bp": val("predict_abs_min_asm_bp")}


def test_picker_skips_higher_busco_unpredictable_candidate():
    busco = {"A": _stats(4_271_025, 746, 7_226, 95.0),   # best BUSCO, small_fragmented
             "B": _stats(5_610_204, 158, 97_321, 80.0)}
    samples = {"A": {"out": "sp_A"}, "B": {"out": "sp_B"}}
    assert prs.pick_representative(["A", "B"], busco, samples, {}) == "A"
    assert prs.pick_representative(["A", "B"], busco, samples, {}, PF) == "B"


def test_picker_keeps_old_choice_when_no_candidate_passes():
    busco = {"A": _stats(4_000_000, 700, 7_000, 95.0),
             "B": _stats(4_100_000, 800, 6_000, 80.0)}
    samples = {"A": {"out": "sp_A"}, "B": {"out": "sp_B"}}
    assert prs.pick_representative(["A", "B"], busco, samples, {}, PF) == "A"


def test_loader_keeps_genomes_without_busco_row(tmp_path):
    import duckdb
    con = duckdb.connect()
    con.execute("CREATE TABLE b AS SELECT 'A' AS ASMID, 90.0 AS complete_pct")
    con.execute("CREATE TABLE a AS SELECT * FROM (VALUES "
                "('A', 100, 1000000, 5000, 20000), ('M', 746, 4271025, 7226, 46788)) "
                "t(ASMID, contig_count, total_length_bp, N50_bp, max_contig_bp)")
    con.execute(f"COPY b TO '{tmp_path}/busco_genome.parquet' (FORMAT parquet)")
    con.execute(f"COPY a TO '{tmp_path}/asm_stats.parquet' (FORMAT parquet)")
    got = prs.load_busco_by_asmid(tmp_path)
    assert got["A"]["complete_pct"] == 90.0 and got["A"]["n50_bp"] == 5000
    # no BUSCO row: old ranking default, but real assembly stats for preflight
    assert got["M"]["complete_pct"] == -1.0 and got["M"]["n50_bp"] == -1
    assert got["M"]["asm_n50"] == 7226
    assert not prs.passes_preflight(got["M"], PF)
