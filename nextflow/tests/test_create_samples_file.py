"""Tests for scripts/create_samples_file.py.

Focus: the results-folder (`out`) machinery and the HAPLOTYPE dedup key, i.e.
the parts that decide whether two genomes stay two genomes. These are the paths
where a mistake is silent -- a dropped row or a shared annotation folder looks
like a smaller samples.csv, not like an error.
"""
import csv
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "create_samples_file.py"

ACC_COLS = ["ACCESSION", "SPECIES", "STRAIN", "NCBI_TAXID", "BIOPROJECT",
            "ASM_NAME", "ASM_FOLDER", "HAPLOTYPE"]
TAX_COLS = ["ASM_ACCESSION", "PHYLUM", "SUBPHYLUM", "CLASS", "SUBCLASS",
            "ORDER", "FAMILY", "GENUS", "SPECIES"]


def _load_module():
    """Import the script as a module so the pure helpers can be called."""
    spec = importlib.util.spec_from_file_location("create_samples_file", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


csf = _load_module()


# --- make_sample_tag: the `out` key -----------------------------------------
# This is a hand port of makeSampleTag/cleanStrain in
# nextflow/modules/common/utils.nf. The pipeline names each results folder with
# the Groovy version; create_samples_file.py only PREDICTS that name, to find
# collisions. If the two drift, the check passes and two genomes still share a
# folder -- so pin the contract, step for step.

# The four worked examples in the utils.nf docstring, treated as golden.
@pytest.mark.parametrize("species,strain,expected", [
    ("Saccharomyces cerevisiae", "CBS 1171", "Saccharomyces_cerevisiae_CBS_1171"),
    ("Aspergillus fumigatus", "Af293; CBS 101", "Aspergillus_fumigatus_Af293"),
    ("Fusarium oxysporum", "", "Fusarium_oxysporum"),
    ("Beauveria bassiana", "ARSEF 2860*", "Beauveria_bassiana_ARSEF_2860"),
])
def test_make_sample_tag_matches_utils_nf_examples(species, strain, expected):
    assert csf.make_sample_tag(species, strain) == expected


@pytest.mark.parametrize("species,strain,expected", [
    ("Aspergillus niger", "", "Aspergillus_niger"),
    ("", "CBS 1", "CBS_1"),
    ("Candida albicans", "A:B", "Candida_albicans_A_B"),               # ':' -> space -> '_'
    ("Saccharomyces 'cerevisiae'", "X", "Saccharomyces_cerevisiae_X"), # quotes dropped
    ("Fusarium sp.", "a/b", "Fusarium_sp._a_b"),                       # '/' -> '_'
    ("Fusarium sp.", "a[1]", "Fusarium_sp._a_1_"),                     # brackets -> '_'
    ("Fusarium sp.", "a*b", "Fusarium_sp._a-b"),                       # inner '*' -> '-'
    ("Fusarium sp.", "*star*", "Fusarium_sp._star"),                   # outer '*' dropped
])
def test_make_sample_tag(species, strain, expected):
    assert csf.make_sample_tag(species, strain) == expected


# Each of these was wrong in the first version of the port. They are the cases
# where a naive reading of the Groovy and the Groovy itself disagree, so they
# are the ones most likely to be broken again by a future "tidy-up".
@pytest.mark.parametrize("species,strain,expected,why", [
    ("Aspergillus niger", "CBS (1)", "Aspergillus_niger_CBS_1",
     "shell metacharacters become a space, then collapse to '_'"),
    ("Aspergillus niger", "a&b", "Aspergillus_niger_a_b", "'&' is a metacharacter"),
    ("Aspergillus niger", "a|b", "Aspergillus_niger_a_b", "'|' is a metacharacter"),
    ("Aspergillus niger", "a<b>c", "Aspergillus_niger_a_b_c", "'<' '>' are metacharacters"),
    ("Aspergillus niger", "O'Brien", "Aspergillus_niger_OBrien",
     "cleanStrain strips quotes from the STRAIN too, not just the species"),
    ("Aspergillus niger", 'X "1"', "Aspergillus_niger_X_1", "double quotes likewise"),
    ("Aspergillus niger", "X *", "Aspergillus_niger_X",
     "a trailing '*' takes the space before it (/\\*+\\s*$/), leaving no trailing '_'"),
    ("Aspergillus niger", "'Af293'; CBS 101", "Aspergillus_niger_Af293",
     "quotes are stripped BEFORE the ';' split, so the first token is clean"),
])
def test_make_sample_tag_groovy_divergences(species, strain, expected, why):
    assert csf.make_sample_tag(species, strain) == expected, why


def test_clean_strain_tag_step_order_matches_the_groovy():
    """':' becomes a space before the '*' rules run, so ':*' is a boundary."""
    assert csf.clean_strain_tag("a:*b") == "a-b"
    assert csf.clean_strain_tag("  *lead") == "lead"
    assert csf.clean_strain_tag("trail*  ") == "trail"


def test_make_sample_tag_separates_haplotype_suffixed_strains():
    """The whole point of out_suffix.csv: suffixed strains give distinct folders."""
    a = csf.make_sample_tag("Amanita phalloides", "gfAmaPhal2_hap1")
    b = csf.make_sample_tag("Amanita phalloides", "gfAmaPhal2_hap2")
    assert a != b


# --- fixture helpers ---------------------------------------------------------
def _acc_row(acc, species, strain, asm, haplotype="", bioproject="PRJNA1", taxid="1"):
    return {"ACCESSION": acc, "SPECIES": species, "STRAIN": strain,
            "NCBI_TAXID": taxid, "BIOPROJECT": bioproject, "ASM_NAME": asm,
            "ASM_FOLDER": f"{acc}_{asm}", "HAPLOTYPE": haplotype}


def _write_csv(path, cols, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow(r)


def _run(tmp_path, acc_rows, *, acc_cols=ACC_COLS, out_suffix=None,
         keep_dupes=None, overrides=None, extra_args=()):
    """Build fixtures, run the generator, return (CompletedProcess, rows, clash)."""
    acc = tmp_path / "accessions.csv"
    tax = tmp_path / "taxonomy.csv"
    out = tmp_path / "samples.csv"
    clash = tmp_path / "clash.csv"
    _write_csv(acc, acc_cols, acc_rows)
    _write_csv(tax, TAX_COLS, [
        {"ASM_ACCESSION": r["ASM_FOLDER"], "PHYLUM": "Ascomycota", "SUBPHYLUM": "Pezizomycotina",
         "CLASS": "Eurotiomycetes", "SUBCLASS": "", "ORDER": "Eurotiales",
         "FAMILY": "Aspergillaceae", "GENUS": r["SPECIES"].split()[0],
         "SPECIES": r["SPECIES"]}
        for r in acc_rows])

    args = [sys.executable, "-I", str(SCRIPT),
            "--accessions", str(acc), "--taxonomy", str(tax),
            "--outfile", str(out), "--clash-report", str(clash),
            # never fall through to the repo's real curation files
            "--exclude", str(tmp_path / "none.txt"),
            "--prefer", str(tmp_path / "none.txt"),
            "--keep-dupes", str(keep_dupes or tmp_path / "none.csv"),
            "--overrides", str(overrides or tmp_path / "none.csv"),
            "--out-suffix", str(out_suffix or tmp_path / "none.csv"),
            *extra_args]
    p = subprocess.run(args, capture_output=True, text=True)
    rows = list(csv.DictReader(open(out))) if out.exists() else []
    clash_rows = list(csv.DictReader(open(clash))) if clash.exists() else []
    return p, rows, clash_rows


# --- the HAPLOTYPE guard -----------------------------------------------------
def test_missing_haplotype_column_is_refused(tmp_path):
    """Without the column a hap1/hap2 pair collapses to one genome. Refuse."""
    cols = [c for c in ACC_COLS if c != "HAPLOTYPE"]
    rows = [{k: v for k, v in _acc_row("GCA_000000001.1", "Aspergillus niger",
                                       "S1", "asmA").items() if k != "HAPLOTYPE"}]
    p, out_rows, _ = _run(tmp_path, rows, acc_cols=cols)
    assert p.returncode != 0
    assert "HAPLOTYPE" in p.stderr
    assert out_rows == []


def test_missing_haplotype_column_allowed_with_flag(tmp_path):
    """extract_hand_edits.py regenerates from old snapshots and needs this."""
    cols = [c for c in ACC_COLS if c != "HAPLOTYPE"]
    rows = [{k: v for k, v in _acc_row("GCA_000000001.1", "Aspergillus niger",
                                       "S1", "asmA").items() if k != "HAPLOTYPE"}]
    p, out_rows, _ = _run(tmp_path, rows, acc_cols=cols,
                          extra_args=["--allow-missing-haplotype"])
    assert p.returncode == 0, p.stderr
    assert len(out_rows) == 1


def test_extract_hand_edits_passes_the_optout():
    """The opt-out must stay wired up, or that script breaks on old snapshots."""
    src = (REPO / "scripts" / "extract_hand_edits.py").read_text()
    assert "--allow-missing-haplotype" in src


# --- the dedup key -----------------------------------------------------------
def test_haplotype_pair_is_not_collapsed_by_dedup(tmp_path):
    """Same species+strain, different haplotype role: two genomes, not one.

    They still share an `out` key, so both are held back by the clash check --
    but the point here is that dedup did not silently drop one of them.
    """
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA", "hap1"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S1", "asmB", "hap2")]
    p, out_rows, clash = _run(tmp_path, rows)
    assert p.returncode == 0, p.stderr
    assert "auto-collapsed 0" in p.stderr
    assert {c["ASMID"] for c in clash} == {"GCA_000000001.1_asmA", "GCA_000000002.1_asmB"}
    assert all(c["STATUS"] == "held" for c in clash)
    assert out_rows == []


def test_plain_duplicate_is_still_collapsed(tmp_path):
    """No haplotype role: the ordinary dedup rule must still fire, newest wins."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S1", "asmB")]
    p, out_rows, _ = _run(tmp_path, rows)
    assert p.returncode == 0, p.stderr
    assert len(out_rows) == 1
    assert out_rows[0]["ASMID"] == "GCA_000000002.1_asmB"


# --- out_suffix.csv ----------------------------------------------------------
def test_out_suffix_separates_a_haplotype_pair(tmp_path):
    """With suffixes the pair gets two folders and both rows survive."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA", "hap1"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S1", "asmB", "hap2")]
    sfx = tmp_path / "out_suffix.csv"
    _write_csv(sfx, ["ASMID", "SUFFIX", "REASON"], [
        {"ASMID": "GCA_000000001.1_asmA", "SUFFIX": "hap1", "REASON": "pair"},
        {"ASMID": "GCA_000000002.1_asmB", "SUFFIX": "hap2", "REASON": "pair"}])
    p, out_rows, clash = _run(tmp_path, rows, out_suffix=sfx)
    assert p.returncode == 0, p.stderr
    assert {r["STRAIN"] for r in out_rows} == {"S1_hap1", "S1_hap2"}
    assert clash == []


def test_out_suffix_on_empty_strain_becomes_the_strain(tmp_path):
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "", "asmA", "hap1")]
    sfx = tmp_path / "out_suffix.csv"
    _write_csv(sfx, ["ASMID", "SUFFIX", "REASON"],
               [{"ASMID": "GCA_000000001.1_asmA", "SUFFIX": "hap1", "REASON": "x"}])
    p, out_rows, _ = _run(tmp_path, rows, out_suffix=sfx)
    assert p.returncode == 0, p.stderr
    assert out_rows[0]["STRAIN"] == "hap1"


def test_hold_keeps_the_row_out_until_someone_chooses(tmp_path):
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA"),
            _acc_row("GCA_000000009.1", "Penicillium rubens", "S9", "asmZ")]
    sfx = tmp_path / "out_suffix.csv"
    _write_csv(sfx, ["ASMID", "SUFFIX", "REASON"],
               [{"ASMID": "GCA_000000001.1_asmA", "SUFFIX": "HOLD", "REASON": "needs a call"}])
    p, out_rows, clash = _run(tmp_path, rows, out_suffix=sfx)
    assert p.returncode == 0, p.stderr
    assert {r["ASMID"] for r in out_rows} == {"GCA_000000009.1_asmZ"}
    held = [c for c in clash if c["ASMID"] == "GCA_000000001.1_asmA"]
    assert len(held) == 1
    assert held[0]["STATUS"].startswith("held (out_suffix HOLD)")
    assert "needs a call" in held[0]["STATUS"]


def test_hold_is_case_insensitive(tmp_path):
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA")]
    sfx = tmp_path / "out_suffix.csv"
    _write_csv(sfx, ["ASMID", "SUFFIX", "REASON"],
               [{"ASMID": "GCA_000000001.1_asmA", "SUFFIX": "hold", "REASON": "x"}])
    p, out_rows, _ = _run(tmp_path, rows, out_suffix=sfx)
    assert p.returncode == 0, p.stderr
    assert out_rows == []


def test_out_suffix_comments_are_skipped(tmp_path):
    """The generated file starts with '#' banner lines."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA")]
    sfx = tmp_path / "out_suffix.csv"
    sfx.write_text("# a banner\n# another\nASMID,SUFFIX,REASON\n"
                   "GCA_000000001.1_asmA,hapX,why\n")
    p, out_rows, _ = _run(tmp_path, rows, out_suffix=sfx)
    assert p.returncode == 0, p.stderr
    assert out_rows[0]["STRAIN"] == "S1_hapX"


# --- the out-clash hold-back -------------------------------------------------
def test_clash_is_held_back_by_default(tmp_path):
    """Two different isolates that sanitize to one folder key: hold both."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S/1", "asmA"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S_1", "asmB")]
    p, out_rows, clash = _run(tmp_path, rows)
    assert p.returncode == 0, p.stderr
    assert out_rows == []
    assert len(clash) == 2
    assert {c["OUT"] for c in clash} == {"Aspergillus_niger_S_1"}
    assert all(c["N_IN_GROUP"] == "2" for c in clash)


def test_allow_out_clash_keeps_them(tmp_path):
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S/1", "asmA"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S_1", "asmB")]
    p, out_rows, clash = _run(tmp_path, rows, extra_args=["--allow-out-clash"])
    assert p.returncode == 0, p.stderr
    assert len(out_rows) == 2
    assert all(c["STATUS"] == "kept" for c in clash)


def test_keep_dupes_isolate_warns_but_is_not_held(tmp_path):
    """Intentional multi-assembly isolates are already in the runs; dropping
    them would delete rows that exist on disk. Warn instead."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S1", "asmB")]
    kd = tmp_path / "keep_dupes.csv"
    _write_csv(kd, ["SPECIES", "STRAIN"], [{"SPECIES": "Aspergillus niger", "STRAIN": "S1"}])
    p, out_rows, clash = _run(tmp_path, rows, keep_dupes=kd)
    assert p.returncode == 0, p.stderr
    assert len(out_rows) == 2, "keep_dupes rows must survive both dedup and the clash check"
    assert "[out-clash] WARNING" in p.stderr
    assert clash == []


def test_clash_report_header(tmp_path):
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA")]
    p, _, _ = _run(tmp_path, rows)
    assert p.returncode == 0, p.stderr
    hdr = (tmp_path / "clash.csv").read_text().splitlines()[0]
    assert hdr == "ASMID,OUT,HAPLOTYPE,N_IN_GROUP,SPECIES,STRAIN,STATUS"


def test_haplotype_is_not_written_to_samples_csv(tmp_path):
    """HAPLOTYPE is a dedup/report key only; samples.csv keeps its column set."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA", "hap1")]
    p, out_rows, _ = _run(tmp_path, rows)
    assert p.returncode == 0, p.stderr
    assert list(out_rows[0].keys()) == csf.OUT_FIELDS


# --- ordering: overrides win, and they move the folder key -------------------
def test_override_is_applied_before_the_clash_check(tmp_path):
    """A SPECIES override renames the folder, so it must land before the check."""
    rows = [_acc_row("GCA_000000001.1", "Aspergillus niger", "S1", "asmA", "hap1"),
            _acc_row("GCA_000000002.1", "Aspergillus niger", "S1", "asmB", "hap2")]
    ov = tmp_path / "overrides.csv"
    _write_csv(ov, ["ASMID", "FIELD", "VALUE"],
               [{"ASMID": "GCA_000000002.1_asmB", "FIELD": "SPECIES",
                 "VALUE": "Aspergillus awamori"}])
    p, out_rows, clash = _run(tmp_path, rows, overrides=ov)
    assert p.returncode == 0, p.stderr
    assert clash == [], "the override gave them different folder keys"
    assert {r["SPECIES"] for r in out_rows} == {"Aspergillus niger", "Aspergillus awamori"}
