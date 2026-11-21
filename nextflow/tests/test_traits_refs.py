"""Tests for the reference files in lib/refs (FungalTraits, growth-form overrides, class-to-subphylum)."""
import csv
from pathlib import Path

REFS = Path(__file__).resolve().parents[2] / "lib" / "refs"


def _rows(name, delim):
    with open(REFS / name, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter=delim))


def test_ft_csv_is_utf8_with_genus_column():
    rows = _rows("fungaltraits_1.2.utf8.csv", ",")
    assert len(rows) >= 10000
    assert {"GENUS", "primary_lifestyle", "jrk_template"} <= set(rows[0])


def test_class_subphylum_has_header_and_unique_classes():
    rows = _rows("class_subphylum.tsv", "\t")
    assert list(rows[0]) == ["CLASS", "SUBPHYLUM", "source", "note"]
    classes = [r["CLASS"] for r in rows]
    assert len(classes) == len(set(classes))
    assert all(r["SUBPHYLUM"] for r in rows)


def test_class_subphylum_sources_are_known():
    rows = _rows("class_subphylum.tsv", "\t")
    assert {r["source"] for r in rows} <= {"taxonomy_py", "atlas", "proposed", "ncbi"}


def test_data_spellings_are_mapped():
    m = {r["CLASS"]: r["SUBPHYLUM"] for r in _rows("class_subphylum.tsv", "\t")}
    assert m["Pneumocystomycetes"] == "Taphrinomycotina"
    assert m["Candelariomycetes"] == "Pezizomycotina"
    assert m["Spiculogloeomycetes"] == "Pucciniomycotina"
    assert m["Moniliellomycetes"] == "Ustilaginomycotina"
    assert "Pneumocystidomycetes" not in m


def test_overrides_file_has_genus_and_growth_form():
    rows = _rows("growth_form_overrides.tsv", "\t")
    assert len(rows) > 100 and {"GENUS", "growth_form"} <= set(rows[0])


def test_kickxellomycota_classes_map_to_kickxellomycotina():
    m = {r["CLASS"]: r["SUBPHYLUM"] for r in _rows("class_subphylum.tsv", "\t")}
    for c in ("Kickxellomycetes", "Harpellomycetes", "Dimargaritomycetes"):
        assert m.get(c) == "Kickxellomycotina", c


def test_phylum_subphylum_file_and_class_rows():
    ph = {r["PHYLUM"]: r["SUBPHYLUM"] for r in _rows("phylum_subphylum.tsv", "\t")}
    assert ph == {"Kickxellomycota": "Kickxellomycotina", "Mortierellomycota": "Mortierellomycotina",
                  "Glomeromycota": "Glomeromycotina", "Basidiobolomycota": "Basidiobolomycotina",
                  "Blastocladiomycota": "Blastocladiomycotina", "Microsporidia": "Microsporidiomycotina",
                  "Rozellomycota": "Rozellomycotina", "Aphelidiomycota": "Aphelidiomycotina"}
    m = {r["CLASS"]: r["SUBPHYLUM"] for r in _rows("class_subphylum.tsv", "\t")}
    assert m["Mortierellomycetes"] == "Mortierellomycotina" and m["Glomeromycetes"] == "Glomeromycotina"


def test_ncbi_classes_and_neocallimastigomycetes_are_mapped():
    rows = {r["CLASS"]: r for r in _rows("class_subphylum.tsv", "\t")}
    want = {"Chytridiomycetes": "Chytridiomycotina", "Entomophthoromycetes": "Entomophthoromycotina",
            "Peribolosporomycetes": "Ustilaginomycotina", "Novakomycetes": "Taphrinomycotina",
            "Xylobotryomycetes": "Pezizomycotina", "Caulochytriomycetes": "Chytridiomycotina",
            "Olpidiomycetes": "Olpidiomycotina", "Calcarisporiellomycetes": "Calcarisporiellomycotina"}
    for c, sub in want.items():
        assert rows[c]["SUBPHYLUM"] == sub and rows[c]["source"] == "ncbi", c
    assert rows["Neocallimastigomycetes"]["SUBPHYLUM"] == "Neocallimastigomycotina"


def test_early_diverging_classes_have_their_own_subphylum():
    m = {r["CLASS"]: r["SUBPHYLUM"] for r in _rows("class_subphylum.tsv", "\t")}
    assert m["Sanchytriomycetes"] == "Sanchytriomycotina"
    assert m["Monoblepharidomycetes"] == "Monoblepharidomycotina"
    assert m["Hyaloraphidiomycetes"] == "Hyaloraphidiomycotina"
