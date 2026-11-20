"""Tests for sql/bfd_traits.sql and the v_genome_labels join, on small synthetic tables."""
import csv
import re
from pathlib import Path

import duckdb

# nextflow/tests/ -> parents[1] is nextflow/, parents[2] is the repo root.
# sql/ lives at the repo root; the helper scripts live in nextflow/bin/.
REPO = Path(__file__).resolve().parents[2]
BIN = Path(__file__).resolve().parents[1] / "bin"
TRAITS_SQL = REPO / "sql" / "bfd_traits.sql"
EXTRA_SQL = REPO / "sql" / "bfd_extra_views.sql"

FT_COLS = ["jrk_template", "GENUS", "primary_lifestyle", "Secondary_lifestyle",
           "Plant_pathogenic_capacity_template", "Decay_substrate_template",
           "Animal_biotrophic_capacity_template", "Growth_form_template"]


def _write(path, header, rows, delim):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter=delim)
        w.writerow(header)
        w.writerows(rows)


def ft(jrk, genus, life=None, growth=None):
    return [jrk, genus, life, None, None, None, None, growth]


def write_refs(tmp_path, ft_rows=(), ovr_rows=(), cls_rows=(), phy_rows=()):
    _write(tmp_path / "fungaltraits_1.2.utf8.csv", FT_COLS, ft_rows, ",")
    _write(tmp_path / "growth_form_overrides.tsv", ["GENUS", "growth_form", "source", "n_genomes", "CLASS"],
           ovr_rows, "\t")
    _write(tmp_path / "class_subphylum.tsv", ["CLASS", "SUBPHYLUM", "source", "note"], cls_rows, "\t")
    _write(tmp_path / "phylum_subphylum.tsv", ["PHYLUM", "SUBPHYLUM", "source", "note"], phy_rows, "\t")


def apply_traits(con, tmp_path):
    sql = TRAITS_SQL.read_text()
    first = sql.splitlines()[0]
    assert first.startswith("-- REF_DIR="), "line 1 of bfd_traits.sql must be '-- REF_DIR=<path>'"
    con.execute(sql.replace(first[len("-- REF_DIR="):], str(tmp_path)))


def load(tmp_path, species, ft_rows=(), ovr_rows=(), cls_rows=(), phy_rows=()):
    """species: list of (LOCUSTAG, GENUS, CLASS, SUBPHYLUM[, PHYLUM]). Returns a connection with the traits SQL applied."""
    write_refs(tmp_path, ft_rows, ovr_rows, cls_rows, phy_rows)
    con = duckdb.connect()
    con.execute("CREATE TABLE species (LOCUSTAG VARCHAR, GENUS VARCHAR, CLASS VARCHAR, SUBPHYLUM VARCHAR, PHYLUM VARCHAR)")
    con.executemany("INSERT INTO species VALUES (?,?,?,?,?)", [tuple(r) + (None,) * (5 - len(r)) for r in species])
    apply_traits(con, tmp_path)
    return con


def one(con, locustag):
    cur = con.execute("SELECT * FROM v_genome_traits WHERE LOCUSTAG = ?", [locustag])
    cols = [d[0] for d in cur.description]
    rows = cur.fetchall()
    assert len(rows) == 1, f"{locustag}: {len(rows)} rows"
    return dict(zip(cols, rows[0]))


def test_duplicate_genus_keeps_first_by_jrk_template(tmp_path):
    con = load(tmp_path, [("A1", "Mucor", None, None), ("A2", "Mucor", None, None)],
               ft_rows=[ft("jrk2", "Mucor", "wrong"), ft("jrk1", "Mucor", "right")])
    assert con.execute("SELECT count(*) FROM v_genome_traits").fetchone()[0] == 2
    assert one(con, "A1")["primary_lifestyle"] == "right"


def test_override_wins_and_sets_source(tmp_path):
    con = load(tmp_path, [("A1", "Candida", None, None), ("A2", "Mucor", None, None)],
               ft_rows=[ft("j1", "Candida", "x", "filamentous_mycelium"), ft("j2", "Mucor", "y", "filamentous_mycelium")],
               ovr_rows=[("Candida", "yeast", "class_rule", "9", "Saccharomycetes"),
                         ("Mucor", "", "class_rule", "1", "Mucoromycetes")])
    a1, a2 = one(con, "A1"), one(con, "A2")
    assert (a1["growth_form"], a1["growth_form_source"]) == ("yeast", "override")
    assert (a2["growth_form"], a2["growth_form_source"]) == ("filamentous_mycelium", "fungaltraits")


def test_unmatched_genus_and_null_genus(tmp_path):
    con = load(tmp_path, [("A1", "Nogenus", None, None), ("A2", None, None, None)], ft_rows=[ft("j1", "Mucor", "y")])
    for lt in ("A1", "A2"):
        r = one(con, lt)
        assert r["in_fungaltraits"] is False
        assert r["primary_lifestyle"] is None and r["growth_form"] is None and r["growth_form_source"] is None


def test_subphylum_from_class_table(tmp_path):
    cls = [("Eurotiomycetes", "Pezizomycotina", "taxonomy_py", "")]
    con = load(tmp_path, [("A1", "G", "Eurotiomycetes", None), ("A2", "G", "Unlisted", None), ("A3", "G", None, None)],
               cls_rows=cls)
    r = one(con, "A1")
    assert (r["SUBPHYLUM"], r["subphylum_source"]) == ("Pezizomycotina", "class_table")
    for lt in ("A2", "A3"):
        r = one(con, lt)
        assert r["SUBPHYLUM"] is None and r["subphylum_source"] is None


def test_species_subphylum_wins(tmp_path):
    cls = [("Eurotiomycetes", "Pezizomycotina", "taxonomy_py", "")]
    con = load(tmp_path, [("A1", "G", "Eurotiomycetes", "Other-sub")], cls_rows=cls)
    r = one(con, "A1")
    assert (r["SUBPHYLUM"], r["subphylum_source"]) == ("Other-sub", "species")


def test_yeast_origin_matches_atlas_rules(tmp_path):
    spec = [("A1", "Gy", "Saccharomycetes", None), ("A2", "Gd", "Eurotiomycetes", None),
            ("A3", "Gf", "Sordariomycetes", None), ("A4", "Gn", "Sordariomycetes", None)]
    con = load(tmp_path, spec, ft_rows=[ft("j1", "Gy", None, "yeast"), ft("j2", "Gd", None, "dimorphic_yeast"),
                                        ft("j3", "Gf", None, "filamentous_mycelium"), ft("j4", "Gn", None, "other_form")])
    got = {lt: one(con, lt)["yeast_origin"] for lt in ("A1", "A2", "A3", "A4")}
    assert got == {"A1": "Saccharomycotina", "A2": "Pezizomycotina-dimorphic", "A3": "filamentous", "A4": "other"}


def test_subphylum_phylum_fallback_order(tmp_path):
    cls = [("Kickxellomycetes", "Kickxellomycotina", "proposed", "")]
    phy = [("Kickxellomycota", "Kickxellomycotina", "proposed", ""), ("Glomeromycota", "Glomeromycotina", "proposed", "")]
    con = load(tmp_path, [("A1", "G", None, None, "Kickxellomycota"),          # no class: phylum table
                          ("A2", "G", "Kickxellomycetes", None, "Kickxellomycota"),   # class table wins
                          ("A3", "G", None, None, "Ascomycota"),               # phylum not listed: NULL
                          ("A4", "G", None, "Given", "Glomeromycota"),          # species value wins
                          ("A5", "G", None, None, None)],                      # no phylum: NULL
               cls_rows=cls, phy_rows=phy)
    got = {lt: (one(con, lt)["SUBPHYLUM"], one(con, lt)["subphylum_source"]) for lt in ("A1", "A2", "A3", "A4", "A5")}
    assert got == {"A1": ("Kickxellomycotina", "phylum_table"), "A2": ("Kickxellomycotina", "class_table"),
                   "A3": (None, None), "A4": ("Given", "species"), "A5": (None, None)}


def _labels_statement():
    sql = EXTRA_SQL.read_text()
    start = sql.index("CREATE OR REPLACE VIEW v_genome_labels AS")
    return sql[start:sql.index(";\n", start) + 1]


def test_v_genome_labels_one_row_per_genome_with_traits(tmp_path):
    write_refs(tmp_path,
               ft_rows=[ft("jrk2", "Mucor", "wrong", "x"), ft("jrk1", "Mucor", "saprotroph", "filamentous_mycelium")],
               cls_rows=[("Mucoromycetes", "Mucoromycotina", "proposed", "")])
    con = duckdb.connect()
    con.execute("""CREATE TABLE species (LOCUSTAG VARCHAR, ASMID VARCHAR, PHYLUM VARCHAR, CLASS VARCHAR, "ORDER" VARCHAR,
                   FAMILY VARCHAR, GENUS VARCHAR, SPECIES VARCHAR, STRAIN VARCHAR, NCBI_TAXONID VARCHAR, SUBPHYLUM VARCHAR)""")
    con.executemany("INSERT INTO species VALUES (?,?,?,?,?,?,?,?,?,?,?)", [
        ("A1", "G1", "Mucoromycota", "Mucoromycetes", "Mucorales", "Mucoraceae", "Mucor", "Mucor a", "s1", "1", None),
        ("A2", "G2", "Mucoromycota", "Mucoromycetes", "Mucorales", "Mucoraceae", "Mucor", "Mucor b", "s2", "2", None),
        ("A3", "G3", None, None, None, None, None, None, None, None, None)])
    con.execute("CREATE TABLE samples (LOCUSTAG VARCHAR, TRANSL_TABLE INTEGER)")
    con.execute("INSERT INTO samples VALUES ('A1',1),('A2',1),('A3',1)")
    con.execute("""CREATE TABLE asm_stats (ASMID VARCHAR, annotation_status VARCHAR, TOTAL_LENGTH BIGINT, contig_count INTEGER,
                   N50_bp BIGINT, GC_PERCENT DOUBLE, masked_pct DOUBLE)""")
    con.execute("INSERT INTO asm_stats VALUES ('G1','annotated',1,1,60000,50,1),('G2','annotated',1,1,60000,50,1),('G3','annotated',1,1,1,50,1)")
    con.execute("""CREATE TABLE busco_genome (ASMID VARCHAR, complete_pct DOUBLE, single_pct DOUBLE, duplicated_pct DOUBLE,
                   fragmented_pct DOUBLE, missing_pct DOUBLE)""")
    con.execute("INSERT INTO busco_genome VALUES ('G1',90,89,1,5,5),('G2',80,79,1,10,10),('G3',50,50,0,25,25)")
    apply_traits(con, tmp_path)
    con.execute(_labels_statement())
    assert con.execute("SELECT count(*) FROM v_genome_labels").fetchone()[0] == 3
    assert con.execute("SELECT count(DISTINCT LOCUSTAG) FROM v_genome_labels").fetchone()[0] == 3
    r = con.execute("""SELECT primary_lifestyle, growth_form, in_fungaltraits, SUBPHYLUM, subphylum_source
                       FROM v_genome_labels WHERE LOCUSTAG = 'A1'""").fetchone()
    assert r == ("saprotroph", "filamentous_mycelium", True, "Mucoromycotina", "class_table")
    r3 = con.execute("SELECT primary_lifestyle, in_fungaltraits, SUBPHYLUM FROM v_genome_labels WHERE LOCUSTAG = 'A3'").fetchone()
    assert r3 == (None, False, None)


# ---- bin/check_traits_vs_atlas.py ----
import subprocess
import sys

CHECK = BIN / "check_traits_vs_atlas.py"


def _check_fixture(tmp_path, db_life, dup_genus=False):
    db = tmp_path / "t.duckdb"
    con = duckdb.connect(str(db))
    con.execute("""CREATE TABLE v_genome_labels (LOCUSTAG VARCHAR, GENUS VARCHAR, PHYLUM VARCHAR, CLASS VARCHAR,
                   growth_form VARCHAR, primary_lifestyle VARCHAR, SUBPHYLUM VARCHAR)""")
    con.executemany("INSERT INTO v_genome_labels VALUES (?,?,?,?,?,?,?)", [
        ("A1", "Mucor", "Mucoromycota", "Mucoromycetes", "f", "saprotroph", None),
        ("A2", "Other", "Ascomycota", "Eurotiomycetes", "f", db_life, "Pezizomycotina")])
    con.close()
    atlas = tmp_path / "atlas.parquet"
    duckdb.execute(f"""COPY (SELECT * FROM (VALUES ('A1','f','saprotroph'),('A2','f','saprotroph'))
                        t(LOCUSTAG, growth_form, primary_lifestyle)) TO '{atlas}' (FORMAT parquet)""")
    ftcsv = tmp_path / "ft.csv"
    rows = [ft("j1", "Mucor", "saprotroph"), ft("j2", "Other", "saprotroph")]
    if dup_genus:
        rows.append(ft("j3", "Other", "different"))
    _write(ftcsv, FT_COLS, rows, ",")
    return db, atlas, ftcsv


def _run_check(db, atlas, ftcsv):
    return subprocess.run([sys.executable, str(CHECK), "--db", str(db), "--atlas", str(atlas), "--ft-csv", str(ftcsv)],
                          capture_output=True, text=True)


def test_check_script_flags_a_mismatch(tmp_path):
    r = _run_check(*_check_fixture(tmp_path, "parasite"))
    assert r.returncode == 1, r.stdout + r.stderr
    assert "A2" in r.stdout


def test_check_script_passes_when_equal(tmp_path):
    r = _run_check(*_check_fixture(tmp_path, "saprotroph"))
    assert r.returncode == 0, r.stdout + r.stderr
    assert "Mucoromycetes" in r.stdout          # unmapped class report: A1 has NULL SUBPHYLUM


def test_check_script_tolerates_mismatch_in_duplicate_genus(tmp_path):
    r = _run_check(*_check_fixture(tmp_path, "parasite", dup_genus=True))
    assert r.returncode == 0, r.stdout + r.stderr
    assert "duplicate" in r.stdout.lower()


def test_check_script_reports_but_tolerates_genome_set_difference(tmp_path):
    db, atlas, ftcsv = _check_fixture(tmp_path, "saprotroph")
    duckdb.execute(f"""COPY (SELECT * FROM (VALUES ('A1','f','saprotroph'),('A2','f','saprotroph'),('A9','f','x'))
                        t(LOCUSTAG, growth_form, primary_lifestyle)) TO '{atlas}' (FORMAT parquet)""")
    r = _run_check(db, atlas, ftcsv)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "only in atlas 1" in r.stdout
