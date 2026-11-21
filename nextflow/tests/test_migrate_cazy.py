"""Fixture tests for scripts/one-off/migrate_cazy_to_hash_buckets.py.

Builds a miniature results/function/cazy/ tree covering every category the
migration can meet, then checks the plan, the move, and that re-running is safe.
No HPCC paths and no real data.
"""
import gzip
import importlib.util
import os
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
MIG = REPO / "scripts" / "one-off" / "migrate_cazy_to_hash_buckets.py"


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


mig = _load(MIG, "migrate_cazy")
TYPES = mig.TYPES

HEADER = "Gene_ID\tEC\tcazyme_fam\tsub_fam\tdiamond_fam\tSubstrate\t#ofTools\n"

SAMPLES_HEADER = ("ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,"
                  "SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,TRANSL_TABLE,LOCUSTAG\n")


def write_genome_dir(cazy, name, locustag, empty=False, header_only=False):
    d = cazy / name
    d.mkdir(parents=True, exist_ok=True)
    if empty:
        return d
    for t in TYPES:
        body = HEADER if header_only else HEADER + "%s_000001\t-\tGH5\t-\t-\t-\t1\n" % locustag
        with gzip.open(d / ("%s.%s.tsv.gz" % (name, t)), "wt") as fh:
            fh.write(body)
    return d


@pytest.fixture
def tree(tmp_path):
    """A runs root with one of every category."""
    root = tmp_path / "runs"
    cazy = root / "results" / "function" / "cazy"
    cazy.mkdir(parents=True)

    rows = [
        # (species, strain, locustag)
        ("Macrocybe indica", "ADT1", "FBE0FD36"),   # renamed: old dir is the loser
        ("Aspergillus nidulans", "FGSC_A4", "FAAAAAA1"),  # plain migrate
    ]
    with open(root / "samples.csv", "w") as fh:
        fh.write(SAMPLES_HEADER)
        for i, (sp, st, lt) in enumerate(rows):
            fh.write("GCA_%08d.1_x,%s,%s,PRJ,1,dikarya,Asco,Pezizo,Cls,,Ord,Fam,Gen,%s,1,%s\n"
                     % (i, sp, st, sp, lt))

    write_genome_dir(cazy, "Aspergillus_nidulans_FGSC_A4", "FAAAAAA1")   # MIGRATE
    write_genome_dir(cazy, "Macrocybe_indica_ADT1", "FBE0FD36")          # MIGRATE (winner)
    write_genome_dir(cazy, "Calocybe_indica_ADT1", "FBE0FD36")           # DUPLICATE_LOSER
    write_genome_dir(cazy, "Paracoccidioides_lutzii_Pb01", "FE07DE89")   # PHANTOM
    write_genome_dir(cazy, "Saccharomyces_bayanus_FM677", "", empty=True)       # EMPTY
    write_genome_dir(cazy, "Weird_genome_X1", "FZZZZZZ9", header_only=True)     # UNRESOLVED
    return root, cazy


def plan_for(root, cazy):
    out_by_lt = mig.load_samples(str(root / "samples.csv"))
    return {r[0]: r for r in mig.build_plan(str(cazy), out_by_lt)}, out_by_lt


def test_every_category_is_classified(tree):
    root, cazy = tree
    plan, _ = plan_for(root, cazy)
    assert plan["Aspergillus_nidulans_FGSC_A4"][2] == "MIGRATE"
    assert plan["Macrocybe_indica_ADT1"][2] == "MIGRATE"
    assert plan["Calocybe_indica_ADT1"][2] == "DUPLICATE_LOSER"
    assert plan["Paracoccidioides_lutzii_Pb01"][2] == "PHANTOM"
    assert plan["Saccharomyces_bayanus_FM677"][2] == "EMPTY"
    assert plan["Weird_genome_X1"][2] == "UNRESOLVED"


def test_duplicate_winner_is_the_current_out_name(tree):
    """The renamed genome keeps the dir matching samples.csv, not the older name."""
    root, cazy = tree
    plan, _ = plan_for(root, cazy)
    assert plan["Calocybe_indica_ADT1"][4].endswith("Macrocybe_indica_ADT1")


def test_bucket_matches_the_shared_helper(tree):
    """The migration must bucket exactly as genome_stats_paths.py does, or every
    path it writes is one Nextflow will not look in."""
    root, cazy = tree
    plan, _ = plan_for(root, cazy)
    for name, lt in [("Aspergillus_nidulans_FGSC_A4", "FAAAAAA1"), ("Macrocybe_indica_ADT1", "FBE0FD36")]:
        assert plan[name][3] == mig.hash_bucket_for_type("cazy", lt)
        assert len(plan[name][3]) == 2


def test_apply_moves_and_quarantines(tree):
    root, cazy = tree
    plan_map, _ = plan_for(root, cazy)
    plan = list(plan_map.values())
    quarantine = root / "results" / "function" / "cazy.superseded_test"
    moved, quarantined, _skipped = mig.apply_plan(plan, str(cazy), str(quarantine))

    assert moved == 2
    assert quarantined == 4

    lt, bucket = "FBE0FD36", mig.hash_bucket_for_type("cazy", "FBE0FD36")
    for t in TYPES:
        assert (cazy / bucket / ("%s.%s.tsv.gz" % (lt, t))).exists()
    # the winner's source dir is gone, the loser is parked out of the glob's reach
    assert not (cazy / "Macrocybe_indica_ADT1").exists()
    assert (quarantine / "Calocybe_indica_ADT1").exists()
    assert (quarantine / "Paracoccidioides_lutzii_Pb01").exists()
    assert (quarantine / "Weird_genome_X1").exists()


def test_verify_passes_after_apply_and_rerun_is_a_noop(tree):
    root, cazy = tree
    plan_map, out_by_lt = plan_for(root, cazy)
    plan = list(plan_map.values())
    quarantine = root / "results" / "function" / "cazy.superseded_test"
    mig.apply_plan(plan, str(cazy), str(quarantine))

    missing, survivors = mig.verify(plan, str(cazy))
    assert missing == []
    assert survivors == []

    # Re-running rebuilds the plan from disk: the bucket dirs must not be
    # mistaken for genome dirs, so there is nothing left to do.
    plan2 = mig.build_plan(str(cazy), out_by_lt)
    assert plan2 == []
    moved2, quarantined2, _ = mig.apply_plan(plan, str(cazy), str(quarantine))
    assert (moved2, quarantined2) == (0, 0)


def test_unresolved_is_never_guessed(tree):
    """A directory whose LOCUSTAG cannot be read must not be assigned one."""
    root, cazy = tree
    plan, _ = plan_for(root, cazy)
    assert plan["Weird_genome_X1"][1] == ""
    assert plan["Weird_genome_X1"][3] == ""


def test_locustag_read_from_file_contents(tmp_path):
    p = tmp_path / "x.overview.tsv.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(HEADER + "F1234567_000123\t-\tGH3\t-\t-\t-\t2\n")
    assert mig.locustag_from_overview(str(p)) == "F1234567"
    empty = tmp_path / "e.overview.tsv.gz"
    with gzip.open(empty, "wt") as fh:
        fh.write(HEADER)
    assert mig.locustag_from_overview(str(empty)) is None
    assert mig.locustag_from_overview(str(tmp_path / "nope.gz")) is None
