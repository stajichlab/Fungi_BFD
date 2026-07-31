#!/usr/bin/env python3
"""
species_reuse_clusters.py — compute ANI-gated ab-initio parameter reuse clusters
and backfill the shared-parameter store from already-completed PREDICT runs.

Design: todo/species_level_abinitio_reuse.md

For each species with >=2 strains, picks one representative strain (highest BUSCO
genome-mode completeness, Scaffold N50 tiebreak, alphabetical `out` final tiebreak),
then marks every other strain of that species `reuse_eligible` iff its ANI to the
representative (from ani.duckdb) is >= --ani-threshold (default 99.0). Missing ANI data
is treated as ineligible (fail-closed), never as eligible.

Writes abinitio_reuse_assignments.csv (species, out, is_representative,
representative_out, ani_to_representative, reuse_eligible) and, for every species
whose representative already has a completed PREDICT (predict_results/<out>.gbk[.gz]
+ predict_misc/ab_initio_parameters/), backfills
<shared-root>/<species_tag>/{parameters.json, provenance.json,
<species_tag_lower>/ (AUGUSTUS species dir), <species_tag>.genemark.mod,
<species_tag>.snap.hmm} by copying (not symlinking) from the representative's own
ab_initio_parameters/. shared-root defaults to a top-level sibling of
genome_annotation/ (gene_prediction_shared_abinitio/), not nested inside it.

AUGUSTUS requires a species directory's basename to exactly match its parameter
files' prefix (species/anidulans/anidulans_parameters.cfg) -- the copied files are
renamed from the representative's own name to the lowercased species tag to satisfy
this, and the resulting directory is also symlinked into --augustus-config's own
species/ directory so it's discoverable there too, without duplicating the shared
cgp/extrinsic/model/parameters/profile/ subdirectories per species.

Usage:
    python3 species_reuse_clusters.py \\
        --samples samples.csv \\
        --ani-db results/ANI/skani/GENUS/ani.duckdb \\
        --busco-genome-dir results/genome_stats/BUSCO_genome \\
        --target genome_annotation \\
        --out abinitio_reuse_assignments.csv \\
        [--shared-root gene_prediction_shared_abinitio] \\
        [--augustus-config lib/augustus/3.5/config] \\
        [--ani-threshold 99.0] [--species "Aspergillus fumigatus"] [--no-backfill]
"""

import argparse
import csv
import json
import re
import shutil
import duckdb
import sys
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path


# ── SampleUtils.groovy port (nextflow/lib/SampleUtils.groovy) ─────────────────
# Keep this in sync with SampleUtils.groovy's cleanStrain/makeSampleTag — this
# is the tag used everywhere in funannotate.nf (`out`), so a Python script that
# recomputes it differently will silently mismatch every downstream directory.

def clean_strain(raw_strain: str) -> str:
    s = (raw_strain or "").strip()
    s = re.sub(r"""['"]""", "", s)
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species: str, raw_strain: str) -> str:
    sp = re.sub(r"""['"]""", "", (raw_species or "").strip())
    st = clean_strain(raw_strain)
    parts = [p for p in (sp, st) if p]
    tag = "_".join(parts)
    tag = re.sub(r"[\s/#\[\]?{}]+", "_", tag)
    return tag


# ── BUSCO genome-mode summary parsing ──────────────────────────────────────────

_COMPLETE_RE = re.compile(r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)")
_N50_RE = re.compile(r"([\d.]+)\s*(Mbp|kbp|bp)\s*Scaffold N50", re.IGNORECASE)
_UNIT_MULT = {"bp": 1, "kbp": 1_000, "Mbp": 1_000_000}


def parse_busco_genome_summary(path: Path):
    text = path.read_text(errors="replace")
    m = _COMPLETE_RE.search(text)
    if not m:
        return None
    n50_bp = 0
    mn = _N50_RE.search(text)
    if mn:
        n50_bp = int(float(mn.group(1)) * _UNIT_MULT[mn.group(2)])
    return {"complete_pct": float(m.group(1)), "n50_bp": n50_bp}


def load_busco_genome_by_asmid(busco_dir: Path) -> dict:
    """Return {asmid: {complete_pct, n50_bp}} keyed by the BUSCO_genome filename stem,
    read directly rather than via the aggregated tables/BUSCO.csv.gz — that table's
    ASMID-basename join has real gaps (verified: only 1/78 Aspergillus fumigatus
    genomes survived the join), so this script reads the per-genome files directly
    as the design (todo/species_level_abinitio_reuse.md S4.1) specifies.

    BUSCO_genome is ASMID-keyed (assembly-level statistic) and hash-bucketed
    (nextflow/modules/BFD/BUSCO_GENOME/main.nf writes
    <bucket>/<asmid>.BUSCO_summary.<lineage>.txt) -- the `*/*` glob covers the
    one bucket-subdirectory level; keyed by ASMID now instead of the old
    Genus_species_strain `out` tag (genome_stats_storage_reorg.md, T-014)."""
    result = {}
    if not busco_dir.is_dir():
        return result
    for p in busco_dir.glob("*/*.BUSCO_summary.*.txt"):
        asmid = p.name.split(".BUSCO_summary.")[0]
        parsed = parse_busco_genome_summary(p)
        if parsed:
            result[asmid] = parsed
    return result


# ── samples.csv ────────────────────────────────────────────────────────────────

def load_samples(samples_path: Path) -> dict:
    """Return {asmid: {species, strain, out, species_tag}}."""
    rows = {}
    with open(samples_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = (row.get("ASMID") or "").strip()
            species = (row.get("SPECIES") or "").strip().strip("'\"")
            strain = (row.get("STRAIN") or "").strip()
            if not asmid or not species:
                continue
            out = make_sample_tag(species, strain)
            species_tag = re.sub(r"\s+", "_", species)
            rows[asmid] = {"species": species, "strain": strain, "out": out, "species_tag": species_tag}
    return rows


# ── ani.duckdb ──────────────────────────────────────────────────────────────────

def load_within_species_ani(ani_db_path: Path) -> dict:
    """Return {(asmid_a, asmid_b): ani} for every same-species pair, both
    (a,b) and (b,a) directions stored so lookups don't need to guess order.
    Asserts non-blank species columns on load (S3.3 of the plan: a known bug
    could silently leave these blank -- fail loudly rather than misclassify
    every strain as ineligible)."""
    con = duckdb.connect(str(ani_db_path))
    result = con.execute(
        "SELECT COUNT(*), SUM(CASE WHEN query_species IS NULL OR query_species='' THEN 1 ELSE 0 END) "
        "FROM ani_pairs"
    ).fetchone()
    if result is None:
        total, blank = 0, 0
    else:
        total = result[0]
        blank = result[1] if result[1] is not None else 0
    if total and blank and blank > 0:
        sys.exit(
            f"ERROR: ani.db has {blank}/{total} rows with a blank query_species -- "
            "this is the exact failure mode from the 2026-07-19 names-glob bug "
            "(.living/learnings.md). Refusing to compute reuse clusters against "
            "unreliable species columns. Re-run compare_ANI.nf/COMBINE_ANI_TABLE first."
        )
    pairs = {}
    rows = con.execute(
        "SELECT query_asmid, ref_asmid, ani FROM ani_pairs WHERE query_species = ref_species"
    ).fetchall()
    for qa, ra, ani in rows:
        if not qa or not ra:
            continue
        pairs[(qa, ra)] = ani
        pairs.setdefault((ra, qa), ani)
    con.close()
    return pairs


def ani_between(pairs: dict, a: str, b: str):
    if a == b:
        return 100.0
    return pairs.get((a, b))


# ── Representative selection + reuse eligibility ────────────────────────────────

def pick_representative(asmids: list, busco_by_asmid: dict, samples: dict, ani_pairs: dict) -> str:
    """Highest BUSCO completeness, N50 tiebreak, alphabetical `out` final tiebreak
    (todo/species_level_abinitio_reuse.md S6.2) -- but ONLY among candidates that
    actually have at least one ANI pair to another strain in this species group.
    A representative with zero ani.db coverage (e.g. a genome added/BUSCO-scored
    after the last compare_ANI.nf run) would win on BUSCO score alone and orphan
    every other strain's reuse eligibility, even when those other strains have
    excellent mutual ANI coverage -- discovered by dry-running against Aspergillus
    fumigatus, where the top-BUSCO strain (F2, 99.3%) had 0 ani.db rows. Falls back
    to the full candidate set (still ranked by BUSCO/N50) only if NONE of them have
    any ANI coverage at all -- that's the genuinely-uncovered-species case (e.g.
    Beauveria bassiana before its SPECIES-level ANI run), where every strain will
    end up fail-closed regardless of which one is nominally "the representative"."""
    asmid_set = set(asmids)
    ani_covered = [a for a in asmids if any((a, other) in ani_pairs for other in asmid_set if other != a)]
    candidates = ani_covered if ani_covered else asmids

    def sort_key(asmid):
        out = samples[asmid]["out"]
        b = busco_by_asmid.get(asmid, {"complete_pct": -1.0, "n50_bp": -1})
        # max() picks the largest tuple; negate `out` alphabetical rank so that,
        # among true ties, the lexicographically SMALLEST `out` wins (a stable,
        # deterministic tiebreak per todo/species_level_abinitio_reuse.md S6.2).
        return (b["complete_pct"], b["n50_bp"], tuple(-ord(c) for c in out))
    return max(candidates, key=sort_key)


def compute_assignments(samples: dict, ani_pairs: dict, busco_by_asmid: dict, ani_threshold: float,
                         species_filter: str = None) -> list:
    by_species = defaultdict(list)
    for asmid, row in samples.items():
        if species_filter and row["species"] != species_filter:
            continue
        by_species[row["species"]].append(asmid)

    assignments = []
    for species, asmids in sorted(by_species.items()):
        if len(asmids) < 2:
            continue
        rep_asmid = pick_representative(asmids, busco_by_asmid, samples, ani_pairs)
        rep_out = samples[rep_asmid]["out"]
        for asmid in sorted(asmids, key=lambda a: samples[a]["out"]):
            out = samples[asmid]["out"]
            is_rep = asmid == rep_asmid
            if is_rep:
                assignments.append({
                    "species": species, "out": out, "is_representative": True,
                    "representative_out": rep_out, "ani_to_representative": 100.0,
                    "reuse_eligible": False,  # the representative doesn't "reuse" itself
                })
                continue
            ani = ani_between(ani_pairs, asmid, rep_asmid)
            # Fail-closed: sparse ANI methods omit low-identity pairs entirely, so
            # "no pair found" (ani is None) must be treated identically to "below
            # threshold" -- never as eligible (todo/species_level_abinitio_reuse.md S4.1).
            eligible = ani is not None and ani >= ani_threshold
            assignments.append({
                "species": species, "out": out, "is_representative": False,
                "representative_out": rep_out,
                "ani_to_representative": ani if ani is not None else "",
                "reuse_eligible": eligible,
            })
    return assignments


def write_assignments(assignments: list, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "species", "out", "is_representative", "representative_out",
            "ani_to_representative", "reuse_eligible",
        ])
        w.writeheader()
        for row in assignments:
            w.writerow(row)


# ── Backfill: copy representative's ab_initio_parameters into a species-level store ──

def find_legacy_predict_dir(target: Path, out: str) -> Path:
    d = target / out
    gbk = d / "predict_results" / f"{out}.gbk"
    gbk_gz = d / "predict_results" / f"{out}.gbk.gz"
    ab_initio = d / "predict_misc" / "ab_initio_parameters"
    if (gbk.exists() or gbk_gz.exists()) and ab_initio.is_dir():
        return d
    return None


def backfill_species_store(species: str, rep_out: str, target: Path, shared_root: Path,
                            ani_threshold: float, augustus_config: Path = None,
                            dry_run: bool = False) -> bool:
    """Copy rep_out's ab_initio_parameters into shared_root/<species_tag>/.

    AUGUSTUS requires a species directory's basename to exactly match the prefix
    of the parameter files inside it (e.g. species/anidulans/anidulans_parameters.cfg)
    -- confirmed by comparing a stock species against an early version of this store
    that kept the representative's own strain-cased file names
    (aspergillus_fumigatus_z5_*.cfg) inside a directory renamed to the species tag
    (Aspergillus_fumigatus/), which broke AUGUSTUS's lookup entirely ("augustus
    --proteinprofile test failed", confirmed during T-004's first real validation
    run). Fixed by renaming the files themselves to <augustus_name>_* on copy, where
    augustus_name = species_tag.lower() (matching AUGUSTUS/funannotate's existing
    all-lowercase species-naming convention, e.g. "anidulans",
    "aspergillus_fumigatus_z5").

    Also symlinks the species dir into augustus_config's own species/ directory
    (the AUGUSTUS_CONFIG_PATH every funannotate predict call already uses) so it's
    discoverable via the simpler --augustus_species <name> flag too, not just via
    the absolute path embedded in parameters.json -- avoids needing to trust that
    funannotate resolves an arbitrary out-of-tree absolute augustus path correctly,
    and reuses the config tree's existing cgp/extrinsic/model/parameters/profile/
    subdirectories rather than duplicating them per species.

    Returns True if a store was written (or would be, under --dry-run)."""
    species_tag = re.sub(r"\s+", "_", species)         # e.g. "Aspergillus_fumigatus"
    augustus_name = species_tag.lower()                 # e.g. "aspergillus_fumigatus"
    rep_dir = find_legacy_predict_dir(target, rep_out)
    if rep_dir is None:
        return False

    ab_initio = rep_dir / "predict_misc" / "ab_initio_parameters"
    lower_out = rep_out.lower()
    store_dir = shared_root / species_tag

    components = {}
    aug_src = ab_initio / "augustus" / "species" / lower_out
    genemark_src = ab_initio / f"{lower_out}.genemark.mod"
    snap_src = ab_initio / f"{lower_out}.snap.hmm"
    if aug_src.is_dir():
        components["augustus"] = aug_src
    if genemark_src.is_file():
        components["genemark"] = genemark_src
    if snap_src.is_file():
        components["snap"] = snap_src

    if not components:
        print(f"[WARN] {species}: representative {rep_out} has no ab-initio components "
              f"to share (checked {ab_initio})", file=sys.stderr)
        return False

    if dry_run:
        print(f"[DRY-RUN] Would backfill {store_dir} from {rep_out}: {sorted(components)} "
              f"(augustus species dir would be named '{augustus_name}')")
        return True

    store_dir.mkdir(parents=True, exist_ok=True)
    params_json = {"augustus": [{}], "genemark": [{}], "snap": [{}], "codingquarry": [{}],
                    "glimmerhmm": [{}], "table": 1}

    # Paths in parameters.json must be absolute: funannotate predict -p <json> is
    # invoked from arbitrary Nextflow work directories, not this script's cwd (real
    # funannotate-produced parameters.json files always use absolute paths too).
    if "augustus" in components:
        dest = store_dir / augustus_name
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True, exist_ok=True)
        # Rename each file's prefix from the representative's own name to
        # augustus_name so it matches this directory's basename (the AUGUSTUS
        # convention this whole function exists to satisfy -- see docstring).
        # *_parameters.cfg additionally hardcodes sibling filenames INSIDE its own
        # content (/BaseCount/weightMatrixFile, /ExonModel/infile, etc. -- the file
        # itself says "change this to your species if at all necessary") -- renaming
        # the file on disk alone leaves those internal references pointing at the
        # representative's old name, which AUGUSTUS then fails to open. Confirmed by
        # a real predict run: it crash-looped every genome chunk on "Couldn't open
        # the file with the weight matrix: <representative>_weightmatrix.txt" despite
        # the file being correctly renamed on disk. Text-substitute the old prefix
        # for the new one in .cfg file content, not just the filename.
        for f in components["augustus"].iterdir():
            suffix = f.name[len(lower_out):] if f.name.startswith(lower_out) else f"_{f.name}"
            dest_f = dest / f"{augustus_name}{suffix}"
            if f.suffix == ".cfg":
                text = f.read_text(errors="replace").replace(lower_out, augustus_name)
                dest_f.write_text(text)
            else:
                shutil.copyfile(f, dest_f)
        params_json["augustus"] = [{"source": "ab-initio-reuse", "representative": rep_out,
                                     "path": str(dest.resolve())}]

        if augustus_config is not None:
            link = augustus_config / "species" / augustus_name
            link.parent.mkdir(parents=True, exist_ok=True)
            if link.is_symlink() or link.exists():
                if link.is_symlink() or link.is_file():
                    link.unlink()
                else:
                    shutil.rmtree(link)
            link.symlink_to(dest.resolve())
            print(f"[INFO] Symlinked {link} -> {dest.resolve()}")

    if "genemark" in components:
        dest = store_dir / f"{species_tag}.genemark.mod"
        shutil.copyfile(components["genemark"], dest)
        params_json["genemark"] = [{"source": "ab-initio-reuse", "representative": rep_out,
                                     "path": str(dest.resolve())}]
    if "snap" in components:
        dest = store_dir / f"{species_tag}.snap.hmm"
        shutil.copyfile(components["snap"], dest)
        params_json["snap"] = [{"source": "ab-initio-reuse", "representative": rep_out,
                                 "path": str(dest.resolve())}]

    # glimmerhmm has weight 0 in FUNANNOTATE_PREDICT (-w glimmerhmm:0) -- it never
    # contributes to the EVM consensus -- but funannotate predict.py still needs to
    # decide a *training method* for it, and without a "path" key it falls back to
    # "busco", which triggers a full genome-mode BUSCO run (funannotate's RunBusco flag
    # is set if ANY predictor needs "busco" mode, regardless of that predictor's weight).
    # In normal TRAIN+PREDICT runs this is masked because PASA data gives glimmerhmm a
    # cheap "pasa" fallback instead -- but a no-TRAIN run (species_level_abinitio_reuse.md
    # T-013) has no PASA data at all, so it silently pays for the expensive BUSCO step
    # just to train a predictor whose output is discarded. An empty stub directory is
    # enough: it registers RunModes["glimmerhmm"]="pretrained" (predict.py only checks
    # for the "path" key's existence, not real trained content) and glimmerhmm's own
    # execution block is separately gated on weight > 0, so with weight 0 it does
    # nothing at all. See .living/learnings.md 2026-07-25 entry.
    glimmerhmm_stub = store_dir / "glimmerhmm_stub"
    glimmerhmm_stub.mkdir(parents=True, exist_ok=True)
    params_json["glimmerhmm"] = [{"source": "suppressed-empty-stub",
                                   "path": str(glimmerhmm_stub.resolve())}]

    (store_dir / "parameters.json").write_text(json.dumps(params_json, indent=2))
    provenance = {
        "representative_out": rep_out,
        "species": species,
        "augustus_species_name": augustus_name,
        "components": sorted(components),
        "ani_reuse_threshold": ani_threshold,
        "date_captured": date.today().isoformat(),
        "generated_at": datetime.now().isoformat(),
        "source_predict_dir": str(rep_dir),
    }
    (store_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))
    print(f"[INFO] Backfilled {store_dir} from representative {rep_out} "
          f"(components: {sorted(components)}, augustus_name={augustus_name})")
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", required=True, help="samples.csv path")
    ap.add_argument("--ani-db", required=True, help="ani.duckdb DuckDB path (compare_ANI.nf output)")
    ap.add_argument("--busco-genome-dir", required=True,
                     help="BUSCO genome-mode results dir (results/genome_stats/BUSCO_genome)")
    ap.add_argument("--target", required=True,
                     help="funannotate PREDICT output root (params.target, e.g. genome_annotation)")
    ap.add_argument("--out", required=True, help="Output abinitio_reuse_assignments.csv path")
    ap.add_argument("--shared-root",
                     default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/gene_prediction_shared_abinitio",
                     help="Shared ab-initio store root -- a top-level sibling of genome_annotation/, "
                          "not nested inside it (default: gene_prediction_shared_abinitio/ at the project root)")
    ap.add_argument("--augustus-config",
                     default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/lib/augustus/3.5/config",
                     help="AUGUSTUS_CONFIG_PATH (params.augustus_config) -- each backfilled species gets "
                          "symlinked into <this>/species/<augustus_name>/ so it's discoverable the same way "
                          "any other AUGUSTUS species is, not just via the absolute path in parameters.json. "
                          "Pass empty string to skip symlinking.")
    ap.add_argument("--ani-threshold", type=float, default=99.0)
    ap.add_argument("--species", default=None, help="Restrict to one species (exact SPECIES match)")
    ap.add_argument("--no-backfill", action="store_true",
                     help="Only write the assignments CSV; skip backfilling the shared store")
    ap.add_argument("--dry-run", action="store_true",
                     help="Print what backfill would do without writing files")
    args = ap.parse_args()

    samples = load_samples(Path(args.samples))
    print(f"[INFO] Loaded {len(samples)} samples.csv rows", file=sys.stderr)

    busco_by_asmid = load_busco_genome_by_asmid(Path(args.busco_genome_dir))
    print(f"[INFO] Loaded {len(busco_by_asmid)} BUSCO genome-mode summaries", file=sys.stderr)

    ani_pairs = load_within_species_ani(Path(args.ani_db))
    print(f"[INFO] Loaded {len(ani_pairs) // 2} unique within-species ANI pairs", file=sys.stderr)

    assignments = compute_assignments(samples, ani_pairs, busco_by_asmid, args.ani_threshold,
                                       species_filter=args.species)
    write_assignments(assignments, Path(args.out))
    n_species = len({a["species"] for a in assignments})
    n_eligible = sum(1 for a in assignments if a["reuse_eligible"])
    n_no_ani = sum(1 for a in assignments if not a["is_representative"] and a["ani_to_representative"] == "")
    print(f"[INFO] Wrote {len(assignments)} assignment rows across {n_species} species "
          f"to {args.out} ({n_eligible} reuse_eligible, {n_no_ani} missing-ANI/fail-closed)",
          file=sys.stderr)

    if args.no_backfill:
        return

    target = Path(args.target)
    shared_root = Path(args.shared_root)
    augustus_config = Path(args.augustus_config) if args.augustus_config else None
    reps = {a["species"]: a["representative_out"] for a in assignments}
    n_backfilled = 0
    for species, rep_out in sorted(reps.items()):
        if backfill_species_store(species, rep_out, target, shared_root, args.ani_threshold,
                                   augustus_config=augustus_config, dry_run=args.dry_run):
            n_backfilled += 1
    print(f"[INFO] Backfilled {n_backfilled}/{len(reps)} species stores", file=sys.stderr)


if __name__ == "__main__":
    main()
