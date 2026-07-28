#!/usr/bin/env python3
"""
pick_representative_strain.py — select representative strain per species from ANI TSV
and backfill shared ab-initio parameter store.

Ported from nextflow/bin/species_reuse_clusters.py to read from the merged ANI TSV
produced by ANI_REPRESENTATIVE_SELECT (instead of DuckDB). Writes the same
abinitio_reuse_assignments.csv format so FUNANNOTATE_PREDICTION's loadAbinitioReuseMap
consumes it without modification.

Usage:
    pick_representative_strain.py \
        --ani-tsv all_pairs_merged.tsv \
        --predict-input predict_input_for_ani.tsv \
        --samples samples.csv \
        --busco-dir results/genome_stats/BUSCO_genome \
        --out-dir genome_annotation/_reuse_assignments \
        --ani-threshold 99.0 \
        [--shared-root gene_prediction_shared_abinitio] \
        [--target genome_annotation] \
        [--augustus-config lib/augustus/3.5/config]
"""

import argparse
import csv
import json
import os
import re
import shutil
import sys
from collections import defaultdict
from datetime import date, datetime
from pathlib import Path
from typing import Optional

_ASMID_EXT_RE = re.compile(r'(\.fasta|\.fna|\.fa|\.faa|\.fas)(\.gz)?$')


def clean_strain(raw_strain: str) -> str:
    s = (raw_strain or "").strip()
    s = re.sub(r"['\"]", "", s)
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(raw_species: str, raw_strain: str) -> str:
    sp = re.sub(r"['\"]", "", (raw_species or "").strip())
    st = clean_strain(raw_strain)
    parts = [p for p in (sp, st) if p]
    tag = "_".join(parts)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def load_samples(samples_path: Path) -> dict:
    rows = {}
    with open(samples_path, newline="") as fh:
        for row in csv.DictReader(fh):
            asmid = (row.get("ASMID") or "").strip()
            species = (row.get("SPECIES") or "").strip().strip("'\"")
            strain = (row.get("STRAIN") or "").strip()
            if not asmid or not species:
                continue
            out = make_sample_tag(species, strain)
            species_tag = re.sub(r"\s+", "_", species)
            rows[asmid] = {
                "species": species, "strain": strain, "out": out,
                "species_tag": species_tag,
            }
    return rows


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
    return {"complete_pct": float(m[1]), "n50_bp": n50_bp}


def load_busco_by_out(busco_dir: Path) -> dict:
    result = {}
    if not busco_dir.is_dir():
        return result
    for p in busco_dir.glob("*.BUSCO_summary.*.txt"):
        out = p.name.split(".BUSCO_summary.")[0]
        parsed = parse_busco_genome_summary(p)
        if parsed:
            result[out] = parsed
    return result


def load_predict_input(path: Path) -> dict:
    """Return {asmid: {out, species, strain, ...}} from the predict_input TSV."""
    result = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            asmid = row.get("asmid", "").strip()
            if asmid:
                result[asmid] = dict(row)
    return result


def load_ani_pairs(ani_tsv_path: Path) -> dict:
    """Return {(q_asmid, r_asmid): ani} for every pair in the merged TSV.
    Infers asmid from filename stem (everything before the first '.').
    Handles both 3-column (q,r,ani) and 5-column (q,r,ani,af_r,af_q) formats."""
    pairs = {}
    if not ani_tsv_path.exists() or ani_tsv_path.name == '/dev/null':
        return pairs
    with open(ani_tsv_path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            return pairs
        for row in reader:
            if len(row) < 3:
                continue
            q_file = row[0].strip()
            r_file = row[1].strip()
            try:
                ani = float(row[2])
            except (ValueError, IndexError):
                continue
            # Strip common sequence file extensions to get the asmid
            q_asmid = _ASMID_EXT_RE.sub('', q_file)
            r_asmid = _ASMID_EXT_RE.sub('', r_file)
            if q_asmid and r_asmid:
                pairs[(q_asmid, r_asmid)] = ani
                pairs.setdefault((r_asmid, q_asmid), ani)
    return pairs


def _sparse_ani_warning(species: str, asmids: list, ani_pairs: dict,
                         rep_asmid: str) -> None:
    """Warn if ANI TSV appears sparse — some intra-species pairs are missing."""
    n = len(asmids)
    expected = n * (n - 1)  # (q,r) + (r,q) stored for every comparison
    asmid_set = set(asmids)
    found = sum(
        1 for (q, r) in ani_pairs
        if q in asmid_set and r in asmid_set and q != r
    )
    missing = expected - found
    if missing > 0:
        print(f"[WARN] {species}: ANI TSV has only {found}/{expected} intra-species "
              f"pairs ({missing} missing). This typically means a sparse ANI method "
              f"(skani/mash/sourmash) missed some comparisons. Strain {rep_asmid} "
              f"was selected as representative on BUSCO+N50 alone.",
              file=sys.stderr)


def pick_representative(asmids: list, busco_by_out: dict, samples: dict,
                         ani_pairs: dict) -> str:
    asmid_set = set(asmids)
    ani_covered = [
        a for a in asmids
        if any((a, other) in ani_pairs for other in asmid_set if other != a)
    ]
    candidates = ani_covered if ani_covered else asmids

    def sort_key(asmid):
        out = samples.get(asmid, {}).get("out", "")
        b = busco_by_out.get(out, {"complete_pct": -1.0, "n50_bp": -1})
        return (b["complete_pct"], b["n50_bp"], tuple(-ord(c) for c in out))

    return max(candidates, key=sort_key)


def compute_assignments(predict_input: dict, samples: dict, ani_pairs: dict,
                         busco_by_out: dict, threshold: float) -> list:
    by_species = defaultdict(list)
    for asmid, row in predict_input.items():
        sp = row.get("species", "").strip()
        if sp:
            by_species[sp].append(asmid)

    assignments = []
    for species, asmids in sorted(by_species.items()):
        if len(asmids) < 2:
            continue
        rep_asmid = pick_representative(asmids, busco_by_out, samples, ani_pairs)
        rep_out = samples.get(rep_asmid, {}).get("out", "")
        _sparse_ani_warning(species, asmids, ani_pairs, rep_asmid)

        for asmid in sorted(asmids, key=lambda a: samples.get(a, {}).get("out", "")):
            out = samples.get(asmid, {}).get("out", "")
            is_rep = (asmid == rep_asmid)
            if is_rep:
                assignments.append({
                    "species": species, "out": out,
                    "is_representative": True,
                    "representative_out": rep_out,
                    "ani_to_representative": 100.0,
                    "reuse_eligible": False,
                })
                continue
            key = (asmid, rep_asmid)
            ani = ani_pairs.get(key)
            eligible = ani is not None and ani >= threshold
            assignments.append({
                "species": species, "out": out,
                "is_representative": False,
                "representative_out": rep_out,
                "ani_to_representative": ani if ani is not None else "",
                "reuse_eligible": eligible,
            })
    return assignments


def write_assignments(assignments: list, out_dir: Path):
    """Write per-species abinitio_reuse_assignments.{species}.csv files and
    a combined abinitio_reuse_assignments.csv for loadAbinitioReuseMap."""
    out_dir.mkdir(parents=True, exist_ok=True)
    by_species = defaultdict(list)
    for row in assignments:
        by_species[row["species"]].append(row)
    
    fieldnames = ["species", "out", "is_representative", "representative_out",
                  "ani_to_representative", "reuse_eligible"]
    
    for species, rows in by_species.items():
        species_tag = re.sub(r"\s+", "_", species)
        out_path = out_dir / f"abinitio_reuse_assignments.{species_tag}.csv"
        with open(out_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            for row in rows:
                w.writerow(row)
    
    combined_path = out_dir / "abinitio_reuse_assignments.csv"
    with open(combined_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in assignments:
            w.writerow(row)


def find_rep_predict_dir(target: Path, out: str) -> Optional[Path]:
    d = target / out
    ab = d / "predict_misc" / "ab_initio_parameters"
    if ab.is_dir() and (d / "predict_results" / f"{out}.gbk").exists():
        return d
    if ab.is_dir() and (d / "predict_results" / f"{out}.gbk.gz").exists():
        return d
    return None


def backfill_species_store(species: str, rep_out: str, target: Path,
                            shared_root: Path, threshold: float,
                            aug_cfg: Optional[Path] = None,
                            dry_run: bool = False) -> bool:
    species_tag = re.sub(r"\s+", "_", species)
    augustus_name = species_tag.lower()
    rep_dir = find_rep_predict_dir(target, rep_out)
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
        print(f"[WARN] {species}: representative {rep_out} has no ab-initio "
              f"components to share (checked {ab_initio})", file=sys.stderr)
        return False

    if dry_run:
        print(f"[DRY-RUN] Would backfill {store_dir} from {rep_out}: "
              f"{sorted(components)}")
        return True

    store_dir.mkdir(parents=True, exist_ok=True)
    params_json_path = store_dir / "parameters.json"

    # Guard: if another project's pipeline already wrote a parameters.json from
    # the same (or a different) representative, don't stomp it.  This prevents
    # concurrent funannotate runs from races where project A removes project B's
    # store just before writing its own — a partial write window that could leave
    # project B's non-representatives with broken shared-params paths.
    if params_json_path.exists():
        existing = {}
        try:
            existing = json.loads(params_json_path.read_text())
        except (json.JSONDecodeError, OSError):
            pass  # corrupted — overwrite safely
        if existing and existing.get("codingquarry") != [{}]:
            print(f"[INFO] {store_dir} already has a parameters.json written by "
                  f"another run (representative={existing.get('representative_out', '?')}) "
                  f"— skipping backfill to avoid stomping it.", file=sys.stderr)
            return False

    params_json = {
        "augustus": [{}], "genemark": [{}], "snap": [{}],
        "codingquarry": [{}], "glimmerhmm": [{}], "table": 1,
    }

    if "augustus" in components:
        dest = store_dir / augustus_name
        if dest.exists():
            shutil.rmtree(dest)
        dest.mkdir(parents=True)
        for f in components["augustus"].iterdir():
            suffix = f.name[len(lower_out):] if f.name.startswith(lower_out) \
                else f"_{f.name}"
            dest_f = dest / f"{augustus_name}{suffix}"
            if f.suffix == ".cfg":
                dest_f.write_text(f.read_text(errors="replace")
                                   .replace(lower_out, augustus_name))
            else:
                shutil.copyfile(f, dest_f)
        params_json["augustus"] = [{
            "source": "ab-initio-reuse", "representative": rep_out,
            "path": str(dest.resolve()),
        }]
        if aug_cfg is not None:
            link = aug_cfg / "species" / augustus_name
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
        params_json["genemark"] = [{
            "source": "ab-initio-reuse", "representative": rep_out,
            "path": str(dest.resolve()),
        }]

    if "snap" in components:
        dest = store_dir / f"{species_tag}.snap.hmm"
        shutil.copyfile(components["snap"], dest)
        params_json["snap"] = [{
            "source": "ab-initio-reuse", "representative": rep_out,
            "path": str(dest.resolve()),
        }]

    glimmerhmm_stub = store_dir / "glimmerhmm_stub"
    glimmerhmm_stub.mkdir(exist_ok=True)
    params_json["glimmerhmm"] = [{
        "source": "suppressed-empty-stub",
        "path": str(glimmerhmm_stub.resolve()),
    }]

    (store_dir / "parameters.json").write_text(
        json.dumps(params_json, indent=2))
    provenance = {
        "representative_out": rep_out,
        "species": species,
        "augustus_species_name": augustus_name,
        "components": sorted(components),
        "ani_reuse_threshold": threshold,
        "date_captured": date.today().isoformat(),
        "generated_at": datetime.now().isoformat(),
        "source_predict_dir": str(rep_dir),
    }
    (store_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2))
    print(f"[INFO] Backfilled {store_dir} from representative {rep_out} "
          f"(components: {sorted(components)}, augustus_name={augustus_name})")
    return True


def write_repr_assignments(assignments: list, out_path: Path):
    """Write a simple TSV for pipeline channel joins:
    out, species, is_representative, representative_out,
    ani_to_representative, reuse_eligible"""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("out\tspecies\tis_representative\t"
                 "representative_out\tani_to_representative\treuse_eligible\n")
        for a in assignments:
            fh.write(f"{a['out']}\t{a['species']}\t"
                     f"{a['is_representative']}\t{a['representative_out']}\t"
                     f"{a['ani_to_representative']}\t{a['reuse_eligible']}\n")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ani-tsv", required=True,
                    help="Merged all-pairs TSV from CONCAT_ANI_TSVS (may be /dev/null)")
    ap.add_argument("--predict-input", required=True,
                    help="predict_input_for_ani.tsv written by ANI_REPRESENTATIVE_SELECT")
    ap.add_argument("--samples", required=True, help="samples.csv")
    ap.add_argument("--busco-dir", required=True,
                    help="results/genome_stats/BUSCO_genome/")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for abinitio_reuse_assignments.csv")
    ap.add_argument("--ani-threshold", type=float, default=99.0)
    ap.add_argument("--shared-root", default="")
    ap.add_argument("--target", default="")
    ap.add_argument("--augustus-config", default="")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    predict_input = load_predict_input(Path(args.predict_input))
    samples = load_samples(Path(args.samples))
    busco_by_out = load_busco_by_out(Path(args.busco_dir))
    ani_pairs = load_ani_pairs(Path(args.ani_tsv))

    print(f"[INFO] predict_input: {len(predict_input)} strains", file=sys.stderr)
    print(f"[INFO] samples: {len(samples)} entries", file=sys.stderr)
    print(f"[INFO] BUSCO summaries: {len(busco_by_out)}", file=sys.stderr)
    print(f"[INFO] ANI pairs: {len(ani_pairs)//2} unique pairs", file=sys.stderr)

    assignments = compute_assignments(predict_input, samples, ani_pairs,
                                       busco_by_out, args.ani_threshold)

    out_dir = Path(args.out_dir)
    write_assignments(assignments, out_dir)

    n_species = len({a["species"] for a in assignments})
    n_eligible = sum(1 for a in assignments if a["reuse_eligible"])
    n_rep = sum(1 for a in assignments if a["is_representative"])
    n_no_ani = sum(1 for a in assignments
                   if not a["is_representative"] and a["ani_to_representative"] == "")
    print(f"[INFO] Wrote {len(assignments)} rows ({n_species} species, "
          f"{n_rep} representative, {n_eligible} reuse_eligible, "
          f"{n_no_ani} fail-closed/missing-ANI) to {out_dir}/", file=sys.stderr)

    # Write repr_assignments.tsv for pipeline channel joins.
    repr_path = Path(args.out_dir) / "repr_assignments.tsv"
    write_repr_assignments(assignments, repr_path)

    if args.dry_run:
        return

    if args.shared_root:
        target = Path(args.target) if args.target else Path.cwd()
        shared_root = Path(args.shared_root)
        aug_cfg = Path(args.augustus_config) if args.augustus_config else None
        reps = {a["species"]: a["representative_out"]
                for a in assignments if a["is_representative"]}
        n_backfilled = 0
        for species, rep_out in sorted(reps.items()):
            if backfill_species_store(species, rep_out, target, shared_root,
                                       args.ani_threshold,
                                       aug_cfg=aug_cfg, dry_run=args.dry_run):
                n_backfilled += 1
        print(f"[INFO] Backfilled {n_backfilled}/{len(reps)} species stores",
              file=sys.stderr)


if __name__ == "__main__":
    main()
