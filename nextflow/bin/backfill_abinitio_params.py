#!/usr/bin/env python3
"""
backfill_abinitio_params.py — copy a representative strain's trained
AUGUSTUS/GeneMark/SNAP ab-initio parameters into the shared per-species
parameter store, so sibling strains can reuse them via
`funannotate predict -p parameters.json`.

One shared implementation, two callers (see .living/decisions.md):
  - pick_representative_strain.py (compare_ani pipeline): backfills right
    after representative selection, if that representative already has a
    prediction on disk from a prior run.
  - BACKFILL_ABINITIO_PARAMS (Nextflow process, nextflow/modules/funannotate/
    predict/BACKFILL_ABINITIO_PARAMS/main.nf): wired with a real channel
    dependency after a representative's FUNANNOTATE_PREDICT in the
    funannotate pipeline's --predict_scope all mode, and reused by the
    standalone --pipeline backfill_abinitio sweep.

Usage as a library:
    from backfill_abinitio_params import backfill_species_store, BackfillOutcome
    result = backfill_species_store(species, rep_out, target, shared_root,
                                     threshold, aug_cfg)
    if result.outcome is BackfillOutcome.REPRESENTATIVE_NOT_PREDICTED:
        ...

Usage as a CLI:
    backfill_abinitio_params.py --species "Neurospora crassa" \
        --representative-out Neurospora_crassa_OR74A \
        --target genome_annotation \
        --shared-root gene_prediction_shared_abinitio \
        --ani-threshold 99.0 \
        [--augustus-config lib/augustus/3.5/config]
"""

import argparse
import hashlib
import json
import os
import re
import shutil
import sys
from dataclasses import dataclass
from datetime import date, datetime
from enum import Enum
from pathlib import Path
from typing import Optional


class BackfillOutcome(str, Enum):
    REPRESENTATIVE_NOT_PREDICTED = "representative_not_predicted"
    UP_TO_DATE = "up_to_date"
    BACKFILLED = "backfilled"


@dataclass
class BackfillResult:
    outcome: BackfillOutcome
    store_dir: Optional[Path] = None
    components: Optional[list] = None


def find_rep_predict_dir(target: Path, out: str) -> Optional[Path]:
    d = target / out
    ab = d / "predict_misc" / "ab_initio_parameters"
    if ab.is_dir() and (d / "predict_results" / f"{out}.gbk").exists():
        return d
    if ab.is_dir() and (d / "predict_results" / f"{out}.gbk.gz").exists():
        return d
    return None


def _tree_hash(root: Path) -> str:
    """Stable content hash of every regular file under root, keyed by relative path."""
    h = hashlib.sha256()
    if not root.exists():
        return h.hexdigest()  # hash of "nothing" -- consistent empty-input sentinel
    for p in sorted(root.rglob("*")):
        if p.is_file():
            h.update(p.relative_to(root).as_posix().encode())
            h.update(p.read_bytes())
    return h.hexdigest()


def _store_content_hash(store_dir: Path, augustus_name: str, species_tag: str) -> str:
    """Hash of the parts of a store that reflect actually-trained content.
    Deliberately excludes parameters.json/provenance.json: those carry
    generated_at timestamps and absolute paths that always differ between
    two backfills of identical training output, which would defeat the
    whole point of comparing content before deciding to write."""
    h = hashlib.sha256()
    aug_dir = store_dir / augustus_name
    if aug_dir.is_dir():
        h.update(_tree_hash(aug_dir).encode())
    for name in (f"{species_tag}.genemark.mod", f"{species_tag}.snap.hmm"):
        f = store_dir / name
        if f.is_file():
            h.update(name.encode())
            h.update(f.read_bytes())
    return h.hexdigest()


def backfill_species_store(species: str, rep_out: str, target: Path,
                            shared_root: Path, threshold: float,
                            aug_cfg: Optional[Path] = None,
                            dry_run: bool = False) -> BackfillResult:
    species_tag = re.sub(r"\s+", "_", species)
    augustus_name = species_tag.lower()
    rep_dir = find_rep_predict_dir(target, rep_out)
    if rep_dir is None:
        return BackfillResult(BackfillOutcome.REPRESENTATIVE_NOT_PREDICTED)

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
        return BackfillResult(BackfillOutcome.REPRESENTATIVE_NOT_PREDICTED)

    if dry_run:
        print(f"[DRY-RUN] Would backfill {store_dir} from {rep_out}: "
              f"{sorted(components)}")
        return BackfillResult(BackfillOutcome.BACKFILLED, store_dir, sorted(components))

    # ── Build the new store content in a staging dir, then swap it in ─────────
    # Every file (augustus tree, genemark.mod, snap.hmm, parameters.json,
    # provenance.json) is written to a private staging dir first. Only once
    # that's fully populated do we touch the live store_dir, so a concurrent
    # reader (a sibling strain's FUNANNOTATE_PREDICT, or a racing backfill from
    # another project directory sharing this store) never observes a half
    # written augustus config tree.
    shared_root.mkdir(parents=True, exist_ok=True)
    staging_dir = shared_root / f".{species_tag}.staging.{os.getpid()}"
    if staging_dir.exists():
        shutil.rmtree(staging_dir)
    staging_dir.mkdir(parents=True)

    try:
        params_json = {
            "augustus": [{}], "genemark": [{}], "snap": [{}],
            "codingquarry": [{}], "glimmerhmm": [{}], "table": 1,
        }

        if "augustus" in components:
            dest = staging_dir / augustus_name
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
                "path": str((store_dir / augustus_name).resolve()),
            }]

        if "genemark" in components:
            dest = staging_dir / f"{species_tag}.genemark.mod"
            shutil.copyfile(components["genemark"], dest)
            params_json["genemark"] = [{
                "source": "ab-initio-reuse", "representative": rep_out,
                "path": str((store_dir / f"{species_tag}.genemark.mod").resolve()),
            }]

        if "snap" in components:
            dest = staging_dir / f"{species_tag}.snap.hmm"
            shutil.copyfile(components["snap"], dest)
            params_json["snap"] = [{
                "source": "ab-initio-reuse", "representative": rep_out,
                "path": str((store_dir / f"{species_tag}.snap.hmm").resolve()),
            }]

        glimmerhmm_stub = staging_dir / "glimmerhmm_stub"
        glimmerhmm_stub.mkdir(exist_ok=True)
        params_json["glimmerhmm"] = [{
            "source": "suppressed-empty-stub",
            "path": str((store_dir / "glimmerhmm_stub").resolve()),
        }]

        # ── Idempotency: skip the swap entirely if content hasn't changed ─────
        # A store's mtime is a real signal downstream (staleSharedParams() in
        # funannotate/utils.nf treats parameters.json being newer than a
        # sibling's GBK as "please re-predict me"), so a no-op sweep run must
        # never touch it -- only a genuine change to the representative's
        # trained parameters should bump it and trigger that cascade.
        new_hash = _store_content_hash(staging_dir, augustus_name, species_tag)
        old_hash = _store_content_hash(store_dir, augustus_name, species_tag)
        if store_dir.is_dir() and new_hash == old_hash:
            print(f"[INFO] {store_dir} already up to date from representative "
                  f"{rep_out} (content hash unchanged) — leaving it untouched")
            return BackfillResult(BackfillOutcome.UP_TO_DATE, store_dir, sorted(components))

        (staging_dir / "parameters.json").write_text(json.dumps(params_json, indent=2))
        provenance = {
            "representative_out": rep_out,
            "species": species,
            "augustus_species_name": augustus_name,
            "components": sorted(components),
            "ani_reuse_threshold": threshold,
            "date_captured": date.today().isoformat(),
            "generated_at": datetime.now().isoformat(),
            "source_predict_dir": str(rep_dir),
            "content_hash": new_hash,
        }
        (staging_dir / "provenance.json").write_text(json.dumps(provenance, indent=2))

        # Atomic-per-step swap: os.rename() of a directory is atomic, but POSIX
        # rename() refuses to replace a non-empty directory in one call, so an
        # existing store has to move out of the way before the new one moves in.
        # That leaves a very small window where store_dir doesn't exist at all
        # (never a torn/half-written one, which is the failure mode this design
        # actually guards against) -- acceptable, see .living/decisions.md.
        backup_dir = None
        if store_dir.exists():
            backup_dir = shared_root / f".{species_tag}.old.{os.getpid()}"
            os.rename(store_dir, backup_dir)
        os.rename(staging_dir, store_dir)
        if backup_dir is not None:
            shutil.rmtree(backup_dir)
    finally:
        if staging_dir.exists():
            shutil.rmtree(staging_dir)

    if aug_cfg is not None and "augustus" in components:
        link = aug_cfg / "species" / augustus_name
        link.parent.mkdir(parents=True, exist_ok=True)
        if link.is_symlink() or link.exists():
            if link.is_symlink() or link.is_file():
                link.unlink()
            else:
                shutil.rmtree(link)
        link.symlink_to((store_dir / augustus_name).resolve())
        print(f"[INFO] Symlinked {link} -> {(store_dir / augustus_name).resolve()}")

    print(f"[INFO] Backfilled {store_dir} from representative {rep_out} "
          f"(components: {sorted(components)}, augustus_name={augustus_name})")
    return BackfillResult(BackfillOutcome.BACKFILLED, store_dir, sorted(components))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--representative-out", required=True)
    ap.add_argument("--target", required=True)
    ap.add_argument("--shared-root", required=True)
    ap.add_argument("--ani-threshold", type=float, default=99.0)
    ap.add_argument("--augustus-config", default="")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    result = backfill_species_store(
        species=args.species,
        rep_out=args.representative_out,
        target=Path(args.target),
        shared_root=Path(args.shared_root),
        threshold=args.ani_threshold,
        aug_cfg=Path(args.augustus_config) if args.augustus_config else None,
        dry_run=args.dry_run,
    )

    print(f"[RESULT] {result.outcome.value}")
    if result.outcome is BackfillOutcome.REPRESENTATIVE_NOT_PREDICTED:
        # Both callers only invoke this script once they believe the
        # representative has already been predicted (the standalone sweep
        # filters on an on-disk GBK; the wired path runs after that
        # representative's own FUNANNOTATE_PREDICT channel emission) -- so
        # reaching this is a real inconsistency, not a normal "not yet" skip.
        print(f"[ERROR] {args.species}: representative {args.representative_out} "
              "has no usable prediction on disk", file=sys.stderr)
        sys.exit(1)
    sys.exit(0)


if __name__ == "__main__":
    main()
