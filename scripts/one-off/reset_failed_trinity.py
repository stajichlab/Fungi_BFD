#!/usr/bin/env python3
"""
reset_failed_trinity.py — reset species whose Trinity assembly failed so a pipeline
re-run re-assembles Trinity WITHOUT re-downloading RNA-seq.

A "failed Trinity" is a species whose cached Trinity-GG FASTA in rnaseq_data/ is
EMPTY (0 bytes) even though its normalized reads in rnaseq_reads/ are non-empty.
Empty Trinity files with empty/absent reads are the legitimate "no RNA-seq" sentinel
(funannotate.nf RNASEQ_PREPARE writes touch <tag>.trinity-GG.fasta in that case) and
are left untouched.

Why this is safe re: re-download
--------------------------------
The pipeline uses two independent storeDir caches:
  * SRA_FETCH      -> storeDir rnaseq_reads/  (downloads + trims + bbnorm-normalizes)
                      outputs <tag>_norm_{R1,R2,SE}.fastq.gz
  * RNASEQ_PREPARE -> storeDir rnaseq_data/   (runs Trinity)
                      output  <tag>.trinity-GG.fasta
Nextflow skips a storeDir process only when ALL its declared outputs already exist.
Deleting an empty rnaseq_data/<tag>.trinity-GG.fasta re-triggers RNASEQ_PREPARE only
(a Trinity re-assembly fed by the already-cached normalized reads). SRA_FETCH stays
cached because its rnaseq_reads/ outputs are intact, so NO NCBI re-download happens.
(A re-download would only be triggered by deleting files in rnaseq_reads/.)

Stale-training short-circuit
----------------------------
RNASEQ_PREPARE has a short-circuit (funannotate.nf:1026): if the representative's
training dir already contains funannotate_train.pasa.gff3 (produced by a prior
FUNANNOTATE_TRAIN run), it re-extracts trinity.fasta from there INSTEAD of re-running
Trinity -- and if that trinity.fasta is missing/empty it just re-touches an empty
output. So for a genuine failure with such a stale dir present, deleting only the
rnaseq_data file would reproduce the empty file. This script detects that risk and,
with --clean-stale-training, also removes the offending training dir to force a full
Trinity re-run.

Usage
-----
    # dry-run (default): report only, change nothing
    python scripts/reset_failed_trinity.py

    # actually remove the empty rnaseq_data Trinity files
    python scripts/reset_failed_trinity.py --apply

    # also remove stale training dirs that would short-circuit back to empty
    python scripts/reset_failed_trinity.py --apply --clean-stale-training

Options
-------
    --apply                 Perform removals (default: dry-run, no changes)
    --clean-stale-training  With --apply, also rm -rf representative training dirs that
                            have pasa.gff3 but no usable trinity.fasta (default: off)
    --rnaseq-data DIR       Trinity FASTA cache       (default: rnaseq_data)
    --rnaseq-reads DIR      Normalized reads cache     (default: rnaseq_reads)
    --training-dir DIR      Training target dir        (default: genome_annotation_training)
    --samples CSV           samples.csv for representative mapping (default: samples.csv)
    --report TSV            Report path (default: misc/trinity_failed_reset_report.tsv)
    --project-dir DIR       Project root (default: current directory)
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from pathlib import Path


# ── Sample-tag helpers: Python ports of nextflow/lib/SampleUtils.groovy ──────────
# Keep in sync with SampleUtils.cleanStrain / makeSampleTag.

def clean_strain(raw_strain: str) -> str:
    s = (raw_strain or "").strip()
    s = re.sub(r"['\"]", "", s)
    s = s.split(";")[0].strip()
    s = s.replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)      # leading '*' -> removed
    s = re.sub(r"\*+\s*$", "", s)      # trailing '*' -> removed
    s = re.sub(r"\s*\*+\s*", "-", s)   # internal '*' -> '-'
    return s.strip()


def make_sample_tag(raw_species: str, raw_strain: str) -> str:
    """Filesystem-safe '{species}_{strain}' tag == the training-dir 'out' name."""
    sp = re.sub(r"['\"]", "", (raw_species or "").strip())
    st = clean_strain(raw_strain)
    tag = "_".join(p for p in (sp, st) if p)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def species_tag_from_species(raw_species: str) -> str:
    """species_tag == trinity filename stem (funannotate.nf maps species \\s+ -> '_')."""
    return re.sub(r"\s+", "_", (raw_species or "").strip())


def species_tag_from_filename(fname: str) -> str:
    return fname.replace(".trinity-GG.fasta", "")


def load_representatives(samples_path: Path):
    """Map species_tag -> representative 'out' (first samples.csv row for that species).

    Returns {} if the file is missing; the caller falls back to glob-based discovery.
    """
    reps: dict[str, str] = {}
    if not samples_path.exists():
        return reps
    with open(samples_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            stag = species_tag_from_species(row.get("SPECIES", ""))
            if stag and stag not in reps:  # first row wins == representative
                reps[stag] = make_sample_tag(row.get("SPECIES", ""), row.get("STRAIN", ""))
    return reps


def reads_nonempty(reads_dir: Path, tag: str) -> tuple[bool, dict]:
    """A species has usable reads if R1 OR SE is non-empty (mirrors funannotate.nf:2097)."""
    sizes = {}
    for suffix in ("R1", "R2", "SE"):
        p = reads_dir / f"{tag}_norm_{suffix}.fastq.gz"
        sizes[suffix] = p.stat().st_size if p.exists() else None
    has = (sizes["R1"] or 0) > 0 or (sizes["SE"] or 0) > 0
    return has, sizes


def candidate_training_dirs(training_dir: Path, tag: str, rep_out: str | None) -> list[Path]:
    """Representative training dir(s). Prefer the exact 'out' from samples.csv; else
    glob '<tag>' and '<tag>_*' as a fallback (covers all strains of the species)."""
    if rep_out:
        d = training_dir / rep_out
        return [d] if d.is_dir() else []
    found = []
    exact = training_dir / tag
    if exact.is_dir():
        found.append(exact)
    found.extend(sorted(d for d in training_dir.glob(f"{tag}_*") if d.is_dir()))
    return found


def training_status(tdir: Path) -> str:
    """Classify a training dir's effect on RNASEQ_PREPARE re-run.

    full_rerun   : no pasa.gff3 -> RNASEQ_PREPARE runs Trinity from scratch (good)
    reextract_ok : pasa.gff3 present AND non-empty trinity.fasta -> re-extract fixes it
    stale_empty  : pasa.gff3 present but trinity.fasta missing/empty -> short-circuits
                   back to an empty file (must remove dir to force a real re-run)
    """
    pasa = tdir / "training" / "funannotate_train.pasa.gff3"
    if not pasa.exists():
        return "full_rerun"
    trinity = tdir / "training" / "trinity.fasta"
    if trinity.exists() and trinity.stat().st_size > 0:
        return "reextract_ok"
    return "stale_empty"


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--apply", action="store_true",
                    help="Perform removals (default: dry-run)")
    ap.add_argument("--clean-stale-training", action="store_true",
                    help="With --apply, also rm -rf stale representative training dirs")
    ap.add_argument("--rnaseq-data", default="rnaseq_data")
    ap.add_argument("--rnaseq-reads", default="rnaseq_reads")
    ap.add_argument("--training-dir", default="genome_annotation_training")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--report", default="misc/trinity_failed_reset_report.tsv")
    ap.add_argument("--project-dir", default=".")
    args = ap.parse_args()

    root = Path(args.project_dir).resolve()
    rnaseq_data = root / args.rnaseq_data
    rnaseq_reads = root / args.rnaseq_reads
    training_dir = root / args.training_dir
    samples_path = root / args.samples
    report_path = root / args.report

    if not rnaseq_data.is_dir():
        sys.exit(f"[ERROR] rnaseq_data not found: {rnaseq_data}")
    if not rnaseq_reads.is_dir():
        sys.exit(f"[ERROR] rnaseq_reads not found: {rnaseq_reads}")

    reps = load_representatives(samples_path)
    if reps:
        print(f"[INFO] Loaded representative mapping for {len(reps)} species from {samples_path}")
    else:
        print(f"[WARN] No samples.csv at {samples_path}; training dirs discovered by glob.")

    empty_trinity = sorted(
        p for p in rnaseq_data.glob("*.trinity-GG.fasta") if p.stat().st_size == 0
    )
    print(f"[INFO] {len(empty_trinity)} empty (0-byte) Trinity files in {rnaseq_data.name}/")

    rows = []
    for fa in empty_trinity:
        tag = species_tag_from_filename(fa.name)
        has_reads, sizes = reads_nonempty(rnaseq_reads, tag)
        if not has_reads:
            continue  # legitimate no-RNA-seq sentinel; skip
        rep_out = reps.get(tag)
        tdirs = candidate_training_dirs(training_dir, tag, rep_out)
        statuses = {str(d.relative_to(root)): training_status(d) for d in tdirs}
        stale = [d for d in tdirs if training_status(d) == "stale_empty"]
        rows.append({
            "species_tag": tag,
            "R1_bytes": sizes["R1"], "R2_bytes": sizes["R2"], "SE_bytes": sizes["SE"],
            "rep_out": rep_out or "(glob)",
            "training_status": ";".join(f"{k}={v}" for k, v in statuses.items()) or "no_dir",
            "stale_dirs": stale,
            "trinity_path": fa,
        })

    print(f"[INFO] {len(rows)} failed-Trinity species (empty Trinity, non-empty reads)\n")

    # ── Report ───────────────────────────────────────────────────────────────────
    report_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species_tag", "R1_bytes", "R2_bytes", "SE_bytes",
              "rep_out", "training_status"]
    with open(report_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow({k: r[k] for k in fields})
    print(f"[INFO] Report written to {report_path}\n")

    # ── Plan / apply ─────────────────────────────────────────────────────────────
    mode = "APPLY" if args.apply else "DRY-RUN"
    n_stale = sum(1 for r in rows if r["stale_dirs"])
    for r in rows:
        print(f"  {r['species_tag']}")
        print(f"      rm {r['trinity_path'].relative_to(root)}")
        for d in r["stale_dirs"]:
            flag = "" if args.clean_stale_training else "  (needs --clean-stale-training)"
            print(f"      rm -rf {d.relative_to(root)}  [stale short-circuit]{flag}")

    print(f"\n[{mode}] {len(rows)} empty Trinity files; "
          f"{n_stale} species have a stale training short-circuit.")

    if not args.apply:
        print("[DRY-RUN] No changes made. Re-run with --apply to remove the empty "
              "Trinity files\n          (add --clean-stale-training to also clear stale "
              "training dirs).")
        return

    removed_fa = 0
    removed_dirs = 0
    for r in rows:
        r["trinity_path"].unlink()
        removed_fa += 1
        if args.clean_stale_training:
            for d in r["stale_dirs"]:
                shutil.rmtree(d)
                removed_dirs += 1
    print(f"\n[APPLY] Removed {removed_fa} empty Trinity files"
          + (f" and {removed_dirs} stale training dirs." if args.clean_stale_training
             else f"; {n_stale} stale training dirs left in place "
                  f"(re-run with --clean-stale-training to remove)."))
    print("[APPLY] Re-run the funannotate pipeline to re-assemble Trinity "
          "(reads are cached; no re-download).")


if __name__ == "__main__":
    main()
