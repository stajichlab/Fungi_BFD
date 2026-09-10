#!/usr/bin/env python3
"""
pick_rnaseq_representative_override.py — detect species whose shared Trinity-GG
assembly (rnaseq_data/<species_tag>.trinity-GG.fasta) is too thin to train on, and
pick a different strain for RNASEQ_PREPARE to rebuild it against.

Why this exists (separate from fix_low_trinity.py):
FUNANNOTATE_RNASEQ.nf normally builds the one shared Trinity-GG assembly for a
species against whichever strain abinitio_reuse_assignments.csv marks
is_representative=True -- an ANI+BUSCO pick optimized for assembly quality, not for
whether the species' actual RNA-seq reads align to that particular genome. The two
occasionally disagree badly: Ascochyta_rabiei's ANI pick (a pks1-deletion construct
genome, BUSCO 98.1%) is a perfectly fine assembly that the species' real RNA-seq
(SRR330019xx) barely aligns to -- HISAT2 -> 9 Trinity clusters -> 7 transcripts total
-- while GCF_004011695.2 (Me14, BUSCO 98.4%, the RefSeq reference for this species)
works. fix_low_trinity.py's own "current" pick (samples.csv row order) usually, but
does not provably, match the actual ANI-picked representative -- this script reads
abinitio_reuse_assignments.csv directly so the "current" it excludes is the strain
that was ACTUALLY tried, and reports the alternate as an `out` value that
FUNANNOTATE_RNASEQ.nf's rnaseqRepOverride can act on directly (see
loadRnaseqRepresentativeOverride() in nextflow/modules/funannotate/utils.nf).

What it does:
  1. Scan rnaseq_data/*.trinity-GG.fasta for transcript counts below --threshold.
  2. For each low-count species, find its current representative `out` from
     abinitio_reuse_assignments.csv (is_representative=True), match it to a
     samples.csv ASMID by strain-name, and pick the best ALTERNATE strain
     (GCF_ reference preferred, then BUSCO complete_pct, then N50) excluding
     that current representative.
  3. Match the alternate's ASMID back to its `out` value via
     abinitio_reuse_assignments.csv (species + normalized-strain match).
  4. Write rnaseq_representative_override.csv (species_tag,out) for Nextflow to
     consume, plus a report TSV with the reasoning, and (unless --dry-run) move the
     stale low-count trinity-GG.fasta out of rnaseq_data/ so RNASEQ_PREPARE's
     storeDir cache doesn't skip rebuilding it.

Species with only one strain, or where the alternate can't be matched to an `out`,
are reported but left out of the override CSV -- funannotate_train's
train_min_trinity_transcripts gate (see FUNANNOTATE_TRAIN/main.nf) is the backstop
for those; there is nothing better to switch to.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
import sys
from pathlib import Path


def normalize_strain(s: str) -> str:
    """Loose approximation of nextflow/modules/common/utils.nf's cleanStrain+makeSampleTag
    normalization, good enough to match a samples.csv STRAIN to an `out` suffix:
    lowercase, first semicolon token, strip everything but alphanumerics."""
    s = (s or "").split(";")[0]
    return re.sub(r"[^a-z0-9]", "", s.lower())


def count_transcripts(fasta_path: Path) -> int:
    count = 0
    with open(fasta_path, errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                count += 1
    return count


def load_samples(path: Path):
    """species_tag -> list of samples.csv rows (species-name-based tag, spaces -> _)."""
    species_map: dict[str, list] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            tag = row["SPECIES"].replace(" ", "_")
            species_map.setdefault(tag, []).append(row)
    return species_map


def load_abinitio_reuse(path: Path):
    """species_tag -> [{out, is_representative}, ...] from abinitio_reuse_assignments.csv."""
    by_tag: dict[str, list] = {}
    if not path.exists():
        return by_tag
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            tag = row["species"].strip().replace(" ", "_")
            by_tag.setdefault(tag, []).append({
                "out": row["out"].strip(),
                "is_representative": row.get("is_representative", "").strip().lower() == "true",
            })
    return by_tag


def load_asm_stats(path: Path) -> dict:
    stats = {}
    if not path.exists():
        return stats
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            asmid = row.get("ASMID", "").strip()
            if asmid:
                try:
                    row["N50_bp"] = int(row.get("N50_bp", 0) or 0)
                except ValueError:
                    row["N50_bp"] = 0
                stats[asmid] = row
    return stats


def load_busco_scores(path: Path) -> dict:
    scores = {}
    if not path.exists():
        return scores
    import gzip
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = row.get("ASMID", "").strip()
            if asmid:
                try:
                    row["complete_pct"] = float(row.get("complete_pct", 0) or 0)
                except ValueError:
                    row["complete_pct"] = 0.0
                scores[asmid] = row
    return scores


def rank_by_gcf_busco_n50(strains: list[dict], asm_stats: dict, busco_scores: dict) -> dict:
    gcf = [s for s in strains if s["ASMID"].startswith("GCF_")]
    pool = gcf if gcf else strains

    def sort_key(s):
        asmid = s["ASMID"]
        busco = busco_scores.get(asmid, {}).get("complete_pct", -1.0) if busco_scores else -1.0
        n50 = asm_stats.get(asmid, {}).get("N50_bp", 0)
        return (busco, n50)

    return max(pool, key=sort_key), ("GCF_reference" if gcf else ("BUSCO" if busco_scores else "N50"))


def find_out_for_asmid(asmid_row: dict, tag: str, reuse_outs: list[dict]) -> str | None:
    """Match a samples.csv row's STRAIN to one of abinitio_reuse's `out` values for
    this species_tag by normalized-strain suffix (out == tag for a blank strain)."""
    target = normalize_strain(asmid_row.get("STRAIN", ""))
    for entry in reuse_outs:
        out = entry["out"]
        suffix = out[len(tag) + 1:] if out.startswith(tag + "_") else ("" if out == tag else None)
        if suffix is None:
            continue
        if normalize_strain(suffix) == target:
            return out
    return None


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--threshold", type=int, default=2000)
    p.add_argument("--rnaseq-data", default="rnaseq_data")
    p.add_argument("--samples", default="samples.csv")
    p.add_argument("--abinitio-reuse-csv", default="genome_annotation/_reuse_assignments/abinitio_reuse_assignments.csv")
    p.add_argument("--asm-stats", default="tables/asm_stats.tsv.gz")
    p.add_argument("--busco-csv", default="tables/BUSCO.csv.gz")
    p.add_argument("--override-csv", default="rnaseq_representative_override.csv")
    p.add_argument("--report", default="misc/rnaseq_representative_override_report.tsv")
    p.add_argument("--poor-dir", default="misc/poor_trinity")
    p.add_argument("--project-dir", default=".")
    p.add_argument("--dry-run", action="store_true",
                    help="Report only -- do not write the override CSV or move stale files")
    args = p.parse_args()

    root = Path(args.project_dir).resolve()
    rnaseq_data = root / args.rnaseq_data
    samples_path = root / args.samples
    abinitio_path = root / args.abinitio_reuse_csv
    override_path = root / args.override_csv
    report_path = root / args.report
    poor_dir = root / args.poor_dir

    species_map = load_samples(samples_path)
    reuse_by_tag = load_abinitio_reuse(abinitio_path)
    asm_stats = load_asm_stats(root / args.asm_stats)
    busco_scores = load_busco_scores(root / args.busco_csv)

    if not reuse_by_tag:
        print(f"[WARN] No ab-initio reuse assignments found at {abinitio_path}; "
              f"cannot identify current representatives, aborting.", file=sys.stderr)
        sys.exit(1)

    fasta_files = sorted(rnaseq_data.glob("*.trinity-GG.fasta"))
    print(f"[INFO] Scanning {len(fasta_files)} trinity-GG fasta files in {rnaseq_data} ...", file=sys.stderr)

    override_rows = []
    report_rows = []
    to_move = []

    for fasta in fasta_files:
        tag = fasta.name[: -len(".trinity-GG.fasta")]
        # A 0-byte file is RNASEQ_PREPARE's "no RNA-seq reads found for this species at
        # all" marker (see its own no-reads branch), not a too-few-transcripts assembly --
        # switching representative strain can't fix "there are no reads", so these must be
        # excluded here rather than folded into n_tx=0 and treated like a real but thin
        # assembly.
        if fasta.stat().st_size == 0:
            continue
        n_tx = count_transcripts(fasta)
        if n_tx >= args.threshold:
            continue

        strains = species_map.get(tag, [])
        reuse_outs = reuse_by_tag.get(tag, [])
        row = {
            "species_tag": tag, "transcript_count": n_tx, "strain_count": len(strains),
            "current_out": "", "alternate_asmid": "", "alternate_strain": "",
            "alternate_out": "", "reason": "",
        }

        if len(strains) < 2:
            row["reason"] = "only_one_strain"
            report_rows.append(row)
            continue

        current_entry = next((e for e in reuse_outs if e["is_representative"]), None)
        if current_entry is None:
            row["reason"] = "no_representative_in_reuse_csv"
            report_rows.append(row)
            continue
        row["current_out"] = current_entry["out"]

        current_asmid = None
        for s in strains:
            if find_out_for_asmid(s, tag, [current_entry]) == current_entry["out"]:
                current_asmid = s["ASMID"]
                break

        candidates = [s for s in strains if s["ASMID"] != current_asmid]
        if not candidates:
            row["reason"] = "no_alternate_strain"
            report_rows.append(row)
            continue

        best, reason = rank_by_gcf_busco_n50(candidates, asm_stats, busco_scores)
        best_out = find_out_for_asmid(best, tag, reuse_outs)
        row["alternate_asmid"] = best["ASMID"]
        row["alternate_strain"] = best["STRAIN"]
        row["alternate_out"] = best_out or ""
        row["reason"] = reason if best_out else f"{reason}_but_unmatched_out"
        report_rows.append(row)

        if best_out:
            override_rows.append({"species_tag": tag, "out": best_out})
            to_move.append(fasta)

    report_path.parent.mkdir(parents=True, exist_ok=True)
    fields = ["species_tag", "transcript_count", "strain_count", "current_out",
              "alternate_asmid", "alternate_strain", "alternate_out", "reason"]
    with open(report_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(report_rows)
    print(f"[INFO] Report written to {report_path} ({len(report_rows)} low-count species)", file=sys.stderr)

    if args.dry_run:
        print(f"[DRY-RUN] Would write {len(override_rows)} overrides to {override_path} "
              f"and move {len(to_move)} stale trinity-GG fastas to {poor_dir}", file=sys.stderr)
        for r in override_rows:
            print(f"  {r['species_tag']} -> {r['out']}", file=sys.stderr)
        return

    if override_rows:
        with open(override_path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["species_tag", "out"])
            w.writeheader()
            w.writerows(override_rows)
        print(f"[INFO] Wrote {len(override_rows)} overrides to {override_path}", file=sys.stderr)

        poor_dir.mkdir(parents=True, exist_ok=True)
        counts_dir = rnaseq_data / "counts"
        n_counts_cleared = 0
        for fasta in to_move:
            shutil.move(str(fasta), str(poor_dir / fasta.name))
            tag = fasta.name[: -len(".trinity-GG.fasta")]
            count_file = counts_dir / f"{tag}.n_transcripts.txt"
            if count_file.exists():
                count_file.unlink()
                n_counts_cleared += 1
        print(f"[INFO] Moved {len(to_move)} stale trinity-GG fastas to {poor_dir} "
              f"(clears RNASEQ_PREPARE's storeDir cache so it rebuilds against the override)",
              file=sys.stderr)
        print(f"[INFO] Cleared {n_counts_cleared} stale cached counts from {counts_dir} "
              f"(clears COUNT_TRINITY_TRANSCRIPTS's storeDir cache so it recounts the rebuild)",
              file=sys.stderr)
    else:
        print("[INFO] No species had a usable alternate strain; nothing written.", file=sys.stderr)


if __name__ == "__main__":
    main()
