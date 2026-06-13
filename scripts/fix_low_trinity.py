#!/usr/bin/env python3
"""
fix_low_trinity.py — identify low-count trinity assemblies and propose better strains.

Usage:
    python scripts/fix_low_trinity.py [options]

Options:
    --threshold INT     Transcript count threshold (default: 2000)
    --rnaseq-dir DIR    Directory with trinity-GG.fasta files (default: rnaseq_data)
    --samples CSV       samples.csv path (default: samples.csv)
    --asm-stats TSV     Assembly stats file (default: tables/asm_stats.tsv.gz)
    --busco-csv CSV     BUSCO scores file (default: tables/BUSCO.csv.gz, optional)
    --report TSV        Output report path (default: misc/trinity_lowcount_report.tsv)
    --rerun-csv CSV     Output rerun samples CSV (default: trinity_rerun.samples.csv)
    --move              Actually move low-count files to misc/poor_trinity/ (default: dry-run)
    --poor-dir DIR      Destination for moved files (default: misc/poor_trinity)
    --project-dir DIR   Project root (default: current directory)
"""

import argparse
import csv
import gzip
import os
import re
import shutil
import sys
from pathlib import Path


def count_transcripts(fasta_path: Path) -> int:
    """Count sequences in a FASTA file by counting '>' lines."""
    count = 0
    opener = gzip.open if str(fasta_path).endswith(".gz") else open
    try:
        with opener(fasta_path, "rt", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    count += 1
    except (OSError, EOFError):
        return -1
    return count


def species_tag_from_filename(fname: str) -> str:
    """Extract species_tag from trinity FASTA filename.
    E.g. 'Saccharomyces_uvarum.trinity-GG.fasta' -> 'Saccharomyces_uvarum'
    """
    return fname.replace(".trinity-GG.fasta", "")


def species_key_from_row(row: dict) -> str:
    """Derive the species_tag key from a samples.csv row (matches trinity filename stem)."""
    return row["SPECIES"].replace(" ", "_")


def find_best_assembly(strains: list[dict], asm_stats: dict, busco_scores: dict) -> tuple[dict, str]:
    """
    Choose the best strain from a list of samples.csv rows.

    Priority:
      1. GCF_ prefix (NCBI reference assembly)
      2. Highest BUSCO complete_pct (if busco_scores available)
      3. Highest N50_bp from asm_stats

    Returns (best_row, reason_string).
    """
    if len(strains) == 1:
        return strains[0], "only_strain"

    # Priority 1: GCF_ reference assemblies
    gcf_strains = [s for s in strains if s["ASMID"].startswith("GCF_")]
    if gcf_strains:
        # Among GCF_ assemblies, still pick best by BUSCO then N50
        best = _rank_by_busco_n50(gcf_strains, asm_stats, busco_scores)
        return best, "GCF_reference"

    # Priority 2 & 3: BUSCO then N50
    best = _rank_by_busco_n50(strains, asm_stats, busco_scores)
    reason = "BUSCO" if busco_scores else "N50"
    return best, reason


def _rank_by_busco_n50(strains: list[dict], asm_stats: dict, busco_scores: dict) -> dict:
    """Sort strains by BUSCO complete_pct (desc) then N50_bp (desc), return best."""

    def sort_key(s):
        asmid = s["ASMID"]
        busco = busco_scores.get(asmid, {}).get("complete_pct", -1.0) if busco_scores else -1.0
        n50 = asm_stats.get(asmid, {}).get("N50_bp", 0)
        return (busco, n50)

    return max(strains, key=sort_key)


def load_asm_stats(path: Path) -> dict:
    """Load asm_stats.tsv.gz into dict keyed by ASMID."""
    stats = {}
    if not path.exists():
        return stats
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
    """Load BUSCO.csv.gz into dict keyed by ASMID. Returns empty dict if file absent."""
    scores = {}
    if not path.exists():
        return scores
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


def load_samples(path: Path) -> tuple[dict, list]:
    """
    Load samples.csv.

    Returns:
        species_map: dict of species_tag -> list[row] (in file order, first = representative)
        all_rows: all rows in order (for writing rerun CSV with correct columns)
    """
    species_map: dict[str, list] = {}
    all_rows: list[dict] = []
    fieldnames: list[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        fieldnames = reader.fieldnames or []
        for row in reader:
            all_rows.append(row)
            key = species_key_from_row(row)
            species_map.setdefault(key, []).append(row)
    return species_map, all_rows, fieldnames


def scan_trinity_files(rnaseq_dir: Path, threshold: int) -> list[tuple[Path, int]]:
    """Return list of (fasta_path, count) for files with count < threshold."""
    low = []
    fasta_files = sorted(rnaseq_dir.glob("*.trinity-GG.fasta"))
    total = len(fasta_files)
    print(f"[INFO] Scanning {total} trinity FASTA files in {rnaseq_dir} ...", flush=True)
    for i, fa in enumerate(fasta_files, 1):
        if i % 500 == 0:
            print(f"  {i}/{total} ...", flush=True)
        n = count_transcripts(fa)
        if n < threshold:
            low.append((fa, n))
    print(f"[INFO] Found {len(low)} files with fewer than {threshold} transcripts.")
    return low


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--threshold", type=int, default=2000)
    parser.add_argument("--rnaseq-dir", default="rnaseq_data")
    parser.add_argument("--samples", default="samples.csv")
    parser.add_argument("--asm-stats", default="tables/asm_stats.tsv.gz")
    parser.add_argument("--busco-csv", default="tables/BUSCO.csv.gz")
    parser.add_argument("--report", default="misc/trinity_lowcount_report.tsv")
    parser.add_argument("--rerun-csv", default="trinity_rerun.samples.csv")
    parser.add_argument("--move", action="store_true", help="Move low-count files to --poor-dir (default: dry-run)")
    parser.add_argument("--poor-dir", default="misc/poor_trinity")
    parser.add_argument("--project-dir", default=".")
    args = parser.parse_args()

    root = Path(args.project_dir).resolve()
    rnaseq_dir = root / args.rnaseq_dir
    samples_path = root / args.samples
    asm_stats_path = root / args.asm_stats
    busco_path = root / args.busco_csv
    report_path = root / args.report
    rerun_path = root / args.rerun_csv
    poor_dir = root / args.poor_dir

    if not rnaseq_dir.is_dir():
        print(f"[ERROR] rnaseq_dir not found: {rnaseq_dir}", file=sys.stderr)
        sys.exit(1)
    if not samples_path.exists():
        print(f"[ERROR] samples.csv not found: {samples_path}", file=sys.stderr)
        sys.exit(1)

    # Load data
    print("[INFO] Loading samples.csv ...")
    species_map, all_rows, fieldnames = load_samples(samples_path)

    print(f"[INFO] Loading assembly stats from {asm_stats_path} ...")
    asm_stats = load_asm_stats(asm_stats_path)
    print(f"  Loaded stats for {len(asm_stats)} assemblies.")

    if busco_path.exists():
        print(f"[INFO] Loading BUSCO scores from {busco_path} ...")
        busco_scores = load_busco_scores(busco_path)
        print(f"  Loaded BUSCO scores for {len(busco_scores)} assemblies.")
    else:
        print(f"[INFO] BUSCO file not found ({busco_path}); skipping BUSCO scoring.")
        busco_scores = {}

    # Scan trinity files
    low_count_files = scan_trinity_files(rnaseq_dir, args.threshold)

    if not low_count_files:
        print("[INFO] No low-count trinity files found. Nothing to do.")
        return

    # Build report
    report_rows = []
    rerun_rows = []

    for fasta_path, tx_count in sorted(low_count_files, key=lambda x: x[1]):
        species_tag = species_tag_from_filename(fasta_path.name)
        strains = species_map.get(species_tag, [])

        if not strains:
            report_rows.append({
                "species_tag": species_tag,
                "trinity_file": fasta_path.name,
                "transcript_count": tx_count,
                "current_asmid": "NOT_IN_SAMPLES_CSV",
                "current_strain": "",
                "strain_count": 0,
                "best_asmid": "",
                "best_strain": "",
                "best_n50": "",
                "best_busco_complete": "",
                "reason": "no_match_in_samples_csv",
            })
            continue

        current = strains[0]  # representative strain (first in file order)
        best, reason = find_best_assembly(strains, asm_stats, busco_scores)

        best_asmid = best["ASMID"]
        best_stats = asm_stats.get(best_asmid, {})
        best_busco = busco_scores.get(best_asmid, {})

        report_rows.append({
            "species_tag": species_tag,
            "trinity_file": fasta_path.name,
            "transcript_count": tx_count,
            "current_asmid": current["ASMID"],
            "current_strain": current["STRAIN"],
            "strain_count": len(strains),
            "best_asmid": best_asmid,
            "best_strain": best["STRAIN"],
            "best_n50": best_stats.get("N50_bp", ""),
            "best_busco_complete": best_busco.get("complete_pct", ""),
            "reason": reason,
        })

        # Add to rerun list only if we have a different/better strain to switch to
        if reason != "only_strain":
            rerun_rows.append(best)

    # Write report
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_fields = [
        "species_tag", "trinity_file", "transcript_count",
        "current_asmid", "current_strain", "strain_count",
        "best_asmid", "best_strain", "best_n50", "best_busco_complete", "reason",
    ]
    with open(report_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=report_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(report_rows)
    print(f"[INFO] Report written to {report_path} ({len(report_rows)} entries)")

    # Write trinity_rerun.samples.csv
    if rerun_rows:
        # Deduplicate by ASMID (a species might appear once anyway, but be safe)
        seen = set()
        deduped = []
        for row in rerun_rows:
            if row["ASMID"] not in seen:
                seen.add(row["ASMID"])
                deduped.append(row)
        with open(rerun_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(deduped)
        print(f"[INFO] Rerun samples CSV written to {rerun_path} ({len(deduped)} rows)")
    else:
        print("[INFO] No strains to rerun (all low-count species have only one strain).")

    # Move / dry-run
    if args.move:
        poor_dir.mkdir(parents=True, exist_ok=True)
        moved = 0
        for fasta_path, tx_count in low_count_files:
            dest = poor_dir / fasta_path.name
            shutil.move(str(fasta_path), str(dest))
            moved += 1
        print(f"[INFO] Moved {moved} files to {poor_dir}")
    else:
        print(f"\n[DRY-RUN] {len(low_count_files)} files would be moved to {poor_dir}")
        print("  Re-run with --move to perform the actual move.")
        for fasta_path, tx_count in sorted(low_count_files, key=lambda x: x[1]):
            print(f"  {tx_count:>6}  {fasta_path.name}")

    # Summary
    only_strain = sum(1 for r in report_rows if r["reason"] == "only_strain")
    has_better = sum(1 for r in report_rows if r["reason"] != "only_strain" and r["reason"] != "no_match_in_samples_csv")
    no_match = sum(1 for r in report_rows if r["reason"] == "no_match_in_samples_csv")
    gcf_wins = sum(1 for r in report_rows if r["reason"] == "GCF_reference")
    print(f"""
Summary:
  Total low-count files (<{args.threshold}):  {len(report_rows)}
  Has better alternative strain:       {has_better}
    - chosen because GCF_ reference:  {gcf_wins}
    - chosen by BUSCO/N50:            {has_better - gcf_wins}
  Only one strain (no alternative):    {only_strain}
  Not found in samples.csv:            {no_match}
""")


if __name__ == "__main__":
    main()
