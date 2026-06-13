#!/usr/bin/env python3
"""
collect_busco_stats.py — collect BUSCO short summary files into tables/BUSCO.csv.gz.

Scans one or more BUSCO output directories and writes tables/BUSCO.csv.gz.

Usage:
    python scripts/collect_busco_stats.py [options]

Options:
    --busco-dir DIR        Legacy BUSCO results root; walks for short_summary*.txt
                           (default: busco_results; skipped if directory absent)
    --genome-busco-dir DIR BFD storeDir for BUSCO_GENOME results (default:
                           results/genome_stats/BUSCO_genome).  Files are named
                           {basename}.BUSCO_summary.{lineage}.txt
    --pep-busco-dir DIR    BFD storeDir for BUSCO_PEP results (default:
                           results/genome_stats/BUSCO_protein).  Same naming.
    --samples CSV          samples.csv path — used to map basename → ASMID
                           (default: samples.csv)
    --output FILE          Output CSV.gz path (default: tables/BUSCO.csv.gz)
    --project-dir DIR      Project root (default: current directory)

Source priority (later sources win):
  1. batch_summary.txt in --busco-dir
  2. per-run short_summary*.txt in --busco-dir subtree
  3. BFD-style {basename}.BUSCO_summary.{lineage}.txt from --genome-busco-dir
  4. BFD-style {basename}.BUSCO_summary.{lineage}.txt from --pep-busco-dir

Output columns:
    ASMID, complete_pct, single_pct, duplicated_pct, fragmented_pct, missing_pct,
    n_markers, lineage
"""

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path


# BUSCO v4/v5 short_summary patterns
_COMPLETE_RE = re.compile(r"C:(\d+\.\d+)%\[S:(\d+\.\d+)%,D:(\d+\.\d+)%\],F:(\d+\.\d+)%,M:(\d+\.\d+)%,n:(\d+)")
_LINEAGE_RE = re.compile(r"The lineage dataset is:\s*(\S+)")

# BFD storeDir filename: {basename}.BUSCO_summary.{lineage}.txt
_BFD_SUMMARY_RE = re.compile(r"^(.+)\.BUSCO_summary\.(.+)\.txt$")


def parse_short_summary(path: Path) -> dict | None:
    """Parse a BUSCO short_summary*.txt file, return dict of metrics or None."""
    text = path.read_text(errors="replace")
    m = _COMPLETE_RE.search(text)
    if not m:
        return None
    lineage = ""
    ml = _LINEAGE_RE.search(text)
    if ml:
        lineage = ml.group(1)
    else:
        parts = path.stem.split(".")
        if len(parts) >= 3:
            lineage = parts[2]
    return {
        "complete_pct": float(m.group(1)),
        "single_pct":   float(m.group(2)),
        "duplicated_pct": float(m.group(3)),
        "fragmented_pct": float(m.group(4)),
        "missing_pct":  float(m.group(5)),
        "n_markers":    int(m.group(6)),
        "lineage":      lineage,
    }


def parse_batch_summary(path: Path) -> dict[str, dict]:
    """Parse a BUSCO 5 batch_summary.txt; returns dict keyed by input file stem."""
    results = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            name = Path(row.get("Input_file", "")).stem
            name = re.sub(r"\.(fa|fna|fasta|masked)$", "", name, flags=re.IGNORECASE)
            try:
                results[name] = {
                    "complete_pct":   float(row.get("C", 0) or 0),
                    "single_pct":     float(row.get("S", 0) or 0),
                    "duplicated_pct": float(row.get("D", 0) or 0),
                    "fragmented_pct": float(row.get("F", 0) or 0),
                    "missing_pct":    float(row.get("M", 0) or 0),
                    "n_markers":      int(float(row.get("n", 0) or 0)),
                    "lineage":        row.get("Lineage_dataset", ""),
                }
            except (ValueError, KeyError):
                continue
    return results


def build_basename_map(samples_path: Path) -> dict[str, str]:
    """
    Build basename -> ASMID map from samples.csv.

    The basename is constructed the same way as the Nextflow workflow:
        species_strain (spaces/slashes/# replaced with _; strain uses first
        semicolon-delimited token with quotes and colons cleaned).
    """
    bmap: dict[str, str] = {}
    if not samples_path.exists():
        return bmap
    with open(samples_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid   = row.get("ASMID", "").strip()
            species = row.get("SPECIES", "").strip()
            strain  = row.get("STRAIN", "").strip().split(";")[0].strip()
            strain  = strain.replace("'", "").replace(":", " ")
            parts   = [p for p in [species, strain] if p]
            basename = "_".join(parts)
            basename = re.sub(r"[\s/\#]+", "_", basename)
            if asmid and basename:
                bmap[basename] = asmid
    return bmap


def load_asmid_map(samples_path: Path) -> dict[str, str]:
    """Map ASMID (and version-only prefix) -> canonical ASMID."""
    amap: dict[str, str] = {}
    if not samples_path.exists():
        return amap
    with open(samples_path) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            asmid = row.get("ASMID", "").strip()
            if not asmid:
                continue
            amap[asmid] = asmid
            ver_only = re.match(r"(GC[AF]_\d+\.\d+)", asmid)
            if ver_only:
                amap.setdefault(ver_only.group(1), asmid)
    return amap


def scan_bfd_storeDir(directory: Path, basename_map: dict[str, str],
                      records: dict[str, dict], label: str) -> int:
    """
    Parse {basename}.BUSCO_summary.{lineage}.txt files from a BFD storeDir.

    Updates `records` in place (basename_map used to resolve ASMID).
    Returns number of new entries added.
    """
    if not directory.is_dir():
        return 0
    added = 0
    for path in sorted(directory.glob("*.BUSCO_summary.*.txt")):
        m = _BFD_SUMMARY_RE.match(path.name)
        if not m:
            continue
        basename = m.group(1)
        stats = parse_short_summary(path)
        if not stats:
            continue
        # Prefer lineage encoded in filename over what's in the file text
        stats["lineage"] = m.group(2)
        asmid = basename_map.get(basename, basename)
        records[asmid] = stats
        added += 1
    if added:
        print(f"  [{label}] {added} entries from {directory}")
    return added


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--busco-dir",        default="busco_results")
    parser.add_argument("--genome-busco-dir", default="results/genome_stats/BUSCO_genome")
    parser.add_argument("--pep-busco-dir",    default="results/genome_stats/BUSCO_protein")
    parser.add_argument("--samples",  default="samples.csv")
    parser.add_argument("--output",   default="tables/BUSCO.csv.gz")
    parser.add_argument("--project-dir", default=".")
    args = parser.parse_args()

    root         = Path(args.project_dir).resolve()
    busco_dir    = root / args.busco_dir
    genome_dir   = root / args.genome_busco_dir
    pep_dir      = root / args.pep_busco_dir
    samples_path = root / args.samples
    output_path  = root / args.output

    asmid_map   = load_asmid_map(samples_path)
    basename_map = build_basename_map(samples_path)

    records: dict[str, dict] = {}

    # ── Source 1 & 2: legacy busco_results/ directory ─────────────────────────
    if busco_dir.is_dir():
        batch_file = busco_dir / "batch_summary.txt"
        if batch_file.exists():
            print(f"[INFO] Parsing batch summary: {batch_file}")
            batch = parse_batch_summary(batch_file)
            for name, stats in batch.items():
                asmid = asmid_map.get(name, name)
                records[asmid] = stats
            print(f"  {len(batch)} entries.")

        added = 0
        for summary_file in sorted(busco_dir.rglob("short_summary*.txt")):
            if summary_file == batch_file:
                continue
            stats = parse_short_summary(summary_file)
            if not stats:
                continue
            parent = summary_file.parent.name
            if parent.startswith("run_"):
                parent = summary_file.parent.parent.name
            asmid = asmid_map.get(parent, parent)
            if asmid not in records:
                records[asmid] = stats
                added += 1
        if added:
            print(f"[INFO] {added} entries from short_summary files in {busco_dir}")
    else:
        print(f"[INFO] --busco-dir not found ({busco_dir}), skipping.")

    # ── Sources 3 & 4: BFD storeDir flat directories ──────────────────────────
    scan_bfd_storeDir(genome_dir, basename_map, records, "BUSCO_genome")
    scan_bfd_storeDir(pep_dir,    basename_map, records, "BUSCO_protein")

    if not records:
        print("[WARN] No BUSCO data found — output will be header-only.", file=sys.stderr)

    print(f"[INFO] Total assemblies with BUSCO data: {len(records)}")

    fieldnames = ["ASMID", "complete_pct", "single_pct", "duplicated_pct",
                  "fragmented_pct", "missing_pct", "n_markers", "lineage"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output_path, "wt", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for asmid, stats in sorted(records.items()):
            writer.writerow({"ASMID": asmid, **stats})

    print(f"[INFO] Written to {output_path}")


if __name__ == "__main__":
    main()
