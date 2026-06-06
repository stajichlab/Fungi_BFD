#!/usr/bin/env python3
"""
Check which genome_annotation folders need (re-)training because:
  1. Non-zero rnaseq_reads/{species_tag}_norm_R1.fastq.gz exists but the
     genome_annotation/{out}/training/ directory is absent, OR
  2. The training/ directory is older than the rnaseq R1 file (new reads
     were generated after training ran).

Species-to-folder mapping is read directly from samples.csv so that
ambiguous prefix matches (e.g. 'Aspergillus sp.' vs 'Aspergillus sp. 2663')
are resolved correctly.

Usage:
    python scripts/check_rnaseq_training.py [--project-dir DIR] [--report-ok] [--tsv]

Output (stdout):
    Tab-separated: out_dir  species_tag  status  detail
    Status values: MISSING_TRAINING | NEW_READS | OK (only with --report-ok)

Exit code is the number of folders that need action (0 = nothing to do).
"""

import argparse
import csv
import os
import re
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--project-dir", default=".", help="Project root (default: CWD)")
    p.add_argument("--report-ok", action="store_true", help="Also print OK entries")
    p.add_argument("--tsv", action="store_true", help="Print TSV header line")
    return p.parse_args()


def make_out(species, strain):
    """Replicate Nextflow: [species, strain].findAll{it}.join('_').replaceAll(spaces,'_').replaceAll(specials,'_')"""
    parts = [x for x in [species, strain] if x]
    s = "_".join(parts)
    s = re.sub(r"\s+", "_", s)
    s = re.sub(r"[\[\]*?{}]", "_", s)
    return s


def make_species_tag(species):
    return re.sub(r"\s+", "_", species)


def load_samples(samples_csv):
    """Return list of (out, species_tag) pairs from samples.csv."""
    result = []
    with open(samples_csv, newline="") as f:
        for row in csv.DictReader(f):
            species = row.get("SPECIES", "").strip().strip("'\"")
            strain = row.get("STRAIN", "").strip().strip("'\"")
            if not species:
                continue
            out = make_out(species, strain)
            tag = make_species_tag(species)
            result.append((out, tag))
    return result


def rnaseq_info(rnaseq_dir, species_tag):
    """Return (path, mtime) for a non-zero R1 file, or None if absent/empty."""
    path = os.path.join(rnaseq_dir, f"{species_tag}_norm_R1.fastq.gz")
    if not os.path.exists(path):
        return None
    if os.path.getsize(path) == 0:
        return None
    return (path, os.path.getmtime(path))


def _fmt_ts(ts):
    import datetime
    return datetime.datetime.fromtimestamp(ts).strftime("%Y-%m-%d %H:%M")


def main():
    args = parse_args()
    project = os.path.realpath(args.project_dir)

    samples_csv = os.path.join(project, "samples.csv")
    rnaseq_dir = os.path.join(project, "rnaseq_reads")
    ann_dir = os.path.join(project, "genome_annotation")

    for path, label in [(samples_csv, "samples.csv"), (rnaseq_dir, "rnaseq_reads/"), (ann_dir, "genome_annotation/")]:
        if not os.path.exists(path):
            print(f"[ERROR] Not found: {path}", file=sys.stderr)
            sys.exit(1)

    entries = load_samples(samples_csv)

    needs_action = 0
    rows = []

    for out, tag in sorted(set(entries)):  # deduplicate in case of repeated rows
        rna = rnaseq_info(rnaseq_dir, tag)
        if rna is None:
            continue  # no usable rnaseq reads for this species

        folder_path = os.path.join(ann_dir, out)
        if not os.path.isdir(folder_path):
            continue  # genome_annotation folder doesn't exist yet

        training_dir = os.path.join(folder_path, "training")
        r1_path, r1_mtime = rna

        if not os.path.isdir(training_dir):
            status = "MISSING_TRAINING"
            detail = "training/ dir absent"
            needs_action += 1
        else:
            train_mtime = os.path.getmtime(training_dir)
            if r1_mtime > train_mtime:
                status = "NEW_READS"
                detail = f"rnaseq {_fmt_ts(r1_mtime)} > training {_fmt_ts(train_mtime)}"
                needs_action += 1
            else:
                status = "OK"
                detail = ""

        if status != "OK" or args.report_ok:
            rows.append((out, tag, status, detail))

    if args.tsv:
        print("out_dir\tspecies_tag\tstatus\tdetail")
    for row in rows:
        print("\t".join(row))

    return needs_action


if __name__ == "__main__":
    sys.exit(min(main(), 125))
