#!/usr/bin/env python3
"""Audit which species' normalized RNA-seq reads can be safely reclaimed.

Safe-deletion predicate (per species_tag):
  1. rnaseq_data/<species_tag>.trinity-GG.fasta exists AND is non-empty, AND
  2. EVERY assembly mapped to that species_tag in samples.csv has a non-empty
     genome_annotation_training/<out>/training/funannotate_train.pasa.gff3
     (accepts .gff3 or .gff3.gz).

Reads are consumed twice: once by RNASEQ_PREPARE (-> Trinity, per species) and
once per strain by FUNANNOTATE_TRAIN (-> PASA/BAM). Both must be complete before
the reads in rnaseq_reads/<species_tag>_norm_{R1,R2,SE}.fastq.gz are disposable.

Read-only. Prints a classification summary and reclaimable byte totals.
"""
import argparse
import csv
import re
import sys
from pathlib import Path


# ── Faithful Python ports of nextflow/lib/SampleUtils.groovy ──────────────────
def clean_strain(raw):
    s = (raw or "").strip().replace("'", "").replace('"', "")
    s = s.split(";")[0].strip().replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)        # leading '*' removed
    s = re.sub(r"\*+\s*$", "", s)        # trailing '*' removed
    s = re.sub(r"\s*\*+\s*", "-", s)     # '*' between words -> '-'
    return s.strip()


def make_sample_tag(raw_species, raw_strain):
    sp = (raw_species or "").strip().replace("'", "").replace('"', "")
    st = clean_strain(raw_strain)
    joined = "_".join([x for x in (sp, st) if x])
    return re.sub(r"[\s/\#\[\]\?\{\}]+", "_", joined)


def species_tag_of(raw_species):
    # Matches funannotate.nf: species (quote-stripped) then \s+ -> '_'
    sp = (raw_species or "").strip().replace("'", "").replace('"', "")
    return re.sub(r"\s+", "_", sp)


def first_existing(*paths):
    for p in paths:
        if p.exists():
            return p
    return None


def human(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f}{unit}"
        n /= 1024


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=".", help="launchDir (default: cwd)")
    ap.add_argument("--samples", default="samples.csv")
    ap.add_argument("--reads-dir", default="rnaseq_reads")
    ap.add_argument("--rnaseq-data", default="rnaseq_data")
    ap.add_argument("--training", default="genome_annotation_training")
    ap.add_argument("--list-safe", action="store_true",
                    help="print every species_tag in the SAFE class")
    ap.add_argument("--list-blocked", action="store_true",
                    help="print Trinity-ready species blocked only by untrained strains")
    args = ap.parse_args()

    root = Path(args.root).resolve()
    samples = root / args.samples
    reads_dir = root / args.reads_dir
    data_dir = root / args.rnaseq_data
    train_dir = root / args.training

    if not samples.exists():
        sys.exit(f"ERROR: samples file not found: {samples}")

    # species_tag -> list of (out, asmid)
    species = {}
    with open(samples, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            sp_in = row.get("SPECIES") or ""
            strain = row.get("STRAIN") or ""
            asmid = (row.get("ASMID") or "").strip()
            if not sp_in.strip() or not asmid:
                continue
            stag = species_tag_of(sp_in)
            out = make_sample_tag(sp_in, strain)
            species.setdefault(stag, []).append((out, asmid))

    classes = {
        "SAFE": [],            # trinity non-empty + all strains trained + reads present
        "ALREADY_RECLAIMED": [],  # trinity ok, all trained, but reads already 0/absent
        "BLOCKED_UNTRAINED": [],  # trinity ok but >=1 strain missing pasa.gff3
        "NO_TRINITY": [],      # trinity missing
        "EMPTY_TRINITY": [],   # trinity exists but 0 bytes (no RNA-seq species)
    }
    reclaimable_bytes = 0

    for stag, members in sorted(species.items()):
        trinity = data_dir / f"{stag}.trinity-GG.fasta"
        r1 = reads_dir / f"{stag}_norm_R1.fastq.gz"
        r2 = reads_dir / f"{stag}_norm_R2.fastq.gz"
        se = reads_dir / f"{stag}_norm_SE.fastq.gz"
        read_files = [p for p in (r1, r2, se) if p.exists()]
        read_bytes = sum(p.stat().st_size for p in read_files)
        has_real_reads = any(p.stat().st_size > 0 for p in read_files)

        if not trinity.exists():
            classes["NO_TRINITY"].append((stag, read_bytes))
            continue
        if trinity.stat().st_size == 0:
            classes["EMPTY_TRINITY"].append((stag, read_bytes))
            continue

        # All strains trained?
        untrained = []
        for out, asmid in members:
            gff = first_existing(
                train_dir / out / "training" / "funannotate_train.pasa.gff3",
                train_dir / out / "training" / "funannotate_train.pasa.gff3.gz",
            )
            if gff is None or gff.stat().st_size == 0:
                untrained.append(out)

        if untrained:
            classes["BLOCKED_UNTRAINED"].append((stag, read_bytes, len(members), len(untrained)))
            continue

        if has_real_reads:
            classes["SAFE"].append((stag, read_bytes))
            reclaimable_bytes += read_bytes
        else:
            classes["ALREADY_RECLAIMED"].append((stag, read_bytes))

    # ── Report ───────────────────────────────────────────────────────────────
    total_species = len(species)
    print(f"\nAudit root: {root}")
    print(f"samples.csv species_tags: {total_species}")
    print("-" * 64)
    print(f"  SAFE to reclaim ............ {len(classes['SAFE']):5d}  "
          f"({human(reclaimable_bytes)} reclaimable)")
    print(f"  Already reclaimed/empty .... {len(classes['ALREADY_RECLAIMED']):5d}")
    print(f"  Blocked: untrained strains . {len(classes['BLOCKED_UNTRAINED']):5d}")
    print(f"  No Trinity yet ............. {len(classes['NO_TRINITY']):5d}")
    print(f"  Empty Trinity (no RNA-seq) . {len(classes['EMPTY_TRINITY']):5d}")
    print("-" * 64)

    blocked_bytes = sum(x[1] for x in classes["BLOCKED_UNTRAINED"])
    print(f"  (blocked-but-Trinity-ready hold {human(blocked_bytes)} pending training)")

    if args.list_safe:
        print("\n# SAFE species_tags (reclaimable bytes):")
        for stag, b in sorted(classes["SAFE"], key=lambda x: -x[1]):
            print(f"  {human(b):>10}  {stag}")
    if args.list_blocked:
        print("\n# BLOCKED (Trinity ready, strains still untrained):")
        for stag, b, n, u in sorted(classes["BLOCKED_UNTRAINED"], key=lambda x: -x[1]):
            print(f"  {human(b):>10}  {stag}  ({u}/{n} strains untrained)")


if __name__ == "__main__":
    main()
