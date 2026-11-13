#!/usr/bin/env python3
"""Summarize the number of predicted proteins per genome and per phylum.

Reads a funannotate output tree (one folder per genome, named like the pipeline's makeSampleTag:
<Species>_<Strain>), counts the proteins in predict_results/, and joins the folder name to a samples
CSV (ASMID, SPECIES, STRAIN, PHYLUM, ...) for the phylum.

Counting methods (--method):
  stats  (default) mRNA count from predict_results/<name>.stats.json (annotation.mRNA). Fast: one small
         file per genome. Falls back to counting FASTA records when the file is missing or unreadable.
  fasta  count '>' records in predict_results/<name>.proteins.fa[.gz]. Slow: reads every protein file.
  --verify N  also count the FASTA records for N randomly chosen genomes and report any difference
         from the stats.json value, so the fast method can be checked on a given tree.

Plain and .gz files are both read.

Example:
  summarize_protein_counts.py --annotation-dir .../genome_annotation.pre_rc6 \\
      --samples .../Fungi_BFD/samples.csv --out-prefix protein_counts_pre_rc6 --verify 200
Writes <prefix>.by_genome.tsv and <prefix>.by_phylum.tsv and prints the summary.
"""
import argparse
import csv
import glob
import gzip
import json
import os
import random
import re
import statistics
import sys
from concurrent.futures import ThreadPoolExecutor


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def clean_strain(raw):
    """Python copy of cleanStrain() in nextflow/modules/common/utils.nf."""
    s = (raw or "").strip()
    s = re.sub(r"['\"]", "", s)
    s = s.split(";")[0].strip().replace(":", " ")
    s = re.sub(r"^\s*\*+", "", s)
    s = re.sub(r"\*+\s*$", "", s)
    s = re.sub(r"\s*\*+\s*", "-", s)
    return s.strip()


def make_sample_tag(species, strain):
    """Python copy of makeSampleTag() in nextflow/modules/common/utils.nf."""
    sp = re.sub(r"['\"]", "", (species or "").strip())
    tag = "_".join(x for x in (sp, clean_strain(strain)) if x)
    return re.sub(r"[\s/#\[\]?{}]+", "_", tag)


def load_samples(path):
    """Map folder tag -> row. Returns (map, number of tags that appear on more than one row)."""
    by_tag, dup = {}, 0
    with open_text(path) as fh:
        for row in csv.DictReader(fh):
            tag = make_sample_tag(row.get("SPECIES"), row.get("STRAIN"))
            if tag in by_tag:
                dup += 1
                continue
            by_tag[tag] = row
    return by_tag, dup


def find_one(folder, patterns):
    for pat in patterns:
        hits = sorted(glob.glob(os.path.join(folder, pat)))
        if hits:
            return hits[0]
    return None


def count_fasta(path):
    n = 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
    return n


def count_stats(path):
    """Return annotation.mRNA from a funannotate stats.json, or None if it cannot be read."""
    try:
        with open_text(path) as fh:
            return int(json.load(fh)["annotation"]["mRNA"])
    except (OSError, ValueError, KeyError, TypeError):
        return None


def count_genome(task):
    folder, method = task
    pr = os.path.join(folder, "predict_results")
    if not os.path.isdir(pr):
        return None, "no_predict_results"
    fasta = find_one(pr, ["*.proteins.fa", "*.proteins.fa.gz"])
    if method in ("stats", "auto"):
        stats = find_one(pr, ["*.stats.json", "*.stats.json.gz"])
        n = count_stats(stats) if stats else None
        if n is not None:
            return n, "stats"
    if fasta:
        return count_fasta(fasta), "fasta"
    return None, "no_proteins"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir", required=True, help="folder with one sub-folder per genome")
    ap.add_argument("--samples", required=True, help="samples CSV with SPECIES, STRAIN and PHYLUM columns")
    ap.add_argument("--method", choices=("stats", "fasta"), default="stats")
    ap.add_argument("--rank", default="PHYLUM", help="samples column to group by (default PHYLUM)")
    ap.add_argument("--verify", type=int, default=0, metavar="N", help="cross-check N random genomes by FASTA count")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out-prefix", default="protein_counts")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    by_tag, dup_tags = load_samples(args.samples)
    folders = sorted(d for d in os.listdir(args.annotation_dir)
                     if os.path.isdir(os.path.join(args.annotation_dir, d)))
    tasks = [(os.path.join(args.annotation_dir, d), args.method) for d in folders]
    with ThreadPoolExecutor(max_workers=args.workers) as ex:
        results = list(ex.map(count_genome, tasks))

    rows, status = [], {}
    for name, (n, src) in zip(folders, results):
        meta = by_tag.get(name)
        group = ((meta.get(args.rank) or "").strip() or "Unknown") if meta else "Unmapped"
        status[src] = status.get(src, 0) + 1
        rows.append({"folder": name, "ASMID": meta["ASMID"] if meta else "", args.rank: group,
                     "SPECIES": (meta.get("SPECIES") or "") if meta else "",
                     "proteins": "" if n is None else n, "source": src})

    with open(args.out_prefix + ".by_genome.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    counted = [r for r in rows if r["proteins"] != ""]
    groups = {}
    for r in counted:
        groups.setdefault(r[args.rank], []).append(r)
    summary = []
    for g, rs in groups.items():
        v = [r["proteins"] for r in rs]
        summary.append({args.rank: g, "genomes": len(v), "species": len({r["SPECIES"] for r in rs if r["SPECIES"]}),
                        "proteins": sum(v), "mean_per_genome": round(sum(v) / len(v), 1),
                        "median_per_genome": int(statistics.median(v)), "min": min(v), "max": max(v)})
    summary.sort(key=lambda x: -x["proteins"])
    with open(args.out_prefix + ".by_phylum.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(summary[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(summary)

    total = sum(r["proteins"] for r in counted)
    print(f"annotation dir : {args.annotation_dir}")
    print(f"genome folders : {len(folders)}; with a protein count: {len(counted)}")
    print("count source   : " + ", ".join(f"{k}={v}" for k, v in sorted(status.items())))
    print(f"unmapped to samples CSV (no {args.rank}): {sum(1 for r in rows if r[args.rank] == 'Unmapped')}"
          f"  (CSV rows sharing a folder tag, ignored: {dup_tags})")
    print(f"\nTOTAL predicted proteins: {total:,} in {len(counted):,} genomes\n")
    print(f"{args.rank:28s} {'genomes':>8s} {'species':>8s} {'proteins':>13s} {'%':>6s} {'mean':>8s} {'median':>7s}")
    for s in summary:
        print(f"{s[args.rank]:28s} {s['genomes']:8,d} {s['species']:8,d} {s['proteins']:13,d} "
              f"{100 * s['proteins'] / total:6.1f} {s['mean_per_genome']:8,.0f} {s['median_per_genome']:7,d}")

    if args.verify:
        random.seed(args.seed)
        pool = [(r["folder"], r["proteins"]) for r in counted if r["source"] == "stats"]
        sample = random.sample(pool, min(args.verify, len(pool)))

        def check(item):
            folder, n = item
            fasta = find_one(os.path.join(args.annotation_dir, folder, "predict_results"),
                             ["*.proteins.fa", "*.proteins.fa.gz"])
            return folder, n, (count_fasta(fasta) if fasta else None)

        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            checked = list(ex.map(check, sample))
        bad = [c for c in checked if c[2] != c[1]]
        print(f"\nverify: {len(checked)} genomes, stats.json vs FASTA records: {len(checked) - len(bad)} equal, "
              f"{len(bad)} different")
        for folder, n, f in bad[:15]:
            print(f"  {folder}: stats.json {n} vs FASTA {f}")


if __name__ == "__main__":
    sys.exit(main())
