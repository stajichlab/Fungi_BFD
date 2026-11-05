#!/usr/bin/env python3
"""Scatter plots of gene structure per genome, colored by phylum (top N phyla).

For each genome folder in a funannotate tree, reads predict_results/<name>.stats.json:
  gene count      annotation.mRNA              (protein-coding transcripts = predicted proteins)
  avg exon length annotation["transcript-level"].avg_exon_length, or the misspelled avg_exon_lenth
                  when the first is 0 (funannotate writes either key; the value is in whichever is nonzero).
                  funannotate averages only the exons of multi-exon transcripts; single-exon genes are left out.
  avg gene length annotation.avg_gene_length   (gene start to end in bp, introns included; mean over every gene
                  record in the GFF3, so tRNA and other non-coding genes are in it too)
  avg protein len annotation["transcript-level"].avg_protein_length  (amino acids, stop codon not counted). Times 3
                  plus 3 it equals the mean total CDS length per transcript (sum of CDS segments, stop included)
Joins the folder to the samples CSV for the phylum (see summarize_protein_counts.py), keeps the top N phyla
ranked by total predicted proteins (rows with no phylum are dropped), and writes:
  <prefix>.tsv                          one row per plotted genome
  <prefix>.avg_exon_length.png/.pdf     avg exon length vs gene count
  <prefix>.avg_gene_length.png/.pdf     avg gene length vs gene count
  <prefix>.avg_protein_length.png/.pdf  avg protein length (aa) vs gene count

Example:
  plot_gene_structure_by_phylum.py --annotation-dir .../genome_annotation.pre_rc6 \\
      --samples .../Fungi_BFD/samples.csv --out-prefix gene_structure_pre_rc6
"""
import argparse
import csv
import gzip
import json
import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from summarize_protein_counts import find_one, load_samples  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def read_stats(folder):
    pr = os.path.join(folder, "predict_results")
    path = find_one(pr, ["*.stats.json", "*.stats.json.gz"]) if os.path.isdir(pr) else None
    if not path:
        return None
    try:
        opener = gzip.open if path.endswith(".gz") else open
        with opener(path, "rt") as fh:
            ann = json.load(fh)["annotation"]
        tl = ann["transcript-level"]
        exon = tl.get("avg_exon_length") or tl.get("avg_exon_lenth") or 0
        aa = tl.get("avg_protein_length")
        return {"genes": int(ann["mRNA"]), "avg_exon_length": float(exon),
                "avg_gene_length": float(ann.get("avg_gene_length") or 0),
                "avg_protein_length": float(aa) if aa else 0.0}
    except (OSError, ValueError, KeyError, TypeError):
        return None


def scatter(rows, phyla, colors, ycol, ylabel, path_base, title, logy=False, logx=True):
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for ph in phyla:
        pts = [r for r in rows if r["PHYLUM"] == ph and r[ycol] > 0]
        ax.scatter([r["genes"] for r in pts], [r[ycol] for r in pts], s=7, alpha=0.4, linewidths=0,
                   color=colors[ph], label=f"{ph} (n={len(pts):,})")
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_xlabel("Predicted protein-coding genes per genome" + (" (log scale)" if logx else ""))
    ax.set_ylabel(ylabel + (" (log scale)" if logy else ""))
    ax.set_title(title)
    ax.grid(alpha=0.25)
    lg = ax.legend(markerscale=3, frameon=False, fontsize=9)
    for h in lg.legend_handles:
        h.set_alpha(1)
    fig.tight_layout()
    fig.savefig(path_base + ".png", dpi=200)
    fig.savefig(path_base + ".pdf")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--annotation-dir")
    ap.add_argument("--samples")
    ap.add_argument("--from-tsv", help="replot from a <prefix>.tsv written by an earlier run (skips reading stats.json)")
    ap.add_argument("--max-genes", type=int, default=0, metavar="N",
                    help="draw genomes with more than N genes on a separate outlier plot (and table) instead of the main plots; 0 = keep all")
    ap.add_argument("--linear", action="store_true", help="linear axes (default: log gene count; log gene length)")
    ap.add_argument("--top", type=int, default=5, help="number of phyla to plot (default 5)")
    ap.add_argument("--rank", default="PHYLUM")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out-prefix", default="gene_structure")
    ap.add_argument("--title-suffix", default="", help="text added to the plot titles, e.g. the tree name")
    args = ap.parse_args()

    if args.from_tsv:
        with open(args.from_tsv) as fh:
            rows = [{"folder": r["folder"], "PHYLUM": r["PHYLUM"], "genes": int(r["genes"]),
                     "avg_exon_length": float(r["avg_exon_length"]), "avg_gene_length": float(r["avg_gene_length"]),
                     "avg_protein_length": float(r.get("avg_protein_length") or 0)}
                    for r in csv.DictReader(fh, delimiter="\t")]
        if rows and not any(r["avg_protein_length"] for r in rows):
            ap.error("--from-tsv file has no avg_protein_length column (older output used avg_cds_length); rerun without --from-tsv")
        totals = {}
        for r in rows:
            totals[r["PHYLUM"]] = totals.get(r["PHYLUM"], 0) + r["genes"]
        phyla = [p for p, _ in sorted(totals.items(), key=lambda x: -x[1])[:args.top]]
        no_stats = unmapped = 0
    else:
        if not (args.annotation_dir and args.samples):
            ap.error("--annotation-dir and --samples are required without --from-tsv")
        by_tag, _ = load_samples(args.samples)
        folders = sorted(d for d in os.listdir(args.annotation_dir)
                         if os.path.isdir(os.path.join(args.annotation_dir, d)))
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            stats = list(ex.map(read_stats, [os.path.join(args.annotation_dir, d) for d in folders]))
        rows, no_stats, unmapped = [], 0, 0
        for name, st in zip(folders, stats):
            if st is None:
                no_stats += 1
                continue
            meta = by_tag.get(name)
            ph = ((meta.get(args.rank) or "").strip() if meta else "")
            if not ph:
                unmapped += 1
                continue
            rows.append({"folder": name, "PHYLUM": ph, **st})
        totals = {}
        for r in rows:
            totals[r["PHYLUM"]] = totals.get(r["PHYLUM"], 0) + r["genes"]
        phyla = [p for p, _ in sorted(totals.items(), key=lambda x: -x[1])[:args.top]]
    rows = [r for r in rows if r["PHYLUM"] in phyla]
    cmap = plt.get_cmap("tab10")
    colors = {p: cmap(i) for i, p in enumerate(phyla)}

    if not args.from_tsv:
        with open(args.out_prefix + ".tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["folder", "PHYLUM", "genes", "avg_exon_length", "avg_gene_length", "avg_protein_length"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    sfx = f" ({args.title_suffix})" if args.title_suffix else ""
    main_rows, out_rows = rows, []
    if args.max_genes:
        main_rows = [r for r in rows if r["genes"] <= args.max_genes]
        out_rows = [r for r in rows if r["genes"] > args.max_genes]
        with open(args.out_prefix + ".outliers.tsv", "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["folder", "PHYLUM", "genes", "avg_exon_length", "avg_gene_length", "avg_protein_length"],
                               delimiter="\t")
            w.writeheader()
            w.writerows(sorted(out_rows, key=lambda r: -r["genes"]))
        limit = f" (genomes with <= {args.max_genes:,} genes)"
    else:
        limit = ""
    scatter(main_rows, phyla, colors, "avg_exon_length", "Average exon length (bp)",
            args.out_prefix + ".avg_exon_length", "Average exon length vs gene count" + limit + sfx,
            logx=not args.linear)
    scatter(main_rows, phyla, colors, "avg_gene_length", "Average gene length (bp)",
            args.out_prefix + ".avg_gene_length", "Average gene length vs gene count" + limit + sfx,
            logy=not args.linear, logx=not args.linear)
    scatter(main_rows, phyla, colors, "avg_protein_length", "Average protein length (aa)",
            args.out_prefix + ".avg_protein_length", "Average protein length vs gene count" + limit + sfx,
            logy=not args.linear, logx=not args.linear)
    if out_rows:
        fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))
        for ax, ycol, yl in ((axes[0], "avg_exon_length", "Average exon length (bp)"),
                             (axes[1], "avg_gene_length", "Average gene length (bp)"),
                             (axes[2], "avg_protein_length", "Average protein length (aa)")):
            for ph in phyla:
                pts = [r for r in out_rows if r["PHYLUM"] == ph and r[ycol] > 0]
                if pts:
                    ax.scatter([r["genes"] for r in pts], [r[ycol] for r in pts], s=25, alpha=0.8, linewidths=0,
                               color=colors[ph], label=f"{ph} (n={len(pts)})")
            ax.set_xlabel("Predicted protein-coding genes per genome")
            ax.set_ylabel(yl)
            ax.grid(alpha=0.25)
        axes[0].legend(frameon=False, fontsize=9)
        fig.suptitle(f"Genomes with more than {args.max_genes:,} predicted genes (n={len(out_rows)})" + sfx)
        fig.tight_layout()
        fig.savefig(args.out_prefix + ".outliers.png", dpi=200)
        fig.savefig(args.out_prefix + ".outliers.pdf")
        plt.close(fig)

    if out_rows:
        print(f"outlier plot: {len(out_rows)} genomes above {args.max_genes:,} genes (main plots show {len(main_rows):,})")
    print(f"genomes plotted: {len(rows):,} in top {len(phyla)} {args.rank.lower()}: {', '.join(phyla)}")
    print(f"skipped: {no_stats:,} without a readable stats.json; {unmapped:,} with no {args.rank} in the samples CSV")
    for p in phyla:
        rs = [r for r in rows if r["PHYLUM"] == p]
        zero_e = sum(1 for r in rs if r["avg_exon_length"] <= 0)
        zero_g = sum(1 for r in rs if r["avg_gene_length"] <= 0)
        zero_p = sum(1 for r in rs if r["avg_protein_length"] <= 0)
        print(f"  {p}: {len(rs):,} genomes ({zero_e} with no avg exon length, {zero_g} with no avg gene length, "
              f"{zero_p} with no avg protein length; not plotted in that panel)")


if __name__ == "__main__":
    sys.exit(main())
