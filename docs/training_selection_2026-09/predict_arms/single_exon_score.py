#!/usr/bin/env python3
"""Single-exon vs multi-exon accuracy of step B predictions on holdout chromosomes
(R6 option b evaluation, DECISIONS D38/D42/D45). Exact CDS-chain matching against
RefSeq protein-coding mRNAs (any isoform).

Usage: single_exon_score.py [arm ...]   (default: se tx2)
Writes predict_arms/single_exon_scores.tsv.
"""
import collections, glob, os, sys

A = os.path.dirname(os.path.abspath(__file__))
X = os.path.dirname(A)
ACC = {"Neurospora_crassa_OR74A": "GCF_000182925.2", "Aspergillus_nidulans_FGSC_A4": "GCF_000011425.1",
       "Botrytis_cinerea_B05.10": "GCF_000143535.2", "Cryptococcus_neoformans_H99": "GCF_000149245.1",
       "Schizophyllum_commune_H4-8": "GCF_000143185.2"}


def chains(gff, keep, coding_only):
    parent, cds, coding = {}, collections.defaultdict(list), set()
    for line in open(gff):
        if line.startswith("#"):
            continue
        c = line.rstrip("\n").split("\t")
        if len(c) < 9 or c[0] not in keep:
            continue
        a = dict(x.split("=", 1) for x in c[8].split(";") if "=" in x)
        if c[2] == "gene" and a.get("gene_biotype") == "protein_coding":
            coding.add(a["ID"])
        elif c[2] in ("mRNA", "transcript") and "ID" in a:
            parent[a["ID"]] = a.get("Parent")
        elif c[2] == "CDS" and "Parent" in a:
            cds[a["Parent"]].append((c[0], c[6], int(c[3]), int(c[4])))
    out = {}  # chain -> gene
    for m, parts in cds.items():
        g = parent.get(m, m)
        if coding_only and g not in coding:
            continue
        parts.sort(key=lambda x: x[2])
        out[(parts[0][0], parts[0][1], tuple((s, e) for _, _, s, e in parts))] = g
    return out


def main(arms):
    rows = []
    for gdir in sorted(glob.glob(f"{A}/*/")):
        genome = os.path.basename(gdir.rstrip("/"))
        if genome not in ACC:
            continue
        hold = set(open(f"{gdir}/split.holdout_chroms.txt").read().split())
        ref = chains(f"{X}/refseq/{ACC[genome]}.gff", hold, True)
        ref_s = {k for k in ref if len(k[2]) == 1}
        ref_m = {k for k in ref if len(k[2]) > 1}
        ref_s_genes = {ref[k] for k in ref_s}
        ref_m_genes = {ref[k] for k in ref_m}
        for arm in arms:
            d = f"{gdir}/{arm}.B.fixed"
            pred_files = [p for p in glob.glob(f"{d}/out/predict_results/*.gff3")]
            if not pred_files:
                continue
            pred = chains(pred_files[0], hold, False)
            ps = {k for k in pred if len(k[2]) == 1}
            pm = {k for k in pred if len(k[2]) > 1}
            r = {"genome": genome, "arm": arm,
                 "ref_single_genes": len(ref_s_genes), "pred_single": len(ps),
                 "single_sn": round(100.0 * len({ref[k] for k in ps & ref_s}) / max(len(ref_s_genes), 1), 1),
                 "single_pr": round(100.0 * len(ps & ref_s) / max(len(ps), 1), 1),
                 "ref_multi_genes": len(ref_m_genes), "pred_multi": len(pm),
                 "multi_sn": round(100.0 * len({ref[k] for k in pm & ref_m}) / max(len(ref_m_genes), 1), 1),
                 "multi_pr": round(100.0 * len(pm & ref_m) / max(len(pm), 1), 1)}
            rows.append(r)
    with open(f"{A}/single_exon_scores.tsv", "w") as f:
        f.write("\t".join(rows[0].keys()) + "\n")
        for r in rows:
            f.write("\t".join(str(v) for v in r.values()) + "\n")
    for r in rows:
        print("\t".join(str(v) for v in r.values()))


if __name__ == "__main__":
    main(sys.argv[1:] or ["tx2", "se"])
