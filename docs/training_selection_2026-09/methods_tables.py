#!/usr/bin/env python3
"""Generate the numeric tables for the training-data-selection methods write-up
directly from the result files, so no number is transcribed by hand.

Output: methods_tables.md (Markdown tables), next to this script.
Sources:
  predict_arms/scorecard.tsv          holdout gffcompare scores per arm
  predict_arms/single_exon_scores.tsv exact single/multi-exon matches (recomputed here)
  refseq_benchmark/benchmark.tsv      training-set / getBestModel vs RefSeq
  refseq_benchmark/rank_benchmark.tsv getBestModel ranking variants vs RefSeq
"""
import csv, os, subprocess

X = os.path.dirname(os.path.abspath(__file__))
A = f"{X}/predict_arms"
subprocess.run(["python3", f"{A}/scorecard.py"], check=True, stdout=subprocess.DEVNULL)
SC = list(csv.DictReader(open(f"{A}/scorecard.tsv"), delimiter="\t"))
NAME = {"Neurospora_crassa_OR74A": "N. crassa OR74A", "Aspergillus_nidulans_FGSC_A4": "A. nidulans FGSC A4",
        "Botrytis_cinerea_B05.10": "B. cinerea B05.10", "Cryptococcus_neoformans_H99": "C. neoformans H99",
        "Schizophyllum_commune_H4-8": "S. commune H4-8"}
ORDER = list(NAME)


def get(g, arm, ev):
    for r in SC:
        if r["genome"] == g and r["arm"] == arm and r["evidence"] == ev and r.get("exit") == "0" and r.get("locus_sn"):
            return r
    return None


def fmt(r, a, b):
    return f"{r[a]} / {r[b]}" if r and r.get(a) else "n/a"


out = []
w = out.append

w("## Table M1. Locus rule for one PASA model per locus (step B with each rule's own PASA evidence)\n")
w("Holdout-chromosome gffcompare, CDS level. Values are locus Sn / Pr (%), predicted genes.\n")
w("| Genome | tx_strand | cds_strand | cds_blind |")
w("|---|---|---|---|")
for g in ORDER:
    cells = []
    for arm in ("tx", "cdsS", "cdsB"):
        r = get(g, arm, "own")
        cells.append(f"{fmt(r, 'locus_sn', 'locus_pr')} ({r['pred_genes']})" if r else "n/a")
    w(f"| {NAME[g]} | " + " | ".join(cells) + " |")

w("\n## Table M2. Training-set construction, same fixed EVM evidence (training effect only)\n")
w("Locus Sn / Pr (%); training models in parentheses. old = funannotate 41a2fd7; oldnew = old getBestModel + new selection; "
  "tx = new getBestModel (tx_strand, complete-first) + new selection; busco = BUSCO-forced training.\n")
w("| Genome | old | oldnew | tx | busco |")
w("|---|---|---|---|---|")
for g in ORDER:
    cells = []
    for arm in ("old", "oldnew", "tx", "busco"):
        r = get(g, arm, "fixed")
        tm = f" ({r['train_models']})" if r and r.get("train_models") else ""
        cells.append(fmt(r, "locus_sn", "locus_pr") + tm if r else "n/a")
    w(f"| {NAME[g]} | " + " | ".join(cells) + " |")

w("\n## Table M3. Ranking inside a locus: complete-first vs guarded (complete only if CDS ≥ 80% of the locus's longest)\n")
w("| Genome | Evidence | complete-first locus Sn / Pr (train models) | guarded locus Sn / Pr (train models) |")
w("|---|---|---|---|")
for g in ORDER:
    for ev in ("fixed", "own"):
        a, b = get(g, "tx", ev), get(g, "txG", ev)
        if a and b:
            w(f"| {NAME[g]} | {ev} | {fmt(a, 'locus_sn', 'locus_pr')} ({a.get('train_models','')}) | "
              f"{fmt(b, 'locus_sn', 'locus_pr')} ({b.get('train_models','')}) |")

# single-exon table from exact matching
subprocess.run(["python3", f"{A}/single_exon_score.py", "tx2", "se", "busco", "old", "tx", "txR1R2"],
               check=True, stdout=subprocess.DEVNULL)
SE = list(csv.DictReader(open(f"{A}/single_exon_scores.tsv"), delimiter="\t"))


def se(g, arm):
    for r in SE:
        if r["genome"] == g and r["arm"] == arm:
            return r
    return None


w("\n## Table M4. Single-exon training genes (R6 option b): se (on) vs tx2 (off), same code, fixed evidence\n")
w("Exact CDS-chain match to RefSeq protein-coding mRNAs on holdout chromosomes. Sn per RefSeq gene, Pr per predicted model.\n")
w("| Genome | RefSeq single / multi genes | Single Sn | Single Pr | Multi Sn | Multi Pr |")
w("|---|---|---|---|---|---|")
for g in ORDER:
    a, b = se(g, "tx2"), se(g, "se")
    if a and b:
        w(f"| {NAME[g]} | {a['ref_single_genes']} / {a['ref_multi_genes']} | {a['single_sn']} → {b['single_sn']} | "
          f"{a['single_pr']} → {b['single_pr']} | {a['multi_sn']} → {b['multi_sn']} | {a['multi_pr']} → {b['multi_pr']} |")

w("\n## Table M5. PASA training vs BUSCO training, divergent reads (N. crassa OR74A; RNA-seq from strain HJDF, median read identity 96.7%)\n")
w("| Arm | PASA input | Single-exon training | Training | Fixed evidence locus Sn / Pr | Own evidence locus Sn / Pr |")
w("|---|---|---|---|---|---|")
g = "Neurospora_crassa_OR74A"
for arm, pasa, sx, tr in [("tx", "rc1 (old minimap2 parser)", "no", "PASA"), ("se", "rc1", "yes", "PASA"),
                          ("txR1", "R1", "no", "PASA"), ("txR1R2", "R1 + R2", "no", "PASA"),
                          ("txID90", "R1 + gmap, 90% identity", "no", "PASA"),
                          ("txRel", "R1, relaxed validation", "no", "PASA"),
                          ("seR1R2", "R1 + R2", "yes", "PASA"), ("busco", "rc1 (evidence only)", "n/a", "BUSCO"),
                          ("buscoR1R2", "R1 + R2 (evidence only)", "n/a", "BUSCO")]:
    f, o = get(g, arm, "fixed"), get(g, arm, "own")
    w(f"| {arm} | {pasa} | {sx} | {tr} | {fmt(f, 'locus_sn', 'locus_pr')} | {fmt(o, 'locus_sn', 'locus_pr')} |")

w("\n## Table M6. Genome with few complete PASA models (S. commune H4-8: 93 complete of 822 PASA models, 30 keepers; divergent reads)\n")
w("| Arm | Locus Sn / Pr | Intron-chain Sn / Pr | Exon Sn / Pr | Proteome BUSCO C (%) |")
w("|---|---|---|---|---|")
g = "Schizophyllum_commune_H4-8"
for arm in ("busco", "old", "txR1R2", "tx"):
    r = get(g, arm, "fixed")
    if r:
        w(f"| {arm} | {fmt(r,'locus_sn','locus_pr')} | {fmt(r,'intron_chain_sn','intron_chain_pr')} | "
          f"{fmt(r,'exon_sn','exon_pr')} | {r.get('busco_C','')} |")

# benchmark tables
B = list(csv.DictReader(open(f"{X}/refseq_benchmark/benchmark.tsv"), delimiter="\t"))
w("\n## Table M7. Training set against RefSeq: old pipeline vs new selection (no filterGeneMark keeper set; upper bound, see text)\n")
w("| Genome | Old: models, exact % | New: models, exact % | Exact RefSeq genes old → new | Redundant old → new |")
w("|---|---|---|---|---|")
for g in ORDER:
    o = [r for r in B if r["genome"] == g and r["variant"] == "old+old_select"]
    n = [r for r in B if r["genome"] == g and r["variant"] == "tx_strand+new_select"]
    if o and n:
        o, n = o[0], n[0]
        w(f"| {NAME[g]} | {o['models']}, {o['exact_pct']} | {n['models']}, {n['exact_pct']} | "
          f"{o['refseq_genes_exact']} → {n['refseq_genes_exact']} | {o['redundant_models']} → {n['redundant_models']} |")

RK = list(csv.DictReader(open(f"{X}/refseq_benchmark/rank_benchmark.tsv"), delimiter="\t"))
w("\n## Table M8. getBestModel ranking against RefSeq (PASA models passed to EVM): exact CDS / exact intron chain\n")
w("| PASA source | Genome | complete-first | structure-first | guarded |")
w("|---|---|---|---|---|")
for src in ("rc1", "R1", "R1R2"):
    for g in ORDER:
        c = {r["rank"]: r for r in RK if r["source"] == src and r["genome"] == g}
        if len(c) == 3:
            w(f"| {src} | {NAME[g]} | " + " | ".join(
                f"{c[k]['exact_cds']} / {c[k]['intron_chain_exact']}" for k in ("complete_first", "struct_first", "guarded")) + " |")

open(f"{X}/methods_tables.md", "w").write("\n".join(out) + "\n")
print(f"wrote {X}/methods_tables.md")
