#!/usr/bin/env python3
"""Collect predict-arm scores into one table: predict_arms/scorecard.tsv.

One row per <genome>/<arm>.B.<evid>: gffcompare holdout scores (score.tsv from
predict_scorer.py) + BUSCO C/S/D/F/M (busco_short_summary.txt) + training-set
size from step A (final_training_models.gff3 gene count) + wall time.
"""
import csv, glob, os, re

A = os.path.dirname(os.path.abspath(__file__))
rows = []
for d in sorted(glob.glob(f"{A}/*/*.B.*")):
    genome, run = d.split("/")[-2], d.split("/")[-1]
    arm = run.split(".B.")[0]
    r = {"genome": genome, "arm": arm, "evidence": run.split(".B.")[1]}
    st = f"{d}/run_status.tsv"
    if not os.path.exists(st):
        continue
    s = dict(l.rstrip("\n").split("\t") for l in open(st))
    r["exit"] = s.get("exit_status", "")
    sc = f"{d}/score.tsv"
    if os.path.exists(sc):
        rec = list(csv.DictReader(open(sc), delimiter="\t"))
        if rec:
            for k, v in rec[-1].items():
                if k not in ("genome", "variant"):
                    r[k] = v
    b = f"{d}/busco_short_summary.txt"
    if os.path.exists(b):
        m = re.search(r"C:([\d.]+)%\[S:([\d.]+)%,D:([\d.]+)%\],F:([\d.]+)%,M:([\d.]+)%", open(b).read())
        if m:
            r.update(busco_C=m[1], busco_D=m[3], busco_F=m[4], busco_M=m[5])
    t = f"{A}/{genome}/{arm}.A/out/predict_misc/final_training_models.gff3"
    if os.path.exists(t):
        r["train_models"] = sum(1 for l in open(t) if "\tgene\t" in l)
    rows.append(r)
cols = []
for r in rows:
    for k in r:
        if k not in cols:
            cols.append(k)
with open(f"{A}/scorecard.tsv", "w") as f:
    w = csv.DictWriter(f, cols, delimiter="\t", restval="")
    w.writeheader()
    w.writerows(rows)
print(f"{len(rows)} rows -> {A}/scorecard.tsv")
