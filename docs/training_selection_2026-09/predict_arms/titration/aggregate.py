#!/usr/bin/env python3
"""Aggregate experiment A titration rows into titration_scores.tsv.

Adds BUSCO-forced comparator rows per genome (same step-B fixed evidence, same
predict_scorer.py, same fungi_odb10 BUSCO): N="busco" draws 1-3 from
busco4_r{1,2,3}.B.fixed (code_new4), and N="busco_code_new" from busco.B.fixed.
Rows whose step A or step B failed are kept with status=failed and empty metrics.
Also writes titration_summary.tsv: mean, SD and count per genome x N (ok rows only).
"""
import csv, glob, os, re, statistics as st

T = os.path.dirname(os.path.abspath(__file__))
A = os.path.dirname(T)
COLS = ["genome", "N", "draw", "status", "n_subset", "n_complete_used", "n_keepers", "n_single_admitted",
        "locus_sn", "locus_pr", "exon_sn", "exon_pr", "intron_chain_sn", "intron_chain_pr", "pred_genes",
        "busco_c", "exit_stepA", "exit_stepB", "runtime_s", "node", "source"]
METRICS = ["locus_sn", "locus_pr", "exon_sn", "exon_pr", "intron_chain_sn", "intron_chain_pr", "busco_c",
           "n_complete_used", "n_keepers", "n_single_admitted", "pred_genes"]

rows = []
for f in sorted(glob.glob(f"{T}/rows/*.tsv")):
    r = list(csv.DictReader(open(f), delimiter="\t"))[-1]
    r["source"] = "titration"
    r["status"] = "ok" if r.get("exit_stepA") == "0" and r.get("exit_stepB") == "0" and r.get("locus_sn") else "failed"
    rows.append(r)

def comparator(g, d, n, draw):
    """One BUSCO-forced step-B row from predict_arms/<g>/<d> (score.tsv + BUSCO summary)."""
    sc = f"{A}/{g}/{d}/score.tsv"
    if not os.path.exists(sc):
        return None
    rec = list(csv.DictReader(open(sc), delimiter="\t"))[-1]
    r = {"genome": g, "N": n, "draw": str(draw), "status": "ok", "source": d}
    r.update({c: rec.get(c, "") for c in METRICS if c in rec})
    b = f"{A}/{g}/{d}/busco_short_summary.txt"
    m = re.search(r"C:([\d.]+)%", open(b).read()) if os.path.exists(b) else None
    r["busco_c"] = m.group(1) if m else ""
    return r


# N="busco": BUSCO-forced repeats run with code_new4, the titration step-A code (D96).
# N="busco_code_new": the earlier single busco.B.fixed run (code_new snapshot), kept for reference.
for g in sorted({r["genome"] for r in rows}):
    for k in (1, 2, 3):
        r = comparator(g, f"busco4_r{k}.B.fixed", "busco", k)
        if r:
            rows.append(r)
    r = comparator(g, "busco.B.fixed", "busco_code_new", 1)
    if r:
        rows.append(r)
    else:
        print("WARNING: no busco.B.fixed comparator for", g)

nkey = lambda v: (1, v) if not v.isdigit() else (0, int(v))
rows.sort(key=lambda r: (r["genome"], nkey(r["N"]), int(r["draw"])))
with open(f"{T}/titration_scores.tsv", "w") as f:
    w = csv.DictWriter(f, COLS, delimiter="\t", extrasaction="ignore", restval="")
    w.writeheader()
    w.writerows(rows)

groups = {}
for r in rows:
    groups.setdefault((r["genome"], r["N"]), []).append(r)
with open(f"{T}/titration_summary.tsv", "w") as f:
    hdr = ["genome", "N", "n_ok", "n_failed"] + [f"{m}_{s}" for m in METRICS for s in ("mean", "sd")]
    f.write("\t".join(hdr) + "\n")
    for (g, n), rs in sorted(groups.items(), key=lambda kv: (kv[0][0], nkey(kv[0][1]))):
        ok = [r for r in rs if r["status"] == "ok"]
        out = [g, n, str(len(ok)), str(len(rs) - len(ok))]
        for m in METRICS:
            v = [float(r[m]) for r in ok if r.get(m) not in (None, "")]
            out += [f"{st.mean(v):.2f}" if v else "", f"{st.stdev(v):.2f}" if len(v) > 1 else ""]
        f.write("\t".join(out) + "\n")
print(f"rows={len(rows)} ok={sum(r['status']=='ok' for r in rows)} failed={sum(r['status']=='failed' for r in rows)}")
