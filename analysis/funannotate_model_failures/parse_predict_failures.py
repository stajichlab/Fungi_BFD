#!/usr/bin/env python3
"""
parse_predict_failures.py — Ground-truth detection of funannotate `predict`
failures, cross-tabulated against the asm_stats pre-screen.

Authoritative failure signal: funannotate emits, verbatim, when it cannot
assemble enough training genes for Augustus:

    Not enough gene models N to train Augustus (30 required), exiting

We scan every persisted funannotate-predict.log, classify each genome's
outcome, join it to taxonomy (samples.csv) and assembly metrics
(results/genome_stats/asm_stats/*.stats.txt), and emit:

  1. <out>/predict_outcomes.tsv      — one row per genome, all fields
  2. <out>/failed_genomes.tsv        — only non-completed runs
  3. <out>/crosstab_summary.txt      — outcome x asm_stats-risk contingency,
                                       i.e. "does the asm_stats screen predict
                                       the real failures?"

Join keys
  log line 1 carries `--name <LOCUSTAG>`; LOCUSTAG is unique in samples.csv,
  so LOCUSTAG -> ASMID -> {PHYLUM, SPECIES} and ASMID -> {total_bp, contigs, N50}.
  (The `-i` path is unreliable: gzipped genomes are inflated to genome_input.fa,
  losing the ASMID, so we never parse ASMID from the input path.)

Usage
  python3 parse_predict_failures.py \
      --log-roots genome_annotation do_annotation/genome_annotation \
      --samples   samples.csv \
      --asm-stats results/genome_stats/asm_stats \
      --outdir    analysis/funannotate_model_failures \
      [--phylum Basidiomycota]    # restrict crosstab to one phylum
"""
import argparse, csv, glob, os, re, sys
from collections import Counter

# funannotate's exact "too few training models" abort line.
RE_TOO_FEW   = re.compile(r"Not enough gene models (\d+) to train Augustus \((\d+) required\)")
RE_NAME      = re.compile(r"--name\s+(\S+)")
RE_BUSCO_VAL = re.compile(r"(\d+)\s+valid BUSCO predictions found")
RE_BUSCO_OK  = re.compile(r"(\d+)\s+BUSCO predictions validated")
RE_PASA_FALL = re.compile(r"Only (\d+) training models found from PASA, less than --min_training_models")
RE_SUMMARY   = re.compile(r"Summary of gene models:.*?'total':\s*(\d+)")
RE_AUG_FAIL  = re.compile(r"AUGUSTUS training failed")
RE_CMD_ERR   = re.compile(r"CMD ERROR")
RE_ZERO_BUSCO= re.compile(r"\b0 valid BUSCO predictions found")


def parse_stats(path):
    d = {}
    with open(path) as fh:
        for line in fh:
            if "=" in line:
                k, v = line.split("=", 1)
                d[k.strip()] = v.strip()
    return d


def asm_risk(total_bp, contigs, n50, nonN_bp):
    """Mirror of the pre-screen rule. Returns ('HIGH'|'MED'|'OK', reason)."""
    if total_bp is None:
        return ("NA", "no asm_stats")
    small = (nonN_bp if nonN_bp is not None else total_bp) < 15_000_000
    fragmented = (n50 is not None and n50 < 20_000) or (contigs is not None and contigs > 500)
    if small and fragmented:
        return ("HIGH", "small+fragmented")
    if fragmented:
        return ("MED", "fragmented_only")
    if small:
        return ("MED", "small_only")
    return ("OK", "")


def classify(text):
    """Return (status, n_training_models_or_None) from full log text."""
    m = RE_TOO_FEW.search(text)
    if m:
        return ("FAILED_too_few_models", int(m.group(1)))
    if RE_AUG_FAIL.search(text):
        return ("FAILED_augustus_training", None)
    if RE_SUMMARY.search(text):
        return ("COMPLETED", None)
    if RE_CMD_ERR.search(text):
        return ("FAILED_cmd_error", None)
    if RE_ZERO_BUSCO.search(text):
        return ("INCOMPLETE_zero_busco", None)
    return ("INCOMPLETE_unknown", None)


def last_meaningful_line(text):
    for line in reversed(text.splitlines()):
        s = line.strip()
        if not s:
            continue
        # skip tbl2asn validation spam and raw gene/exon dump lines
        if "out_of_bounds" in s or s.startswith("gene") or "errors(" in s:
            continue
        return s[:120]
    return ""


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log-roots", nargs="+", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--asm-stats", required=True, help="dir of <ASMID>.stats.txt")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--phylum", default=None, help="restrict crosstab to this phylum")
    args = ap.parse_args()

    # --- taxonomy: LOCUSTAG -> (asmid, phylum, species, genus) ---
    by_locustag = {}
    with open(args.samples, newline="") as fh:
        for row in csv.DictReader(fh):
            lt = (row.get("LOCUSTAG") or "").strip()
            if lt:
                by_locustag[lt] = (
                    row.get("ASMID", "").strip(),
                    row.get("PHYLUM", "").strip(),
                    (row.get("SPECIES_IN") or "").strip(),
                    (row.get("GENUS") or "").strip(),
                )

    # --- asm_stats: ASMID -> dict of metrics ---
    asm = {}
    for p in glob.glob(os.path.join(args.asm_stats, "*.stats.txt")):
        asmid = os.path.basename(p)[: -len(".stats.txt")]
        d = parse_stats(p)
        def gi(k):
            try:
                return int(d[k])
            except (KeyError, ValueError):
                return None
        total = gi("TOTAL LENGTH")
        nb = gi("TOTAL N BASES")
        asm[asmid] = dict(
            total_bp=total,
            contigs=gi("CONTIG COUNT"),
            n50=gi("N50"),
            nonN_bp=(total - nb) if (total is not None and nb is not None) else total,
        )

    # --- scan logs (dedupe by realpath; prefer most recent mtime per locustag) ---
    seen = set()
    best = {}  # locustag -> (mtime, record)
    n_logs = 0
    for root in args.log_roots:
        for log in glob.glob(os.path.join(root, "*", "logfiles", "funannotate-predict.log")):
            rp = os.path.realpath(log)
            if rp in seen:
                continue
            seen.add(rp)
            n_logs += 1
            try:
                with open(log, errors="replace") as fh:
                    text = fh.read()
            except OSError:
                continue
            first = text.split("\n", 1)[0]
            mname = RE_NAME.search(first)
            locustag = mname.group(1) if mname else None
            out_name = log.split(os.sep)[-3]
            status, n_models = classify(text)
            mb = RE_BUSCO_VAL.search(text)
            mo = RE_BUSCO_OK.search(text)
            mp = RE_PASA_FALL.search(text)
            ms = RE_SUMMARY.search(text)
            rec = dict(
                out=out_name,
                locustag=locustag or "",
                status=status,
                n_training_models=n_models if n_models is not None else "",
                busco_valid=mb.group(1) if mb else "",
                busco_validated=mo.group(1) if mo else "",
                pasa_fallback_models=mp.group(1) if mp else "",
                total_genes=ms.group(1) if ms else "",
                mtime=os.path.getmtime(log),
                last_line=last_meaningful_line(text),
                logpath=log,
            )
            key = locustag or out_name
            if key not in best or rec["mtime"] > best[key][0]:
                best[key] = (rec["mtime"], rec)

    records = [r for _, r in best.values()]

    # --- enrich with taxonomy + asm_stats ---
    for r in records:
        asmid, phylum, species, genus = by_locustag.get(r["locustag"], ("", "", "", ""))
        r["asmid"], r["phylum"], r["species"], r["genus"] = asmid, phylum, species, genus
        a = asm.get(asmid, {})
        r["total_bp"] = a.get("total_bp", "")
        r["contigs"] = a.get("contigs", "")
        r["n50"] = a.get("n50", "")
        risk, reason = asm_risk(a.get("total_bp"), a.get("contigs"), a.get("n50"), a.get("nonN_bp"))
        r["asm_risk"], r["asm_risk_reason"] = risk, reason

    os.makedirs(args.outdir, exist_ok=True)
    cols = ["asmid", "out", "locustag", "phylum", "genus", "species", "status",
            "n_training_models", "busco_valid", "busco_validated",
            "pasa_fallback_models", "total_genes", "total_bp", "contigs", "n50",
            "asm_risk", "asm_risk_reason", "last_line", "logpath"]

    def write_tsv(path, rows):
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for r in rows:
                w.writerow(r)

    all_path = os.path.join(args.outdir, "predict_outcomes.tsv")
    write_tsv(all_path, sorted(records, key=lambda r: (r["status"], r["phylum"], r["out"])))
    failed = [r for r in records if not r["status"].startswith("COMPLETED")]
    write_tsv(os.path.join(args.outdir, "failed_genomes.tsv"),
              sorted(failed, key=lambda r: (r["status"], r["out"])))

    # --- summary + crosstab ---
    lines = []
    lines.append(f"Logs scanned (unique): {n_logs}")
    lines.append(f"Genomes (deduped by locustag): {len(records)}")
    lines.append("")
    lines.append("Outcome counts (all phyla):")
    for st, c in sorted(Counter(r["status"] for r in records).items(), key=lambda x: -x[1]):
        lines.append(f"  {st:<28} {c}")
    lines.append("")

    sub = records
    if args.phylum:
        sub = [r for r in records if r["phylum"].lower().startswith(args.phylum.lower())]
        lines.append(f"--- restricted to phylum startswith '{args.phylum}': {len(sub)} genomes ---")
        lines.append("")

    # Crosstab: asm_risk (rows) x failed? (cols). Tests whether the pre-screen
    # actually catches real failures.
    lines.append("Crosstab  asm_risk x outcome  (does the asm_stats screen predict failures?)")
    lines.append(f"  {'asm_risk':<8} {'FAILED/INCOMPLETE':>18} {'COMPLETED':>10} {'total':>7}  failure_rate")
    risks = ["HIGH", "MED", "OK", "NA"]
    grid = {rk: [0, 0] for rk in risks}
    for r in sub:
        rk = r["asm_risk"]
        bad = not r["status"].startswith("COMPLETED")
        grid.setdefault(rk, [0, 0])[0 if bad else 1] += 1
    for rk in risks:
        bad, ok = grid.get(rk, [0, 0])
        tot = bad + ok
        rate = f"{100*bad/tot:.2f}%" if tot else "-"
        lines.append(f"  {rk:<8} {bad:>18} {ok:>10} {tot:>7}  {rate}")
    lines.append("")
    lines.append("Failed/incomplete genomes (this phylum subset):")
    lines.append(f"  {'out':<46}{'status':<26}{'models':>7}{'risk':>6}  asmid")
    for r in sorted(sub, key=lambda r: r["status"]):
        if r["status"].startswith("COMPLETED"):
            continue
        lines.append(f"  {r['out'][:45]:<46}{r['status']:<26}{str(r['n_training_models']):>7}{r['asm_risk']:>6}  {r['asmid']}")

    summary = "\n".join(lines)
    with open(os.path.join(args.outdir, "crosstab_summary.txt"), "w") as fh:
        fh.write(summary + "\n")
    print(summary)
    print(f"\nWrote:\n  {all_path}\n  {os.path.join(args.outdir,'failed_genomes.tsv')}\n  {os.path.join(args.outdir,'crosstab_summary.txt')}")


if __name__ == "__main__":
    main()
