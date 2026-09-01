#!/usr/bin/env python3
"""
Profile the paralogoscope pilot run from the Nextflow trace + storeDir artifacts.

For each pilot genome, extracts:
  - wgd dmd / ksd wall time (realtime, ms from trace)
  - peak RSS (rss, bytes from trace)
  - artifact file sizes (storeDir outputs under --outdir wgd_dmd/wgd_ksd subdirs)

Then fits a linear model (metrics ~ genes) and extrapolates to the full
4,365-genome dataset as "per genome per 5k genes" rates.

Usage:
  python3 profile_pilot.py [--trace LOGS/nextflow/paralogoscope_trace.xxx.txt]
                           [--counts analysis/WGD_PERFORMANCE_ANALYSIS/outputs/pilot_genome_counts.tsv]
                           [--allcounts /tmp/cds_counts.tsv]
                           [--outdir /bigdata/.../paralogoscope_pilot]
                           [--out analysis/WGD_PERFORMANCE_ANALYSIS/outputs/pilot_profile]
"""
import argparse, glob, os, re, sys

def parse_counts(path):
    m = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        assert "sampletag" in head and "genes" in head, f"bad counts header: {head}"
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            d = dict(zip(head, cols))
            m[d["sampletag"]] = int(d["genes"])
    return m

def parse_duration(s):
    """'11s' | '12m 34s' | '1h 1m 19s' -> seconds (float)."""
    if not s:
        return 0.0
    tot = 0.0
    for tok in str(s).strip().split():
        m = re.match(r"([0-9.]+)([smh])$", tok)
        if m:
            v = float(m.group(1))
            tot += {"s": 1, "m": 60, "h": 3600}[m.group(2)] * v
    return tot

def parse_memory(s):
    """'366.4 MB' | '2.1 GB' | raw bytes int -> bytes (float)."""
    if not s:
        return 0.0
    m = re.match(r"([0-9.]+)\s*(KB|MB|GB|TB)?$", str(s).strip())
    if not m:
        return float(s)
    v = float(m.group(1))
    unit = (m.group(2) or "B").upper()
    return v * {"B": 1, "KB": 1e3, "MB": 1e6, "GB": 1e9, "TB": 1e12}[unit]

def parse_trace(path):
    rows = []  # (name, sampletag, status, realtime_s, rss_bytes, exit)
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < len(head):
                continue
            d = dict(zip(head, cols))
            status = d.get("status", "")
            if status not in ("COMPLETED",):
                continue
            tag = d.get("tag", "").strip()
            rows.append({
                "name": d.get("name", ""),
                "tag": tag,
                "realtime": parse_duration(d.get("realtime", "")),
                "rss": parse_memory(d.get("rss", "")),
                "exit": d.get("exit", ""),
            })
    return rows

def classify(name):
    n = name.upper()
    if "WGD_DMD" in n:   return "dmd"
    if "WGD_KSD" in n:   return "ksd"
    if "WGD_SYN" in n:   return "syn"
    return "other"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--trace", default=None)
    ap.add_argument("--counts", default="analysis/WGD_PERFORMANCE_ANALYSIS/outputs/pilot_genome_counts.tsv")
    ap.add_argument("--allcounts", default="/tmp/cds_counts.tsv")
    ap.add_argument("--outdir", default="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/paralogoscope_pilot")
    ap.add_argument("--out", default="analysis/WGD_PERFORMANCE_ANALYSIS/outputs/pilot_profile")
    args = ap.parse_args()

    trace = args.trace or sorted(glob.glob("logs/nextflow/paralogoscope_trace.*.txt"))[-1]
    print(f"# trace : {trace}", file=sys.stderr)
    counts = parse_counts(args.counts)
    rows = parse_trace(trace)
    print(f"# rows  : {len(rows)} completed tasks", file=sys.stderr)

    # per-genome aggregation
    per = {}  # sampletag -> dict(stage -> [ms]), peak rss
    for r in rows:
        stage = classify(r["name"])
        tag = r["tag"]
        per.setdefault(tag, {"stage": {}, "rss": 0})
        per[tag]["stage"].setdefault(stage, []).append(r["realtime"])
        per[tag]["rss"] = max(per[tag]["rss"], r["rss"])

    # artifact sizes from storeDir
    def dirsize(pattern):
        tot = 0.0
        for p in glob.glob(pattern, recursive=True):
            if os.path.isfile(p):
                tot += os.path.getsize(p)
        return tot

    outdir = args.outdir.rstrip("/")
    sizes = {}
    for tag in counts:
        sizes[tag] = {
            "dmd": dirsize(f"{outdir}/wgd_dmd/*/{tag}.cds-transcripts.fa*"),
            "ksd": dirsize(f"{outdir}/wgd_ksd/*/{tag}.cds-transcripts.fa.tsv.ks.*"),
            "total": 0.0,
        }
        sizes[tag]["total"] = sizes[tag]["dmd"] + sizes[tag]["ksd"]

    # assemble measurement table
    recs = []
    for tag, g in counts.items():
        if tag not in per:
            continue
        st = per[tag]["stage"]
        colors = "rgba(30,80,180,0.3) 16"  # placeholder
        recs.append({
            "sampletag": tag, "genes": g,
            "dmd_min": sum(st.get("dmd", [])) / 60.0,
            "ksd_min": sum(st.get("ksd", [])) / 60.0,
            "tot_min": (sum(st.get("dmd", [])) + sum(st.get("ksd", []))) / 60.0,
            "peak_rss_gb": per[tag]["rss"] / 1e9,
            "dmd_mb": sizes[tag]["dmd"] / 1e6,
            "ksd_mb": sizes[tag]["ksd"] / 1e6,
            "art_mb": sizes[tag]["total"] / 1e6,
        })
    recs.sort(key=lambda r: r["genes"])

    import statistics as stt
    def linfit(xs, ys):
        n = len(xs)
        mx, my = stt.mean(xs), stt.mean(ys)
        sxy = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        sxx = sum((x - mx) ** 2 for x in xs)
        b = sxy / sxx if sxx else 0.0
        a = my - b * mx
        # R2
        sst = sum((y - my) ** 2 for y in ys)
        ssr = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
        r2 = 1 - ssr / sst if sst else float("nan")
        return a, b, r2

    genes = [r["genes"] for r in recs]
    tot_min = [r["tot_min"] for r in recs]
    rss_gb = [r["peak_rss_gb"] for r in recs]
    art_mb = [r["art_mb"] for r in recs]

    def model(metric, label, unit, per5k_unit):
        a, b, r2 = linfit(genes, metric)
        per5k_mean = stt.mean([v * 5000 / g for v, g in zip(metric, genes)]) if genes else float("nan")
        print(f"\n## {label}  (n={len(recs)} genomes)")
        print(f"  linear fit : y = {a:.3f} + {b:.6f} * genes   (R2={r2:.3f})")
        print(f"  per 5,000 genes (slope)   : {b*5000:.3f} {per5k_unit}")
        print(f"  per 5,000 genes (mean)    : {per5k_mean:.3f} {per5k_unit}")
        print(f"  median observed rate      : {stt.median([v*5000/g for v,g in zip(metric,genes)]):.3f} {per5k_unit}")
        return a, b, r2

    # full-dataset projection
    all_counts = {}
    if os.path.exists(args.allcounts):
        with open(args.allcounts) as fh:
            for line in fh:
                t, c = line.rstrip("\n").split("\t")
                try:
                    all_counts[t] = int(c)
                except ValueError:
                    pass
    total_genes = sum(all_counts.values())
    n_all = len(all_counts)
    med_all = stt.median(all_counts.values())
    print(f"\n## Full-dataset (from {args.allcounts})")
    print(f"  genomes={n_all}  total_genes={total_genes:,}  median_genes={med_all:,}")

    for name, metric, unit, p5u in [
        ("total runtime", tot_min, "min/genome", "min"),
        ("peak RSS", rss_gb, "GB/genome", "GB"),
        ("artifact size", art_mb, "MB/genome", "MB"),
    ]:
        a, b, r2 = model(metric, name, unit, p5u)
        per5k = b * 5000
        per_genome_med = (a + b * med_all) if genes else float("nan")
        if name == "total runtime":
            tot = per_genome_med * n_all / 60
            print(f"  => full run ({n_all} genomes @ median {med_all} genes): {tot:.0f} task-hours")
            for C in (16, 32, 56):
                print(f"     at concurrency {C}: {tot/C:.1f} h")
        else:
            print(f"  => full pipeline {name}: {per_genome_med:.1f} {unit}/genome @ median genes; "
                  f"total ≈ {per_genome_med*n_all:,.0f} {p5u}")

    # table
    out = args.out
    os.makedirs(os.path.dirname(out), exist_ok=True)
    with open(out + ".tsv", "w") as fh:
        fh.write("sampletag\tgenes\tdmd_min\tksd_min\ttot_min\tpeak_rss_gb\tdmd_mb\tksd_mb\tart_mb\n")
        for r in recs:
            fh.write("\t".join(str(r[k]) for k in
                     ["sampletag","genes","dmd_min","ksd_min","tot_min","peak_rss_gb","dmd_mb","ksd_mb","art_mb"]) + "\n")
    print(f"\nper-genome table -> {out}.tsv")
    print("sampletag\tgenes\ttot_min\trss_GB\tart_MB")
    for r in recs:
        print(f"{r['sampletag']}\t{r['genes']}\t{r['tot_min']:.1f}\t{r['peak_rss_gb']:.2f}\t{r['art_mb']:.1f}")

if __name__ == "__main__":
    main()
