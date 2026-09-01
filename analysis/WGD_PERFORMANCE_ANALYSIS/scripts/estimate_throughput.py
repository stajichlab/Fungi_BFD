#!/usr/bin/env python3
"""Extrapolate per-genome WGD timing to the full function_neurospora run.

Reads the measurements TSV (from collect_measurements.py) plus a few model
assumptions and emits wgd_throughput_estimate.tsv summarising per-step and
whole-set wall/CPU time at a range of SLURM concurrency levels.

Assumptions (documented in WGD_PERFORMANCE_ANALYSIS.md; tune as you like):
  - ksd wall @cpus_base scales ~linearly with multi-copy family count proxy
    (gene count) and ~1.7x when doubling threads (8 vs 4 cpus).
  - median genome gene count (from the real input set) vs the measured
    genome's gene count normalises the per-genome ksd time.
  - dmd is ~linear in genes; measured at ~1 min for a 12-13 k-gene genome.
  - syn is UNMEASURED at scale: treated as 0 here (opt-in); flagged in docs.

Usage:
    estimate_throughput.py MEASUREMENTS_TSV N_GENOMES --median-genes 8040 \
        --out OUT_TSV [--cpus [4 8]] [--thread-scaling 1.7]
"""

import argparse
import os
import sys

def parse_tsv(p):
    rows = []
    with open(p) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append(dict(zip(header, f)))
    return rows

def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("measurements_tsv")
    ap.add_argument("n_genomes", type=int)
    ap.add_argument("--median-genes", type=float, default=8040.0,
                    help="median gene count across the full input set (measured)")
    ap.add_argument("--out", default=None)
    ap.add_argument("--cpus", type=int, nargs="+", default=[8],
                    help="cpus per task to compute per-genome times for")
    ap.add_argument("--thread-scaling", type=float, default=1.7,
                    help="ksd speedup when doubling threads")
    ap.add_argument("--concurrency", type=int, nargs="+", default=[16, 24, 32, 40, 56])
    args = ap.parse_args()

    rows = parse_tsv(args.measurements_tsv)
    # completed tasks only, with begin+end
    done = [r for r in rows if r["exitcode"] == "0" and r["begin"]]
    import datetime
    ksd_done, dmd_done, syn_done = [], [], []
    for r in done:
        if r["step"] == "ksd" and r["end_log_mtime"]:
            b = datetime.datetime.strptime(r["begin"], "%Y-%m-%d %H:%M:%S")
            e = datetime.datetime.strptime(r["end_log_mtime"], "%Y-%m-%d %H:%M:%S")
            ksd_done.append((r, (e - b).total_seconds() / 60.0))
        elif r["step"] == "dmd" and r["end_log_mtime"]:
            b = datetime.datetime.strptime(r["begin"], "%Y-%m-%d %H:%M:%S")
            e = datetime.datetime.strptime(r["end_log_mtime"], "%Y-%m-%d %H:%M:%S")
            dmd_done.append((r, (e - b).total_seconds() / 60.0))

    def est(step_measure, scale, base_cpus, tgt_cpus, thread_scaling):
        if not step_measure:
            return None
        # normalise by measured genome genes -> median-genome genes
        genes = [num(r["n_genes"]) for r, _ in step_measure]
        g = [x for x in genes if x]
        gene_ratio = args.median_genes / (sum(g) / len(g)) if g else 1.0
        wall = [t for _, t in step_measure]
        base = sum(wall) / len(wall) * gene_ratio
        if tgt_cpus > base_cpus:
            base /= thread_scaling
        elif tgt_cpus < base_cpus:
            base *= thread_scaling
        return base

    print(f"measured completed tasks: dmd={len(dmd_done)} ksd={len(ksd_done)}")
    lines = []
    header = ["step", "cpus_per_task", "per_genome_min", "n_genomes",
              "total_task_h", "total_cpu_h", "concurrent", "wall_days"]
    for tgt in args.cpus:
        for step, meas, base_cpus in (("wgd_dmd", dmd_done, 4), ("wgd_ksd", ksd_done, 4)):
            per = est(meas, 1.0, base_cpus, tgt, args.thread_scaling)
            if per is None:
                per = 0.0
            total_min = per * args.n_genomes
            total_cpu_h = total_min * tgt / 60.0
            for c in args.concurrency:
                wall_days = total_min / 60.0 / c / 24.0
                lines.append([step, tgt, round(per, 1), args.n_genomes,
                              round(total_min / 60.0, 0), round(total_cpu_h, 0),
                              c, round(wall_days, 1)])
    out = args.out or os.path.join(os.path.dirname(args.measurements_tsv),
                                   "wgd_throughput_estimate.tsv")
    with open(out, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for l in lines:
            fh.write("\t".join(str(x) for x in l) + "\n")
    print(f"wrote {out}")
    # quick printed summary for dmd/ksd at 8 cpus
    for step in ("wgd_dmd", "wgd_ksd"):
        rows_s = [l for l in lines if l[0] == step and l[1] == args.cpus[0]]
        if rows_s:
            p = rows_s[0]
            print(f"{step} @ {p[1]}cpus: {p[2]} min/genome -> {p[3]} genomes = "
                  f"{p[4]} task-h ({p[5]} cpu-h); @{p[6]} concurrent = {p[7]} days")

if __name__ == "__main__":
    main()
