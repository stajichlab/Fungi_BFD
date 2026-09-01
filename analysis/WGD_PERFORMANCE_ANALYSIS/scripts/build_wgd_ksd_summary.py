#!/usr/bin/env python3
"""Build wgd_ksd_summary.parquet + wgd_ksd_density.parquet from the merged
wgd Ks tables, using wgd's own mixture models for per-genome Ks peaks.

Answers the per-genome question "how many duplicates, and what are the mean
Ks of the peaks?" from the merge surface alone:

  * number of duplicates  -> already in tables/wgd.ks.summary.parquet
    (n_pairs, n_pairs_with_ds, n_families); carried through unchanged.
  * mean Ks of the peaks  -> per genome, fit GMM/BGMM to the filtered Ks
    distribution exactly as `wgd mix` does (wgd.mix.filter_group_data +
    wgd.mix.fit_gmm, min-BIC model selection, components ordered by mean).
    The wgd CLI never exports the component means, so each genome's fit is
    done by the container worker wgd_ks_mix_fit.py inside the wgd SIF; the
    best model's components become the peaks (mean Ks, sd, weight), and the
    n_ks_peaks column is the selected component count.

Outputs (zstd parquet, written to --outdir):
  * wgd_ksd_summary.parquet -- one row per genome: counts, fit meta, and
    wide per-peak columns peak1..peak{max-components} (mean_ks, sd, weight),
    ordered by increasing mean Ks; fit_ok marks worker failures.
  * wgd_ksd_density.parquet -- long form genome x shared Ks grid of the
    per-genome KDE (silverman), the artifact the aggregate figures consume.

Usage:
    build_wgd_ksd_summary.py --sif PATH/TO/wgd.sif [--ks-parquet tables/wgd.ks.parquet]
        [--summary-parquet tables/wgd.ks.summary.parquet]
        [--outdir analysis/WGD_PERFORMANCE_ANALYSIS/tables]
        [--jobs 8] [--tmpdir DIR] [--keep-tmp]
        [--method gmm|bgmm] [--min-components 1] [--max-components 4]
        [--n-init 200] [--max-iter 200] [--seed 2352890]
        [--filters 300] [--ks-min 0] [--ks-max 5]
        [--density-grid 200] [--density-max 5]
"""

import argparse
import concurrent.futures as cf
import glob
import json
import os
import shutil
import subprocess
import sys
import tempfile

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def find_apptainer():
    """Locate the apptainer binary (PATH, or the HPCC module install)."""
    exe = shutil.which("apptainer")
    if exe:
        return exe
    hits = sorted(glob.glob(
        "/opt/linux/rocky/8.x/x86_64/pkgs/apptainer/*/bin/apptainer"))
    if hits:
        return hits[-1]
    return None


def write_pseudo_ks(ks_df, genome, path):
    """Minimal per-genome Ks TSV the wgd worker reads (index=pair).

    family + node are required: wgd.mix.filter_group_data node-averages the
    pairs with df.groupby(['family', 'node']).mean() before fitting.
    """
    sub = ks_df.loc[ks_df["genome"] == genome,
                    ["pair", "family", "node",
                     "alignmentlength", "dS"]]
    sub.to_csv(path, sep="\t", index=False)


def run_worker(sif, worker, tsv, out_json, fit_args, bind_dirs):
    binds = "".join(" --bind {}".format(d) for d in bind_dirs)
    cmd = [find_apptainer(), "exec"] + binds.split() + [sif, "python3", worker,
           tsv, "--out", out_json] + fit_args
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                          text=True)
    return proc


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sif", required=True, help="path to the wgd singularity image")
    ap.add_argument("--ks-parquet", default="tables/wgd.ks.parquet")
    ap.add_argument("--summary-parquet", default="tables/wgd.ks.summary.parquet")
    ap.add_argument("--outdir", default=os.path.join(HERE, "..", "tables"))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--tmpdir", default=None,
                    help="scratch for per-genome TSVs/JSONs (default: mkdtemp)")
    ap.add_argument("--keep-tmp", action="store_true")
    ap.add_argument("--worker", default=os.path.join(HERE, "wgd_ks_mix_fit.py"))
    ap.add_argument("--method", choices=["gmm", "bgmm"], default="gmm")
    ap.add_argument("--min-components", type=int, default=1)
    ap.add_argument("--max-components", type=int, default=4)
    ap.add_argument("--n-init", type=int, default=200)
    ap.add_argument("--max-iter", type=int, default=200)
    ap.add_argument("--seed", type=int, default=2352890)
    ap.add_argument("--filters", type=int, default=300)
    ap.add_argument("--ks-min", type=float, default=0.0)
    ap.add_argument("--ks-max", type=float, default=5.0)
    ap.add_argument("--density-grid", type=int, default=200)
    ap.add_argument("--density-max", type=float, default=5.0)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    ks = pd.read_parquet(args.ks_parquet,
                         columns=["genome", "species_prefix", "pair", "family",
                                  "node", "alignmentlength", "dS"])
    genomes = sorted(ks["genome"].unique())
    print("{} genomes in {}".format(len(genomes), args.ks_parquet),
          file=sys.stderr)

    summary = pd.read_parquet(args.summary_parquet)
    counts = summary.set_index("genome").reindex(genomes).reset_index()

    tmp = args.tmpdir or tempfile.mkdtemp(prefix="wgd_ksd_")
    if not os.path.isdir(tmp):
        os.makedirs(tmp, exist_ok=True)
    print("tmpdir: {}".format(tmp), file=sys.stderr)

    fit_args = ["--method", args.method,
                "--min-components", str(args.min_components),
                "--max-components", str(args.max_components),
                "--n-init", str(args.n_init),
                "--max-iter", str(args.max_iter),
                "--seed", str(args.seed),
                "--filters", str(args.filters),
                "--ks-min", str(args.ks_min),
                "--ks-max", str(args.ks_max),
                "--density-grid", str(args.density_grid),
                "--density-max", str(args.density_max)]
    bind_dirs = sorted({os.path.dirname(os.path.abspath(args.worker)),
                        os.path.abspath(tmp)})

    def fit_one(genome):
        tsv = os.path.join(tmp, "{}.ks.tmp.tsv".format(genome))
        out_json = os.path.join(tmp, "{}.fit.json".format(genome))
        write_pseudo_ks(ks, genome, tsv)
        proc = run_worker(args.sif, os.path.abspath(args.worker),
                          tsv, out_json, fit_args, bind_dirs)
        if proc.returncode != 0:
            err = (proc.stderr or proc.stdout)[-2000:]
            return genome, None, err
        with open(out_json) as fh:
            return genome, json.load(fh), None

    results = {}
    with cf.ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futs = {ex.submit(fit_one, g): g for g in genomes}
        done = 0
        for fut in cf.as_completed(futs):
            g, fit, err = fut.result()
            results[g] = (fit, err)
            done += 1
            if done % 100 == 0 or done == len(genomes):
                print("  fitted {}/{}".format(done, len(genomes)),
                      file=sys.stderr)

    n_peaks_max = args.max_components
    peak_cols = []
    for i in range(1, n_peaks_max + 1):
        peak_cols += [
            "peak{}_mean_ks".format(i),
            "peak{}_sd".format(i),
            "peak{}_weight".format(i),
        ]

    rows = []
    density_rows = []
    for g in genomes:
        fit, err = results[g]
        row = {"genome": g}
        if fit is None:
            row.update({"fit_ok": False, "fit_error": err})
            for c in peak_cols:
                row[c] = None
            rows.append(row)
            continue
        best = fit["best"]
        row.update({
            "fit_ok": True,
            "method": fit["method"],
            "n_ks": fit["n_ks"],
            "n_ks_peaks": best["n_components"] if best else 0,
            "bic_best": best["bic"] if best else None,
            "aic_best": best["aic"] if best else None,
        })
        comps = (best or {}).get("components", []) or []
        for i, c in enumerate(comps, start=1):
            if i > n_peaks_max:
                break
            row["peak{}_mean_ks".format(i)] = c["mean_ks"]
            row["peak{}_sd".format(i)] = c["sd"]
            row["peak{}_weight".format(i)] = c["weight"]
        for c in peak_cols:
            row.setdefault(c, None)
        rows.append(row)

        if fit.get("density") and fit["density"]["grid"]:
            grid = fit["density"]["grid"]
            y = fit["density"]["y"]
            for gx, gy in zip(grid, y):
                density_rows.append({"genome": g, "ks": gx, "density": gy})

    out_sum = pd.DataFrame(rows)
    out_sum = out_sum.merge(
        counts[["genome", "species_prefix", "n_pairs", "n_pairs_with_ds",
                "n_families"]], on="genome", how="left")
    cols = ["genome", "species_prefix", "n_pairs", "n_pairs_with_ds",
            "n_families", "n_ks", "n_ks_peaks", "bic_best", "aic_best",
            "fit_ok", "method"] + peak_cols + ["fit_error"]
    out_sum = out_sum[[c for c in cols if c in out_sum.columns]]
    sum_path = os.path.join(args.outdir, "wgd_ksd_summary.parquet")
    out_sum.to_parquet(sum_path, index=False, compression="zstd")
    print("wrote {} ({} genomes, {} cols)".format(
        sum_path, len(out_sum), len(out_sum.columns)), file=sys.stderr)

    if density_rows:
        out_den = pd.DataFrame(density_rows)
        den_path = os.path.join(args.outdir, "wgd_ksd_density.parquet")
        out_den.to_parquet(den_path, index=False, compression="zstd")
        print("wrote {} ({} rows)".format(den_path, len(out_den)),
              file=sys.stderr)
    else:
        print("no density rows (density-grid=0 or no fits)", file=sys.stderr)

    if not args.keep_tmp:
        import shutil
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
