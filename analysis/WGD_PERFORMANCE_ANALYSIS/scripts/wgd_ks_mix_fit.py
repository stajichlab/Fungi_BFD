#!/usr/bin/env python3
"""Fit wgd-style Ks mixture models to one genome and emit peak summaries.

Runs INSIDE the wgd container (hyphaltip_wgd2_complete SIF) so it reuses the
wgd package's own mixture-model code (wgd.mix.filter_group_data +
wgd.mix.fit_gmm / fit_bgmm -- the exact machinery behind `wgd mix`). The wgd
CLI only writes per-pair posterior-probability TSVs and plots; the component
means/weights ("mean Ks of the peaks") are not exported anywhere, so this
worker captures them and emits a JSON that the host driver
(build_wgd_ksd_summary.py) folds into wgd_ksd_summary.parquet.

Input is a minimal per-genome Ks TSV produced by the host driver from
tables/wgd.ks.parquet: column 1 = pair id, then columns `alignmentlength`
and `dS` (the only fields wgd.mix actually reads). Output is a single JSON on
stdout or to --out.

Semantics mirror `wgd mix --method gmm` exactly:
  - filter: alignmentlength >= --filters, min_ks < dS <= max_ks,
  - X = log(dS) of the surviving pairs (raw pair rows; wgd mixes every pair,
    it does not node-average here),
  - GMM (full covariance) fitted for every k in [min,max] components,
    n_init restarts, max_iter EM iterations,
  - best model = min BIC (GMM); BGMM keeps the max-component model (wgd
    behavior -- weights are then informative but model count is not chosen).

Also emits (optionally) a KDE of the filtered log-Ks on a linear Ks grid,
the planned wgd_ksd_density.parquet artifact this framework feeds.

Usage:
    wgd_ks_mix_fit.py KSIN --out OUT.json [--method gmm|bgmm]
        [--min-components 1] [--max-components 4] [--n-init 200]
        [--max-iter 200] [--filters 300] [--ks-min 0] [--ks-max 5]
        [--density-grid 200] [--density-max 5]
"""

import argparse
import json
import sys

import numpy as np
import pandas as pd


def components_sorted_by_mean(model):
    """Return [(mean_ks, sd, weight)] ordered by increasing mean Ks."""
    comps = []
    for j in range(len(model.means_)):
        mean_ks = float(np.exp(model.means_[j][0]))
        sd = float(np.sqrt(model.covariances_[j][0][0]))
        weight = float(model.weights_[j])
        comps.append((mean_ks, sd, weight))
    comps.sort(key=lambda c: c[0])
    return [{"mean_ks": c[0], "sd": c[1], "weight": c[2]} for c in comps]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("ksin", help="per-genome Ks TSV (index=pair, alignmentlength, dS)")
    ap.add_argument("--out", help="output JSON (default: stdout)")
    ap.add_argument("--method", choices=["gmm", "bgmm"], default="gmm")
    ap.add_argument("--min-components", type=int, default=1)
    ap.add_argument("--max-components", type=int, default=4)
    ap.add_argument("--n-init", type=int, default=200)
    ap.add_argument("--max-iter", type=int, default=200)
    ap.add_argument("--seed", type=int, default=2352890,
                    help="random_state for the GMM (wgd peak's default seed)")
    ap.add_argument("--filters", type=int, default=300)
    ap.add_argument("--ks-min", type=float, default=0.0)
    ap.add_argument("--ks-max", type=float, default=5.0)
    ap.add_argument("--density-grid", type=int, default=200,
                    help="KDE grid points on [0, density-max] (0 disables)")
    ap.add_argument("--density-max", type=float, default=5.0)
    args = ap.parse_args()

    from wgd.mix import filter_group_data, get_array_for_mixture, fit_gmm, fit_bgmm

    df = pd.read_csv(args.ksin, index_col=0, sep="\t")
    n_rows_input = int(df.shape[0])
    df = filter_group_data(df, args.filters, args.ks_min, args.ks_max)
    X = get_array_for_mixture(df)
    n_ks = int(X.shape[0])

    result = {
        "n_rows_input": n_rows_input,
        "n_ks": n_ks,
        "method": args.method,
        "filters": args.filters,
        "ks_range": [args.ks_min, args.ks_max],
        "n_init": args.n_init,
        "max_iter": args.max_iter,
    }

    if n_ks == 0:
        result["models"] = []
        result["best"] = None
        if args.density_grid > 0:
            result["density"] = {"grid": [], "y": []}
        payload = json.dumps(result)
        open(args.out, "w").write(payload) if args.out else sys.stdout.write(payload)
        return

    # Only fit up to min(max_components, distinct usable values); a full
    # covariance model needs at least 1 observation per component but the
    # range is already bounded by the caller (1..4).
    models = []
    best = None
    if args.method == "gmm":
        fitted, bic, aic, best = fit_gmm(
            X, args.min_components, args.max_components,
            max_iter=args.max_iter, n_init=args.n_init,
            random_state=args.seed)
        models = []
        for k, m in enumerate(fitted):
            nc = args.min_components + k
            models.append({
                "n_components": nc,
                "bic": float(bic[k]),
                "aic": float(aic[k]),
                "components": components_sorted_by_mean(m),
            })
    else:
        fitted = fit_bgmm(
            X, args.min_components, args.max_components, gamma=0.001,
            max_iter=args.max_iter, n_init=args.n_init)
        for k, m in enumerate(fitted):
            nc = args.min_components + k
            models.append({
                "n_components": nc,
                "bic": None,
                "aic": None,
                "components": components_sorted_by_mean(m),
            })
        best = fitted[-1]

    best_rec = None
    if best is not None:
        best_rec = {
            "n_components": int(best.n_components),
            "bic": float(best.bic(X)) if args.method == "gmm" else None,
            "aic": float(best.aic(X)) if args.method == "gmm" else None,
            "components": components_sorted_by_mean(best),
        }
    result["models"] = models
    result["best"] = best_rec

    if args.density_grid > 0:
        from scipy.stats import gaussian_kde
        grid = np.linspace(0.0, args.density_max, args.density_grid)
        if X.shape[0] > 1:
            kde = gaussian_kde(X[:, 0], bw_method="silverman")
            y = np.exp(kde(np.log(np.maximum(grid, 1e-9)))).tolist()
        else:
            y = [np.nan] * args.density_grid
        result["density"] = {"grid": grid.tolist(), "y": y}

    payload = json.dumps(result)
    if args.out:
        with open(args.out, "w") as fh:
            fh.write(payload)
    else:
        sys.stdout.write(payload)


if __name__ == "__main__":
    main()
