#!/usr/bin/env python3
"""
plot_runtime_relationships.py — runtime vs. genome-size / gene-count /
repeat-content / scaffold-count, faceted by pipeline process, with outlier
flagging, from the table built by build_runtime_metadata_table.py.

Reads: analysis/funannotate_runtime_metadata.parquet
Writes (into --outdir):
  runtime_vs_genome_size.png
  runtime_vs_gene_count.png
  runtime_vs_repeat_pct.png
  runtime_vs_num_scaffolds.png
  runtime_outliers.csv   — per-facet log-log (or linear, for repeat %) fit
                           residual outliers, for follow-up investigation.

Design notes:
  - Facets are small multiples, one subplot per process, sharing the y axis
    (runtime in minutes, log scale) so magnitudes are comparable across
    processes at a glance.
  - Color encodes task status (COMPLETED / FAILED / other) — a fixed,
    non-cycled mapping — so failed/aborted attempts stay visible rather than
    being silently dropped; this table intentionally includes aborted runs
    as part of a broader funannotate benchmarking record.
  - An OLS fit (log-log for size/gene-count/scaffolds, linear for repeat %)
    is drawn per facet over COMPLETED points only, and the top-|residual|
    points per facet are labeled directly on the plot and written to
    runtime_outliers.csv for follow-up.
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STATUS_COLORS = {
    'COMPLETED': '#4c72b0',
    'FAILED': '#c44e52',
    'CACHED': '#8c8c8c',
    'SKIP_EXIT': '#55a868',
}
DEFAULT_COLOR = '#dd8452'  # any other status (ABORTED, killed, etc.)

FACET_PROCESSES = ['FUNANNOTATE_TRAIN', 'FUNANNOTATE_PREDICT', 'GENEMARK_RUN',
                    'FUNANNOTATE_PREDICT_SIB', 'GENEMARK_RUN_SIB']

METRICS = [
    # (column, output name, x label, log-x, min n per facet to fit/plot)
    ('genome_size_bp', 'runtime_vs_genome_size', 'Genome size (bp)', True, 5),
    ('gene_count', 'runtime_vs_gene_count', 'Predicted gene count', True, 5),
    ('repeat_masked_pct', 'runtime_vs_repeat_pct', 'Repeat-masked (%)', False, 5),
    ('num_scaffolds', 'runtime_vs_num_scaffolds', 'Scaffold count', True, 5),
]

N_LABEL_OUTLIERS = 5  # labeled directly on each facet


def status_color(status):
    return STATUS_COLORS.get(status, DEFAULT_COLOR)


def fit_and_residuals(x, y, log_x):
    """OLS fit of y (already log10 minutes) vs x (log10'd here if log_x).
    Returns (slope, intercept, x_for_fit, predicted, residuals)."""
    xv = np.log10(x) if log_x else x.astype(float)
    yv = y
    coeffs = np.polyfit(xv, yv, 1)
    predicted = np.polyval(coeffs, xv)
    residuals = yv - predicted
    return coeffs, xv, predicted, residuals


def plot_metric(df, column, out_name, x_label, log_x, min_n, outdir, outlier_rows,
                 success_only=False, collect_outliers=True):
    """success_only=True drops FAILED/CACHED/SKIP_EXIT rows before plotting entirely
    (not just from the fit) — for "how long does a real run take" figures. The OLS fit
    and flagged outliers are unaffected either way (both already fit COMPLETED,
    non-skip-exit rows only), so pass collect_outliers=False on this pass to avoid
    duplicating outlier_rows when both variants are generated from the same df."""
    facets = [p for p in FACET_PROCESSES if p in df['process'].unique()]
    if not facets:
        return
    ncols = 3
    nrows = int(np.ceil(len(facets) / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.6 * nrows), squeeze=False)
    axes_flat = axes.flatten()

    for i, proc in enumerate(facets):
        ax = axes_flat[i]
        sub = df[(df['process'] == proc) & df[column].notna() & df['realtime_s'].notna()].copy()
        if success_only:
            sub = sub[(sub['status'] == 'COMPLETED') & (~sub.get('likely_skip_exit', False))]
        sub = sub[sub[column] > 0]
        sub = sub[sub['realtime_s'] > 0]
        if sub.empty:
            ax.set_visible(False)
            continue

        sub['runtime_min'] = sub['realtime_s'] / 60.0
        sub['log_runtime_min'] = np.log10(sub['runtime_min'])
        # SKIP_EXIT: funannotate's own existing-output check short-circuited the task
        # (confirmed 2026-09-11) -- a benign no-op, not a real execution, so it gets its
        # own display bucket rather than being counted as a genuine COMPLETED runtime.
        sub['display_status'] = np.where(sub.get('likely_skip_exit', False),
                                          'SKIP_EXIT', sub['status'])

        for status, grp in sub.groupby('display_status'):
            ax.scatter(grp[column], grp['runtime_min'], s=14, alpha=0.5,
                       color=status_color(status), label=status, edgecolors='none')

        completed = sub[(sub['status'] == 'COMPLETED') & (~sub.get('likely_skip_exit', False))]
        if len(completed) >= min_n:
            coeffs, xv, predicted, residuals = fit_and_residuals(
                completed[column].values, completed['log_runtime_min'].values, log_x)
            order = np.argsort(xv)
            x_sorted = completed[column].values[order]
            fit_line_y = 10 ** predicted[order]
            ax.plot(x_sorted, fit_line_y, color='black', lw=1.2, ls='--',
                    label='OLS fit (COMPLETED)')

            completed = completed.assign(residual=residuals, fit_xv=xv)
            top_outliers = completed.reindex(
                completed['residual'].abs().sort_values(ascending=False).index
            ).head(N_LABEL_OUTLIERS)
            legend_lines = []
            for j, (_, row) in enumerate(top_outliers.iterrows(), start=1):
                ax.scatter([row[column]], [row['runtime_min']], s=28, facecolors='none',
                           edgecolors='black', linewidths=1.0, zorder=5)
                ax.annotate(str(j), (row[column], row['runtime_min']),
                            fontsize=7, fontweight='bold', color='black',
                            xytext=(4, 4), textcoords='offset points', zorder=6)
                short = row['species'] if len(row['species']) <= 28 else row['species'][:26] + '…'
                legend_lines.append(f"{j}. {short}")
                if collect_outliers:
                    outlier_rows.append({
                        'process': proc,
                        'metric': column,
                        'root': row['root'],
                        'species': row['species'],
                        'hash': row['hash'],
                        'status': row['status'],
                        'metric_value': row[column],
                        'runtime_minutes': row['runtime_min'],
                        'predicted_runtime_minutes': 10 ** (row['log_runtime_min'] - row['residual']),
                        'log10_residual': row['residual'],
                        'workdir': row['workdir'],
                    })
            if legend_lines:
                ax.text(0.98, 0.02, '\n'.join(legend_lines), transform=ax.transAxes,
                        fontsize=6, ha='right', va='bottom', family='monospace',
                        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                  edgecolor='0.7', alpha=0.85), zorder=7)

        ax.set_title(f'{proc}  (n={len(sub)})', fontsize=10)
        ax.set_xlabel(x_label, fontsize=9)
        if log_x:
            ax.set_xscale('log')
        ax.set_yscale('log')
        ax.tick_params(labelsize=8)
        ax.legend(fontsize=7, loc='upper left')

    for j in range(len(facets), len(axes_flat)):
        axes_flat[j].set_visible(False)

    fig.supylabel('Runtime (minutes, log scale)', fontsize=10)
    title_suffix = ' — successful runs only (FAILED/CACHED/SKIP_EXIT omitted)' if success_only else ''
    fig.suptitle(f'Runtime vs. {x_label}{title_suffix}', fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    png_path = outdir / f'{out_name}.png'
    pdf_path = outdir / f'{out_name}.pdf'
    fig.savefig(png_path, dpi=150)
    fig.savefig(pdf_path)
    plt.close(fig)
    print(f"Wrote {png_path} and {pdf_path}", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--table',
                     default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis/funannotate_runtime_metadata.parquet')
    ap.add_argument('--outdir',
                     default='/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/analysis')
    args = ap.parse_args()

    df = pd.read_parquet(args.table)
    print(f"Loaded {len(df)} rows from {args.table}", file=sys.stderr)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    outlier_rows = []
    for column, out_name, x_label, log_x, min_n in METRICS:
        plot_metric(df, column, out_name, x_label, log_x, min_n, outdir, outlier_rows)
        plot_metric(df, column, f'{out_name}_success_only', x_label, log_x, min_n, outdir,
                    outlier_rows, success_only=True, collect_outliers=False)

    if outlier_rows:
        out_df = pd.DataFrame(outlier_rows)
        out_df = out_df.sort_values('log10_residual', key=lambda s: s.abs(), ascending=False)
        out_path = outdir / 'runtime_outliers.csv'
        out_df.to_csv(out_path, index=False)
        print(f"Wrote {len(out_df)} outlier rows -> {out_path}", file=sys.stderr)
    else:
        print("No outliers flagged (insufficient data per facet).", file=sys.stderr)

    return 0


if __name__ == '__main__':
    sys.exit(main())
