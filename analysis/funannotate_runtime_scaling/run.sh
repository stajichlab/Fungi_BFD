#!/bin/bash
# Reproduce the funannotate runtime-scaling analysis (runtime vs. genome size /
# gene count / repeat content / scaffold count, faceted by pipeline process).
#
# Scripts live in Fungi_BFD_runs (the pipeline repo), not here — this analysis
# is a benchmarking view over that repo's trace files, run against the shared
# results/genome_stats_by_name DB, so it's kept in Fungi_BFD (the mycelium
# analysis home) rather than duplicated per pipeline root.
set -euo pipefail
RUNS_ROOT=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
cd "$(dirname "$0")" || exit 1

/usr/bin/python3.12 "$RUNS_ROOT/do_annotation_asco/scripts/build_runtime_metadata_table.py" \
  --roots "$RUNS_ROOT/do_annotation" "$RUNS_ROOT/do_annotation_asco" "$RUNS_ROOT/do_annotation_Rhod" \
  --outdir "$RUNS_ROOT/analysis" \
  --also-csv --no-merge

/usr/bin/python3.12 "$RUNS_ROOT/do_annotation_asco/scripts/plot_runtime_relationships.py" \
  --table "$RUNS_ROOT/analysis/funannotate_runtime_metadata.parquet" \
  --outdir "$RUNS_ROOT/analysis"

/usr/bin/python3.12 - "$RUNS_ROOT/analysis/funannotate_runtime_metadata.parquet" outputs/fit_summary.csv <<'PYEOF'
import sys
import numpy as np
import pandas as pd

table_path, out_path = sys.argv[1], sys.argv[2]
df = pd.read_parquet(table_path)

metrics = [('genome_size_bp', True), ('gene_count', True),
           ('repeat_masked_pct', False), ('num_scaffolds', True)]
procs = ['FUNANNOTATE_TRAIN', 'FUNANNOTATE_PREDICT', 'GENEMARK_RUN',
         'FUNANNOTATE_PREDICT_SIB', 'GENEMARK_RUN_SIB']

rows = []
for col, logx in metrics:
    for proc in procs:
        sub = df[(df.process == proc) & (df.status == 'COMPLETED')
                 & (~df['likely_skip_exit']) & df[col].notna() & df.realtime_s.notna()]
        sub = sub[(sub[col] > 0) & (sub.realtime_s > 0)]
        if len(sub) < 10:
            continue
        x = np.log10(sub[col].values) if logx else sub[col].values.astype(float)
        y = np.log10(sub.realtime_s.values / 60.0)
        coeffs = np.polyfit(x, y, 1)
        pred = np.polyval(coeffs, x)
        ss_res = np.sum((y - pred) ** 2)
        ss_tot = np.sum((y - y.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
        rows.append(dict(process=proc, metric=col, log_x=logx, n_completed=len(sub),
                          r_squared=round(r2, 4), slope=round(coeffs[0], 4),
                          intercept=round(coeffs[1], 4)))
pd.DataFrame(rows).to_csv(out_path, index=False)
PYEOF

cp "$RUNS_ROOT"/analysis/runtime_vs_*.png "$RUNS_ROOT"/analysis/runtime_vs_*.pdf \
   "$RUNS_ROOT"/analysis/runtime_outliers.csv \
   outputs/
