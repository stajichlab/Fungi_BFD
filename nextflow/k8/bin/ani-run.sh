#!/usr/bin/env bash
# k8/bin/ani-run.sh — launch one compute-phase (run_ani_compute.nf) run on
# Nautilus, detached, from your laptop.
#
# Only does the SKANI/mash/sourmash/fastani sketch+compare — no REPORT_ANI or
# COMBINE_ANI_TABLE (run k8/bin/ani-gather.sh for those, once you're ready to
# aggregate). Runs as a Kubernetes Job (k8/bin/nf-job.sh) whose command is
# `nextflow run`, so nothing stays running after the pipeline ends (NRP
# prohibits idle head pods). Logs: `kubectl logs job/nf-ani-<name>` and
# /workspace/logs/cli-runs/ani-<name>.log on the PVC.
#
# Usage:
#   ani-run.sh --taxon GENUS:Yarrowia --compare SPECIES [--name yarrowia] \
#              [--method skani] [-- --n_test 2 --skani_preset fast]
#
# Anything after a bare `--` is passed straight through as extra nextflow
# params (e.g. -- --n_test 2), letting you override any default without
# editing this script.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

TAXON=""
COMPARE="GENUS"
METHOD="skani"
NAME=""
EXTRA_ARGS=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxon) TAXON="$2"; shift 2 ;;
    --compare) COMPARE="$2"; shift 2 ;;
    --method) METHOD="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    --) shift; EXTRA_ARGS="$*"; break ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$TAXON" ]]; then
  echo "Usage: $0 --taxon RANK:VALUE --compare LEVEL [--name slug] [--method skani|mash|sourmash|fastani] [-- extra nextflow args]" >&2
  exit 1
fi

if [[ -z "$NAME" ]]; then
  NAME=$(echo "$TAXON" | sed -E 's/.*://' | tr '[:upper:]' '[:lower:]' | tr -c 'a-z0-9' '-' | sed -E 's/-+$//')
fi

RUN_DIR="/workspace/runs/ani-${NAME}"
PARAMS_FILE="$(mktemp -t "ani-${NAME}.XXXX").yaml"
trap 'rm -f "$PARAMS_FILE"' EXIT

echo "==> [${NAME}] taxon=${TAXON} compare=${COMPARE} method=${METHOD}"

cat > "$PARAMS_FILE" <<EOF
samples:               "samples.csv"
genome_name_style:     "asmid"
genome_dir:            "s3://stajichlab/BFD/input_clean_genomes"
genome_suffix:         ".masked.fasta.gz"
compare:               "${COMPARE}"
taxon:                 "${TAXON}"
ani_method:            "${METHOD}"
sketch_cache:          "work/ANI/sketch_cache"
ani_cluster_threshold: 95.0
ani_outlier_threshold: 90.0
min_group_size:        2
outdir:                "s3://stajichlab/BFD/results/ANI"
skani_preset:          "medium"
skani_min_af:          15
skani_compression:     0
skani_sketch_chunk:    50
n_test:                0
EOF

# shellcheck disable=SC2086  # EXTRA_ARGS is deliberately word-split into nextflow args
"${HERE}/nf-job.sh" --name "ani-${NAME}" --run-dir "$RUN_DIR" --params "$PARAMS_FILE" -- \
  /workspace/repo/nextflow/run_ani_compute.nf \
  -c /workspace/repo/nextflow/nextflow.config \
  -profile compare_ani_k8s ${EXTRA_ARGS}
echo "==> [${NAME}] Check with: ani-status.sh ani-${NAME}"
