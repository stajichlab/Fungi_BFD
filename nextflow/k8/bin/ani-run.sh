#!/usr/bin/env bash
# k8/bin/ani-run.sh — launch one compute-phase (run_ani_compute.nf) run on
# Nautilus, detached, from your laptop.
#
# Only does the SKANI/mash/sourmash/fastani sketch+compare — no REPORT_ANI or
# COMBINE_ANI_TABLE (run k8/bin/ani-gather.sh for those, once you're ready to
# aggregate). Launches with `nohup ... &` inside the head pod so it survives
# your kubectl exec session ending; logs land durably on the PVC.
#
# Usage:
#   ani-run.sh --taxon GENUS:Yarrowia --compare SPECIES [--name yarrowia] \
#              [--method skani] [-- --n_test 2 --skani_preset fast]
#
# Anything after a bare `--` is passed straight through as extra nextflow
# params (e.g. -- --n_test 2), letting you override any default without
# editing this script.
set -euo pipefail

NAMESPACE=ucr-stajichlab
POD=bfd-nextflow-head

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

RUN_DIR="/root/runs/${NAME}"
LOG_DIR="/workspace/logs/cli-runs"
LOG_FILE="${LOG_DIR}/${NAME}.log"
PARAMS_FILE="/workspace/runs/${NAME}/params.yaml"

echo "==> [${NAME}] taxon=${TAXON} compare=${COMPARE} method=${METHOD}"
echo "==> [${NAME}] writing params to ${PARAMS_FILE}"

kubectl exec -i -n "$NAMESPACE" "$POD" -- sh -c "mkdir -p \$(dirname '${PARAMS_FILE}') && cat > '${PARAMS_FILE}'" <<EOF
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

echo "==> [${NAME}] launching detached; log: ${LOG_FILE}"
kubectl exec -n "$NAMESPACE" "$POD" -- sh -c "
  mkdir -p '${RUN_DIR}' '${LOG_DIR}'
  cp /workspace/repo/samples.csv '${RUN_DIR}/samples.csv'
  cd '${RUN_DIR}'
  nohup nextflow run /workspace/repo/nextflow/run_ani_compute.nf \
    -c /workspace/repo/nextflow/nextflow.config \
    -profile compare_ani_k8s \
    -params-file '${PARAMS_FILE}' \
    -resume ${EXTRA_ARGS} \
    > '${LOG_FILE}' 2>&1 < /dev/null &
  disown
  sleep 1
"
echo "==> [${NAME}] started. Check with: ani-status.sh ${NAME}"
