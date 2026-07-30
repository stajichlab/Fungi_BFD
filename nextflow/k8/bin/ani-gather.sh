#!/usr/bin/env bash
# k8/bin/ani-gather.sh — run REPORT_ANI + COMBINE_ANI_TABLE for a taxon whose
# compute phase (ani-run.sh) has already published *.ani.tsv files to S3.
#
# Safe to run repeatedly / before all groups finish: run_ani_gather.nf only
# picks up groups that have actually published so far, skipping the rest.
# Occasionally needs a second pass to fully settle — S3 listing consistency
# can lag a fresh publish by a few seconds (seen live; harmless, just rerun
# this script again if all_pairs.csv looks incomplete right after a big batch
# finishes).
#
# Usage: ani-gather.sh --taxon GENUS:Yarrowia --compare SPECIES [--name yarrowia] [--method skani]
set -euo pipefail

NAMESPACE=ucr-stajichlab
POD=bfd-nextflow-head

TAXON=""
COMPARE="GENUS"
METHOD="skani"
NAME=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --taxon) TAXON="$2"; shift 2 ;;
    --compare) COMPARE="$2"; shift 2 ;;
    --method) METHOD="$2"; shift 2 ;;
    --name) NAME="$2"; shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ -z "$TAXON" ]]; then
  echo "Usage: $0 --taxon RANK:VALUE --compare LEVEL [--name slug] [--method skani|mash|sourmash|fastani]" >&2
  exit 1
fi

if [[ -z "$NAME" ]]; then
  NAME=$(echo "$TAXON" | sed -E 's/.*://' | tr '[:upper:]' '[:lower:]' | tr -c 'a-z0-9' '-' | sed -E 's/-+$//')
fi

RUN_DIR="/root/runs/${NAME}-gather"
LOG_FILE="/workspace/logs/cli-runs/${NAME}-gather.log"
PARAMS_FILE="/workspace/runs/${NAME}/params.yaml"

echo "==> [${NAME}] gathering (foreground — this is cheap, seconds to low minutes)"
kubectl exec -n "$NAMESPACE" "$POD" -- sh -c "
  mkdir -p '${RUN_DIR}' \$(dirname '${LOG_FILE}')
  cp /workspace/repo/samples.csv '${RUN_DIR}/samples.csv'
  cd '${RUN_DIR}'
  nextflow run /workspace/repo/nextflow/run_ani_gather.nf \
    -c /workspace/repo/nextflow/nextflow.config \
    -profile compare_ani_k8s \
    -params-file '${PARAMS_FILE}' \
    -resume 2>&1 | tee '${LOG_FILE}'
"
echo "==> [${NAME}] done. Results: s3://stajichlab/BFD/results/ANI/${METHOD}/${COMPARE}/"
