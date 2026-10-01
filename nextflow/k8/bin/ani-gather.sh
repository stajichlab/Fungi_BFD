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

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NAMESPACE=ucr-stajichlab

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

RUN_DIR="/workspace/runs/ani-${NAME}-gather"
# Same params as the compute run (ani-run.sh), fetched from its ConfigMap so
# the two phases can't drift apart.
PARAMS_FILE="$(mktemp -t "ani-${NAME}-gather.XXXX").yaml"
trap 'rm -f "$PARAMS_FILE"' EXIT
JOB_COMPUTE="nf-$(echo "ani-${NAME}" | tr '[:upper:]_' '[:lower:]-' | tr -cd 'a-z0-9-' | cut -c1-55 | sed -E 's/-+$//')"
if ! kubectl get configmap "${JOB_COMPUTE}-params" -n "$NAMESPACE" -o jsonpath='{.data.params\.yaml}' > "$PARAMS_FILE" || [ ! -s "$PARAMS_FILE" ]; then
  echo "ERROR: no params for the compute run (${JOB_COMPUTE}-params); run ani-run.sh --taxon ${TAXON} first" >&2
  exit 1
fi

echo "==> [${NAME}] gathering (foreground — this is cheap, seconds to low minutes)"
"${HERE}/nf-job.sh" --name "ani-${NAME}-gather" --run-dir "$RUN_DIR" --params "$PARAMS_FILE" --foreground -- \
  /workspace/repo/nextflow/run_ani_gather.nf \
  -c /workspace/repo/nextflow/nextflow.config \
  -profile compare_ani_k8s
echo "==> [${NAME}] done. Results: s3://stajichlab/BFD/results/ANI/${METHOD}/${COMPARE}/"
