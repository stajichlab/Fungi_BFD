#!/usr/bin/env bash
# k8/bin/ani-status.sh — check on CLI-launched ANI runs.
#
# Usage:
#   ani-status.sh                 # list all run logs + active task pods
#   ani-status.sh <name>          # tail one run's log (e.g. "yarrowia")
set -euo pipefail

NAMESPACE=ucr-stajichlab
POD=bfd-nextflow-head

if [[ $# -eq 0 ]]; then
  echo "=== Run logs (/workspace/logs/cli-runs/) ==="
  kubectl exec -n "$NAMESPACE" "$POD" -- sh -c 'ls -la /workspace/logs/cli-runs/ 2>/dev/null || echo "(none yet)"'
  echo
  echo "=== Active nextflow head processes (one per in-flight ani-run.sh/ani-gather.sh) ==="
  kubectl exec -n "$NAMESPACE" "$POD" -- sh -c "ps aux | grep '[n]extflow run' || echo '(none running)'"
  echo
  echo "=== k8s task pods (skani/mash/sourmash/fastani/report — one per group) ==="
  kubectl get pods -n "$NAMESPACE" --no-headers 2>/dev/null | grep '^nf-' || echo "(none)"
else
  NAME="$1"
  echo "=== Tail of ${NAME}.log ==="
  kubectl exec -n "$NAMESPACE" "$POD" -- tail -n 60 "/workspace/logs/cli-runs/${NAME}.log" 2>/dev/null \
    || echo "No log at /workspace/logs/cli-runs/${NAME}.log — check the name (see ani-status.sh with no args)"
fi
