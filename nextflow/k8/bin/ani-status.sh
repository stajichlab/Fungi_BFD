#!/usr/bin/env bash
# k8/bin/ani-status.sh — check on ANI (and other nf-job.sh) runs.
#
# Usage:
#   ani-status.sh                 # pipeline Jobs + active task pods
#   ani-status.sh <name>          # tail one run's log (e.g. "ani-yarrowia")
set -euo pipefail

NAMESPACE="${NAMESPACE:-ucr-stajichlab}"
NAME="${1:-}"

if [[ -z "$NAME" ]]; then
  echo "=== Pipeline Jobs (one per ani-run.sh / ani-gather.sh / ips6-run.sh launch) ==="
  kubectl get jobs -n "$NAMESPACE" -l app=nf-run 2>/dev/null || echo "(none)"
  echo
  echo "=== Active task pods ==="
  kubectl get pods -n "$NAMESPACE" --no-headers 2>/dev/null | grep '^nf-' | grep -v -- '-[a-z0-9]\{5\}$' || true
  echo "($(kubectl get pods -n "$NAMESPACE" --no-headers 2>/dev/null | grep -c '^nf-[0-9a-f]\{32\}' || true) task pods)"
else
  JOB="nf-$(echo "$NAME" | tr '[:upper:]_' '[:lower:]-' | tr -cd 'a-z0-9-' | cut -c1-55 | sed -E 's/-+$//')"
  kubectl get job "$JOB" -n "$NAMESPACE" 2>/dev/null || echo "No Job $JOB (finished Jobs are kept 7 days)"
  echo "---"
  kubectl logs -n "$NAMESPACE" "job/$JOB" --tail=60 2>/dev/null \
    || echo "No log from job/$JOB. The PVC copy is /workspace/logs/cli-runs/${NAME}.log (read it with the shell pod)."
fi
