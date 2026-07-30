#!/usr/bin/env bash
# k8/stage_repo.sh — one-time (and re-run-safe) setup for the compare_ANI k8s
# pilot: PVC, RBAC, S3-credentials Secret, head pod, and a repo checkout on
# the PVC that the head pod launches Nextflow from.
#
# samples.csv and all pipeline code come from `git clone`/`git pull` — nothing
# needs staging via `kubectl cp`. Genome inputs and results stay on S3
# (s3://stajichlab/BFD/...); only the repo checkout and Nextflow's own workDir
# live on the PVC, because the k8s executor's workDir must be a POSIX volume
# shared between the head process and every task pod (see the profile config
# for why). Run this from a machine with `kubectl` pointed at the Nautilus
# 'ucr-stajichlab' namespace.
set -euo pipefail

NAMESPACE=ucr-stajichlab
S3CFG="${S3CFG:-$HOME/.s3cfg}"
REPO_URL="${REPO_URL:-https://github.com/stajichlab/Fungi_BFD.git}"
REPO_BRANCH="${REPO_BRANCH:-$(git -C "$(dirname "$0")/../.." rev-parse --abbrev-ref HEAD 2>/dev/null || echo main)}"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "==> Applying PVC"
kubectl apply -f "$here/pvc.yaml"

echo "==> Applying RBAC (ServiceAccount + Role + RoleBinding)"
kubectl apply -f "$here/rbac.yaml"

echo "==> Creating/updating nrp-s3-creds Secret from $S3CFG"
if [ ! -f "$S3CFG" ]; then
  echo "ERROR: $S3CFG not found. Set S3CFG=/path/to/.s3cfg or create ~/.s3cfg" >&2
  echo "  (NRP User Portal -> S3 Tokens page for west-pool credentials)" >&2
  exit 1
fi
AK=$(awk -F'=' '/^access_key/{gsub(/ /,"",$2); print $2}' "$S3CFG")
SK=$(awk -F'=' '/^secret_key/{gsub(/ /,"",$2); print $2}' "$S3CFG")
if [ -z "$AK" ] || [ -z "$SK" ]; then
  echo "ERROR: could not parse access_key/secret_key from $S3CFG" >&2
  exit 1
fi
kubectl create secret generic nrp-s3-creds -n "$NAMESPACE" \
  --from-literal=AWS_ACCESS_KEY_ID="$AK" \
  --from-literal=AWS_SECRET_ACCESS_KEY="$SK" \
  --dry-run=client -o yaml | kubectl apply -f -
unset AK SK

echo "==> Applying head pod"
kubectl apply -f "$here/head-pod.yaml"
kubectl wait --for=condition=Ready pod/bfd-nextflow-head -n "$NAMESPACE" --timeout=120s

echo "==> Cloning/updating $REPO_URL (branch: $REPO_BRANCH) onto the PVC"
kubectl exec -n "$NAMESPACE" bfd-nextflow-head -- sh -c "
  set -e
  if [ -d /workspace/repo/.git ]; then
    cd /workspace/repo && git fetch origin && git checkout '$REPO_BRANCH' && git pull --ff-only
  else
    git clone --branch '$REPO_BRANCH' '$REPO_URL' /workspace/repo
  fi
"

echo "==> Done. Launch a run with e.g.:"
# NOTE: the launch dir must NOT be on the PVC (/workspace/...) — rook-cephfs
# doesn't support the file locking Nextflow's resume-cache DB needs, and
# `nextflow run` fails outright if .nextflow/ ends up there. Use a directory
# local to the head pod's own container filesystem instead (/root/runs/...),
# same as k8/README_compare_ani.md and the ani-*.sh wrappers do.
cat <<'EOF'
kubectl exec -it -n ucr-stajichlab bfd-nextflow-head -- sh -c '
  mkdir -p /root/runs/ANI && cd /root/runs/ANI
  cp /workspace/repo/samples.csv .
  nextflow run /workspace/repo/nextflow/main.nf \
    -c /workspace/repo/nextflow/nextflow.config \
    -profile compare_ani_k8s \
    --pipeline compare_ani \
    -params-file /path/to/params_ani.yaml \
    -resume
'
EOF
echo "NOTE: main.nf currently fails to parse on this branch regardless of"
echo "--pipeline (see k8/README_compare_ani.md's 'What this doesn't solve"
echo "for you'). Use k8/bin/ani-run.sh + k8/bin/ani-gather.sh instead until"
echo "that's fixed."
