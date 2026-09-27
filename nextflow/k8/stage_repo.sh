#!/usr/bin/env bash
# k8/stage_repo.sh — one-time (and re-run-safe) setup for the compare_ANI k8s
# pilot: PVC, RBAC, S3-credentials Secret, and a repo checkout on the PVC that
# the pipeline Jobs (k8/bin/nf-job.sh) run Nextflow from. There is no head pod:
# NRP prohibits idle `sleep infinity` pods, so the checkout is updated by a
# finite git Job, and each pipeline run is its own Job.
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

echo "==> Cloning/updating $REPO_URL (branch: $REPO_BRANCH) onto the PVC (finite Job)"
kubectl delete job repo-sync -n "$NAMESPACE" --ignore-not-found --wait=true >/dev/null
kubectl apply -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: repo-sync
  namespace: ${NAMESPACE}
spec:
  backoffLimit: 1
  ttlSecondsAfterFinished: 3600
  template:
    spec:
      restartPolicy: Never
      containers:
        - name: git
          image: alpine/git:2.45.2
          command: ["sh", "-c"]
          args:
            - |
              set -e
              if [ -d /workspace/repo/.git ]; then
                cd /workspace/repo && git fetch origin && git checkout '${REPO_BRANCH}' && git pull --ff-only
              else
                git clone --branch '${REPO_BRANCH}' '${REPO_URL}' /workspace/repo
              fi
              git -C /workspace/repo log --oneline -1
          resources:
            requests: {cpu: "1", memory: 1Gi}
            limits: {cpu: "1", memory: 1Gi}
          volumeMounts: [{name: work, mountPath: /workspace}]
      volumes:
        - name: work
          persistentVolumeClaim: {claimName: bfd-work-pvc}
EOF
kubectl wait --for=condition=complete job/repo-sync -n "$NAMESPACE" --timeout=300s >/dev/null
kubectl logs -n "$NAMESPACE" job/repo-sync | tail -1

echo "==> Done. Launch runs with k8/bin/ani-run.sh (or k8/bin/nf-job.sh directly), e.g.:"
echo "    k8/bin/ani-run.sh --taxon GENUS:Yarrowia --compare SPECIES"
