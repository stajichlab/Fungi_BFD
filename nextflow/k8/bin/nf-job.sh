#!/usr/bin/env bash
# k8/bin/nf-job.sh — run one Nextflow pipeline on Nautilus as a Kubernetes Job.
#
# NRP prohibits idle interactive pods ("sleep infinity ... can be banned",
# https://nrp.ai/documentation/userdocs/running/jobs/), so there is no head pod:
# each run is a Job whose command IS `nextflow run`, and it exits when the
# pipeline does. Everything that must survive lives on the PVC (bfd-work-pvc,
# /workspace): the repo checkout (/workspace/repo), the run's launch dir
# (.nextflow resume cache — verified to work on rook-cephfs with Nextflow
# 25.10.7), workDir, logs. Relaunching the same --name resumes (-resume).
#
# The params file is passed from your laptop as a ConfigMap (nf-<name>-params),
# mounted at /config/params.yaml.
#
# Usage:
#   nf-job.sh --name NAME --run-dir /workspace/runs/NAME --params local.yaml \
#             [--cvmfs] [--foreground] -- <nextflow run args...>
# e.g.
#   nf-job.sh --name ani-yarrowia --run-dir /workspace/runs/ani-yarrowia \
#             --params /tmp/p.yaml -- /workspace/repo/nextflow/run_ani_compute.nf \
#             -c /workspace/repo/nextflow/nextflow.config -profile compare_ani_k8s
#
#   --cvmfs       also mount the CVMFS PVC at /cvmfs (InterProScan 6 data)
#   --foreground  wait for the Job to finish and stream its log
#
# Watch / stop:
#   kubectl logs -f -n ucr-stajichlab job/nf-NAME   (also /workspace/logs/cli-runs/NAME.log)
#   kubectl delete job -n ucr-stajichlab nf-NAME    (Nextflow gets SIGTERM and cleans up its task pods)
set -euo pipefail

NAMESPACE="${NAMESPACE:-ucr-stajichlab}"
NXF_IMAGE="${NXF_IMAGE:-nextflow/nextflow:25.10.7}"
NRP_PROJECT="${NRP_PROJECT:-stajichlab-fungi-bfd}"

NAME="" RUN_DIR="" PARAMS="" CVMFS=0 FOREGROUND=0
while [[ $# -gt 0 ]]; do
  case "$1" in
    --name) NAME="$2"; shift 2 ;;
    --run-dir) RUN_DIR="$2"; shift 2 ;;
    --params) PARAMS="$2"; shift 2 ;;
    --cvmfs) CVMFS=1; shift ;;
    --foreground) FOREGROUND=1; shift ;;
    --) shift; break ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done
if [[ -z "$NAME" || -z "$RUN_DIR" || -z "$PARAMS" || $# -eq 0 ]]; then
  echo "Usage: $0 --name NAME --run-dir DIR --params FILE [--cvmfs] [--foreground] -- <nextflow run args>" >&2
  exit 1
fi
[[ -f "$PARAMS" ]] || { echo "ERROR: params file $PARAMS not found" >&2; exit 1; }

# Kubernetes names: lowercase alphanumerics and '-', max 63 chars.
JOB="nf-$(echo "$NAME" | tr '[:upper:]_' '[:lower:]-' | tr -cd 'a-z0-9-' | cut -c1-55 | sed -E 's/-+$//')"

# One Job per name: refuse to start a second copy of a running run; replace a
# finished one (a Job's pod template is immutable, so relaunch = delete + create).
if kubectl get job "$JOB" -n "$NAMESPACE" >/dev/null 2>&1; then
  if [[ "$(kubectl get job "$JOB" -n "$NAMESPACE" -o jsonpath='{.status.active}')" == "1" ]]; then
    echo "ERROR: $JOB is still running. Stop it first: kubectl delete job -n $NAMESPACE $JOB" >&2
    exit 1
  fi
  kubectl delete job "$JOB" -n "$NAMESPACE" --wait=true >/dev/null
fi

kubectl create configmap "${JOB}-params" -n "$NAMESPACE" --from-file=params.yaml="$PARAMS" \
  --dry-run=client -o yaml | kubectl apply -f - >/dev/null

# Quote the nextflow args for embedding in the container's bash -c script.
NF_ARGS=$(printf ' %q' "$@")

CVMFS_ENV="" CVMFS_MOUNT="" CVMFS_VOL=""
if [[ $CVMFS -eq 1 ]]; then
  CVMFS_ENV='
            - {name: CVMFS_SERVER_URL, value: "http://www.ebi.ac.uk/cvmfs/geo2.embl.de"}
            - {name: CVMFS_REPOSITORIES, value: "geo2.embl.de"}
            - {name: CVMFS_HTTP_PROXY, value: "DIRECT"}
            - {name: CVMFS_KEYS_DIR, value: "/etc/cvmfs/keys.d"}'
  # HostToContainer: the CSI driver's FUSE mount must propagate into the pod.
  CVMFS_MOUNT='
            - {name: my-cvmfs, mountPath: /cvmfs, mountPropagation: HostToContainer}'
  CVMFS_VOL='
        - name: my-cvmfs
          persistentVolumeClaim: {claimName: cvmfs}'
fi

kubectl apply -f - <<EOF
apiVersion: batch/v1
kind: Job
metadata:
  name: ${JOB}
  namespace: ${NAMESPACE}
  labels: {app: nf-run, run: ${JOB}}
spec:
  # Retries cover node reboots / eviction of the Nextflow pod; each resumes.
  backoffLimit: 2
  ttlSecondsAfterFinished: 604800
  template:
    metadata:
      labels: {app: nf-run, run: ${JOB}}
      annotations: {nrp-nautilus.io/project: ${NRP_PROJECT}}
    spec:
      restartPolicy: Never
      serviceAccountName: nextflow-runner
      containers:
        - name: nextflow
          image: ${NXF_IMAGE}
          command: ["bash", "-c"]
          args:
            - |
              set -euo pipefail
              mkdir -p '${RUN_DIR}' /workspace/logs/cli-runs
              cd '${RUN_DIR}'
              [ -f samples.csv ] || cp /workspace/repo/samples.csv .
              LOG=/workspace/logs/cli-runs/${NAME}.log
              [ -f "\$LOG" ] && mv "\$LOG" "\$LOG.\$(date -u +%Y%m%dT%H%M%S)"
              echo "[\$(date -u)] ${JOB}: nextflow run${NF_ARGS} -params-file /config/params.yaml -resume" | tee "\$LOG"
              nextflow run${NF_ARGS} -params-file /config/params.yaml -resume 2>&1 | tee -a "\$LOG"
          env:
            - {name: NXF_HOME, value: /workspace/.nextflow}
            - name: AWS_ACCESS_KEY_ID
              valueFrom: {secretKeyRef: {name: nrp-s3-creds, key: AWS_ACCESS_KEY_ID}}
            - name: AWS_SECRET_ACCESS_KEY
              valueFrom: {secretKeyRef: {name: nrp-s3-creds, key: AWS_SECRET_ACCESS_KEY}}${CVMFS_ENV}
          # NRP admission: cpu limit/request ratio must be <= 1.2.
          resources:
            requests: {cpu: "2", memory: 4Gi}
            limits: {cpu: "2", memory: 4Gi}
          volumeMounts:
            - {name: work, mountPath: /workspace}
            - {name: params, mountPath: /config}${CVMFS_MOUNT}
      volumes:
        - name: work
          persistentVolumeClaim: {claimName: bfd-work-pvc}
        - name: params
          configMap: {name: ${JOB}-params}${CVMFS_VOL}
EOF

echo "==> ${JOB} started (run dir ${RUN_DIR})"
echo "    follow: kubectl logs -f -n ${NAMESPACE} job/${JOB}"
echo "    stop:   kubectl delete job -n ${NAMESPACE} ${JOB}"
if [[ $FOREGROUND -eq 1 ]]; then
  kubectl wait --for=condition=Ready pod -l job-name="${JOB}" -n "$NAMESPACE" --timeout=600s >/dev/null
  kubectl logs -f -n "$NAMESPACE" "job/${JOB}"
  kubectl wait --for=condition=complete "job/${JOB}" -n "$NAMESPACE" --timeout=60s >/dev/null
fi
