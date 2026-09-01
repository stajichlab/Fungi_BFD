# interproscan6 on Nautilus (Kubernetes)

Files in this directory adapt `nextflow/interproscan6.nf` to run on the
[NRP Nautilus](https://nationalresearchplatform.org/) Kubernetes cluster.

## Files

| File | Purpose |
|---|---|
| `interproscan6_k8s.nf` | Drop into `nextflow/`; replaces nested singularity NF run with direct `interproscan.sh` call |
| `conf/profile_interproscan6_k8s.config` | Drop into `nextflow/conf/`; K8s executor, Docker, PVC paths |
| `k8s/pvc.yaml` | One-time PVC creation on Nautilus |
| `k8s/stage_data.sh` | One-time data staging script (IPS6 data + genomes → PVC) |
| `nextflow.config.patch` | One-line diff to add the profile to `nextflow/nextflow.config` |

## Key change from HPC version

The HPC `interproscan6.nf` runs a **nested Nextflow pipeline** inside each Slurm job:
```
module load singularity
nextflow run ebi-pf-team/interproscan6 -profile singularity,local ...
```

This does not work on Kubernetes because:
- No Lmod / `module load`
- Nested Nextflow inside a pod cannot submit further pods
- Singularity is not available (containerd is the runtime)

The K8s version instead calls `interproscan.sh` directly inside the official
Docker image (`docker.io/interpro/interproscan:6.0.0`), which is the image's
own entrypoint. One pod = one genome.

## Setup (one-time)

```bash
# 1. Apply PVC
kubectl apply -f k8s/pvc.yaml -n stajichlab

# 2. Stage data onto PVC (run from HPC login node with Nautilus kubeconfig)
bash k8s/stage_data.sh \
    --namespace stajichlab \
    --pvc stajichlab-pvc \
    --ips6-data /bigdata/stajichlab/shared/lib/interproscan6_data \
    --genomes genome_annotation/ \
    --samples samples.csv

# 3. Patch nextflow.config (add one profile line — see nextflow.config.patch)
```

## Running

```bash
nextflow run nextflow/interproscan6_k8s.nf \
    -c nextflow/nextflow.config \
    -profile interproscan6_k8s \
    -resume
```

Stub/dry run to verify channel logic without submitting pods:
```bash
nextflow run nextflow/interproscan6_k8s.nf \
    -c nextflow/nextflow.config \
    -profile interproscan6_k8s \
    -stub-run --n_test 2
```

## Tuning

- **`queueSize`** in the profile config controls max concurrent pods. Start with
  `10`, check Nautilus quota (`kubectl describe resourcequota -n stajichlab`),
  then increase. Each pod requests 16 CPU + 64 GB.
- **`cleanup = false`** is intentional initially — keep it off until `-resume`
  is confirmed working across sessions, then flip to `true` to reclaim PVC space.
- **`pullPolicy = IfNotPresent`** avoids re-pulling the ~8 GB IPS6 image on every
  pod. First pod on each node will pull; subsequent ones use the cached image.

## Retrieving results

```bash
# Spin up a loader pod (same as staging) and kubectl cp results back
kubectl cp stajichlab/bfd-data-loader:/workspace/genome_annotation/ ./genome_annotation/
```

Or use `rclone` with the Nautilus CEPH S3 endpoint if you set that up.
