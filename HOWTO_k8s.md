# Running Nextflow Jobs on Kubernetes with Container Images

## Key Difference: K8s Uses OCI/Docker Images, Not SIF

Kubernetes has its own container runtime (containerd or CRI-O) and pulls **Docker/OCI images**
directly — it does not use `.sif` files. The same `container = 'docker://quay.io/...'` URIs
used with Singularity work on k8s too, just through a different runtime. If the k8s nodes have
Apptainer installed as the container runtime (uncommon but possible in HPC k8s clusters), SIF
files can work, but the standard path is OCI images.

---

## Step 1: Kubeconfig Access

You need a `kubeconfig` pointing at the remote cluster. The remote system admin provides this;
it typically lives at `~/.kube/config` or set it explicitly:

```bash
export KUBECONFIG=/path/to/remote-cluster.yaml
```

Test connectivity: `kubectl get nodes`

---

## Step 2: Shared Storage — PersistentVolumeClaim

Nextflow's k8s executor requires a **PersistentVolumeClaim (PVC)** that all pods can read/write.
This replaces the local `work/` directory. The k8s admin creates this, or you create it:

```yaml
# pvc.yaml
apiVersion: v1
kind: PersistentVolumeClaim
metadata:
  name: fungi-bfd-work
spec:
  accessModes: [ReadWriteMany]   # must be RWX for multi-pod access
  resources:
    requests:
      storage: 500Gi
  storageClassName: nfs-client    # whatever your cluster provides
```

```bash
kubectl apply -f pvc.yaml
```

Your input data also needs to be reachable — either on the PVC, an NFS mount the pods can see,
or in object storage (S3/GCS) via Fusion (see Step 4).

---

## Step 3: Add a `k8s` Profile to `nextflow.config`

```groovy
// nextflow.config — add inside profiles { ... }
profiles {
    // ... existing profiles ...

    k8s {
        process.executor = 'k8s'

        k8s {
            namespace        = 'stajichlab'         // k8s namespace
            serviceAccount   = 'nextflow'            // service account with pod/PVC permissions
            storageClaimName = 'fungi-bfd-work'      // PVC name from above
            storageMountPath = '/workspace'          // where PVC mounts inside pods
            // storageSubPath = 'fungi_bfd'          // optional: subdirectory on the PVC
        }

        workDir = '/workspace/work'                  // must be on the PVC mount

        // Container images for each process label.
        // These are OCI images — same tool versions as the Singularity URIs but no
        // 'docker://' prefix; k8s pulls them natively via containerd/CRI-O.
        process {
            withLabel: 'pfam' {
                container = 'quay.io/biocontainers/hmmer:3.4--hdbdd923_1'
                cpus      = 4
                memory    = '32 GB'
            }
            withLabel: 'cazy' {
                container = 'quay.io/biocontainers/dbcan:4.1.4--pyhdfd78af_0'
                cpus      = 16
                memory    = '8 GB'
            }
            withLabel: 'merops' {
                container = 'quay.io/biocontainers/blast:2.16.0--h66d330f_2'
                cpus      = 8
                memory    = '4 GB'
            }
            // add further labels as needed
        }
    }
}
```

---

## Step 4 (Optional): Fusion Filesystem for S3/Object Storage

If the remote k8s cluster can reach S3 or GCS, Nextflow's **Fusion** filesystem lets pods
read/write object storage as if it were local — no NFS PVC needed for the work directory:

```groovy
k8s {
    // namespace / serviceAccount as above
    // omit storageClaimName / storageMountPath when using Fusion
}

fusion {
    enabled = true
}

wave {
    enabled = true    // Fusion requires Wave for container augmentation
}

workDir = 's3://your-bucket/fungi-bfd/work'
```

Fusion + Wave require a Nextflow license key for production use (free for open-source /
academic projects).

---

## Step 5: Splitting Jobs Between SLURM and K8s

Route specific processes to k8s while keeping the rest on SLURM by using a combined profile.
Add a new config file `conf/profile_k8s_hybrid.config`:

```groovy
// conf/profile_k8s_hybrid.config
process {
    // Default: stay on SLURM (inherited from nextflow.config executor block)

    // Move selected labels to k8s:
    withLabel: 'pfam' {
        executor  = 'k8s'
        container = 'quay.io/biocontainers/hmmer:3.4--hdbdd923_1'
    }
    withLabel: 'cazy' {
        executor  = 'k8s'
        container = 'quay.io/biocontainers/dbcan:4.1.4--pyhdfd78af_0'
    }
}

k8s {
    namespace        = 'stajichlab'
    serviceAccount   = 'nextflow'
    storageClaimName = 'fungi-bfd-work'
    storageMountPath = '/workspace'
}
```

Register it in `nextflow.config`:

```groovy
profiles {
    // ... existing profiles ...
    k8s_hybrid { includeConfig 'conf/profile_k8s_hybrid.config' }
}
```

Run with both profiles active:

```bash
nextflow run nextflow/BFD.nf -profile BFD,k8s_hybrid
```

---

## Step 6: Running

```bash
# Point at the remote cluster
export KUBECONFIG=~/.kube/remote-cluster.yaml

# Launch — Nextflow submits pods directly; no sbatch wrapper needed
nextflow run nextflow/BFD.nf \
    -profile k8s \
    --samples samples.csv \
    -resume
```

Nextflow itself runs on the local/head node and submits one pod per task. Progress monitoring,
`-resume`, HTML reports, and `nextflow log` all work the same way as with the SLURM executor.

---

## Practical Notes for This Project

| Concern | Detail |
|---|---|
| **Nextflow version** | Running 26.04.3 — k8s executor is fully supported; Fusion and Wave are available |
| **SignalP / AIUpred (GPU)** | K8s GPU scheduling uses resource limits. Set via `pod = [[resourceLimit: 'nvidia.com/gpu', value: '1']]` in the process block alongside `container =` |
| **MariaDB (funannotate)** | The `singularity instance start` pattern in `funannotate.nf` does not translate to k8s. Replace with a MariaDB sidecar container or a dedicated k8s Service |
| **Input data visibility** | Proteins/GFFs under `/bigdata/stajichlab/...` are not visible to k8s pods unless the PVC is backed by the same NFS share, or data is copied to S3/the PVC first. Confirm this with the k8s cluster admin before migrating jobs |
| **Container image tag pinning** | Always pin to an explicit tag (e.g. `hmmer:3.4--hdbdd923_1`) rather than `latest` for reproducibility |
| **`beforeScript` removal** | When migrating a process label to k8s, remove the `beforeScript` module-load block for that label — it has no effect inside a container and may cause errors |

The biggest practical hurdle is usually **data locality** — getting `input/` protein FASTAs
visible to the pods is the first thing to resolve with the k8s cluster admin.
