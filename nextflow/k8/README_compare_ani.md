# compare_ANI on Nautilus (Kubernetes)

Files here adapt `workflows/compare_ANI.nf` (run via `main.nf --pipeline compare_ani`)
to run on the [NRP Nautilus](https://nrp.ai/) Kubernetes cluster.

## Files

| File | Purpose |
|---|---|
| `profile_compare_ani_k8s.config` | Drop into `nextflow/conf/`; k8s executor, Docker images, PVC paths |
| `stage_data_ani.sh` | One-time staging: repo checkout + genomes + samples.csv → PVC |
| `pvc.yaml` | Reused from the IPS6 setup — same PVC, ANI lives under `storageSubPath: ani` |

No new `.nf` workflow file is needed. Every ANI process already runs a BioContainers
OCI image with no `module load` and no nested Nextflow call, so `workflows/compare_ANI.nf`
runs unchanged under the k8s executor — unlike `interproscan6_k8s.nf`, which had to be
rewritten because the HPC version shells out to `module load singularity` and launches a
*second* Nextflow run inside each Slurm job (neither works inside a pod).

## The one real gotcha: `${projectDir}/bin/*.py`

`process.shell = ['/bin/bash', '-l']` is set globally in `nextflow.config`. The login
shell wipes the PATH Nextflow would otherwise use to auto-stage its `bin/` directory
into each task, so six ANI report/component scripts call themselves via the absolute
path `${projectDir}/bin/<script>.py` instead of a bare filename:

- `combine_ani_table.py`, `pick_representative_strain.py`, `report_ani.py`,
  `report_query_ani.py`, `combine_query_calls.py`, `mash_components.py`,
  `sourmash_matrix_to_long.py`

`${projectDir}` is interpolated **once**, on the machine that launches Nextflow, into
the literal path baked into every task's rendered `.command.sh`. On Slurm that's a
`/bigdata/...` path all nodes see over NFS — fine. On k8s it is **not** fine unless that
same absolute path also exists inside every pod. `stage_data_ani.sh` step 3 handles
this by copying the whole `nextflow/` checkout onto the PVC; you then have to **launch
Nextflow from that staged copy**, not from your `/bigdata` working copy, so
`${projectDir}` resolves identically on both sides.

(The alternative — dropping the `-l` from `process.shell` and switching the six
scripts to bare filenames so Nextflow's normal `bin/`-injection works again — is a
real code change across profiles shared with BFD/funannotate on HPC, so it's out of
scope for a first k8s pilot. Worth doing later if k8s becomes the primary path.)

## Setup (one-time)

```bash
# 1. PVC already exists from the IPS6 setup? Skip. Otherwise:
kubectl apply -f k8/pvc.yaml -n stajichlab

# 2. Stage repo + genomes + samples.csv onto the PVC
bash k8/stage_data_ani.sh \
    --namespace stajichlab \
    --pvc stajichlab-pvc \
    --repo nextflow/ \
    --genomes input/dna \
    --samples samples.csv

# 3. Register the profile (one line — see the diff below)
```

Add to `nextflow/nextflow.config`, in the existing `profiles { }` block:

```diff
 profiles {
     BFD           { includeConfig 'conf/profile_BFD.config'           }
     ...
     ani           { includeConfig 'conf/profile_ANI.config'           }
+    compare_ani_k8s { includeConfig 'conf/profile_compare_ani_k8s.config' }
     test          { includeConfig 'conf/test.config'                  }
 }
```

(Copy `profile_compare_ani_k8s.config` from `k8/` into `nextflow/conf/` first, or point
`includeConfig` at its path in `k8/` directly — either works.)

## Running

Because of the `${projectDir}/bin` gotcha above, run from the **staged** checkout, not
your `/bigdata` working copy. Easiest path: a Nautilus client mount (`ceph-fuse`, or any
node with the PVC bind-mounted at `/workspace`) or a throwaway pod with a shell and the
PVC mounted, matching `stage_data_ani.sh`'s loader pod but kept alive:

```bash
cd /workspace/nextflow/..     # wherever the PVC's ani/ subpath lands locally
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile compare_ani_k8s \
    --pipeline compare_ani \
    --compare GENUS \
    -resume
```

Stub/dry run first to check channel logic without submitting real pods:

```bash
nextflow run nextflow/main.nf \
    -c nextflow/nextflow.config \
    -profile compare_ani_k8s,stub \
    --pipeline compare_ani \
    -stub-run
```

## Tuning

- **`queueSize`** (30 by default here) caps concurrent pods — check
  `kubectl describe resourcequota -n stajichlab` before raising it.
- **`ani_batch_size`** / **`fastani_prefilter`**: same params as the HPC profile.
  For a first k8s pilot, prefer `skani` (default `ani_method`) over `fastani` —
  skani's per-group job count is much lower, which matters more on k8s where each
  task is a full pod-schedule round trip, not a cheap Slurm array step.
- **`cleanup = false`** initially, same reasoning as the IPS6 setup: keep work dirs
  until `-resume` is confirmed working, then flip to `true`.
- **Preemption**: `errorStrategy` retries once on exit 137/143 (pod eviction/OOM).
  ANI tasks are short, so this tolerates Nautilus's opportunistic scheduling well —
  much better than a multi-hour funannotate/IPS6 job would.

## What this doesn't solve for you

- Namespace/serviceAccount access on Nautilus — request that from the NRP admins,
  edit the `k8s { namespace = ...; serviceAccount = ... }` block once you have it.
- Actually confirming the PVC's storage class (`rook-cephfs`) and quota in your
  specific namespace — `kubectl get pvc`, `kubectl describe resourcequota`.
- Running/validating any of this against the live cluster — these files are
  unexercised until you (or an agent with your kubeconfig) run them against Nautilus.
