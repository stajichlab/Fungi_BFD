# compare_ANI on Nautilus (Kubernetes)

Files here adapt `workflows/compare_ANI.nf` (run via `main.nf --pipeline compare_ani`)
to run on the [NRP Nautilus](https://nrp.ai/) Kubernetes cluster, in the
`ucr-stajichlab` namespace. **Live-tested**: real skani comparisons of real
genomes staged from `s3://stajichlab/BFD`, through real k8s Job pods,
publishing results back to S3 — not just written and hoped for.

## Files

| File | Purpose |
|---|---|
| `profile_compare_ani_k8s.config` | The `compare_ani_k8s` profile (see `nextflow.config`): k8s executor, Docker images, PVC + S3 paths |
| `pvc.yaml` | `bfd-work-pvc` — dedicated RWX PVC (`rook-cephfs`) for workDir + the repo checkout |
| `rbac.yaml` | ServiceAccount `nextflow-runner` + Role/RoleBinding it needs to submit/watch task pods |
| `head-pod.yaml` | Long-lived pod `nextflow run` is actually launched from |
| `stage_repo.sh` | One-time (and re-run-safe): applies the above, creates the S3-credentials Secret, clones/updates this repo onto the PVC |

No new `.nf` workflow file is needed. Every ANI process already runs a BioContainers
OCI image with no `module load` and no nested Nextflow call, so `workflows/compare_ANI.nf`
runs unchanged under the k8s executor — unlike `interproscan6_k8s.nf`, which had to be
rewritten because the HPC version shells out to `module load singularity` and launches a
*second* Nextflow run inside each Slurm job (neither works inside a pod).

## Why the head process runs inside the cluster, not your laptop

`kuberun` (which let you drive a k8s-executor run from an outside machine while it
submitted pods into the cluster) is deprecated and unmaintained upstream. Without
Fusion (a paid Seqera add-on, not used here), the k8s executor's `workDir` must be a
POSIX volume mounted into both the head process and every task pod — that only works
inside the cluster. So `nextflow run` is launched from inside `head-pod.yaml`, which
you reach via `kubectl exec`.

Genome inputs and results don't need this treatment: they stay at
`s3://stajichlab/BFD/...` as plain params (see `profile_compare_ani_k8s.config`'s
`aws {}` block), and Nextflow's own head process stages those to/from the PVC-backed
workDir directly — task containers never need AWS credentials or the AWS CLI baked in.

## Two real gotchas, both verified live

**1. `${projectDir}/bin/*.py`.** `process.shell = ['/bin/bash', '-l']` is set globally
in `nextflow.config`. The login shell wipes the PATH Nextflow would otherwise use to
auto-stage its `bin/` directory into each task, so six ANI report/component scripts
call themselves via the absolute path `${projectDir}/bin/<script>.py` instead of a
bare filename: `combine_ani_table.py`, `pick_representative_strain.py`,
`report_ani.py`, `report_query_ani.py`, `combine_query_calls.py`,
`mash_components.py`, `sourmash_matrix_to_long.py`.

`${projectDir}` is interpolated **once**, on whatever launched Nextflow — the *head*
process — into the literal path baked into every task's rendered `.command.sh`. On
Slurm that's a `/bigdata/...` path all nodes see over NFS — fine. On k8s it's only
fine if `${projectDir}` also exists, at the identical path, inside every task pod.
Confirmed live: this requires launching from a checkout of this repo that lives
*on the same PVC* mounted into every pod (`stage_repo.sh` clones it to
`/workspace/repo`) — not from a laptop or `/bigdata` copy. Also confirmed: you must
launch via `main.nf`, not by invoking `workflows/compare_ANI.nf` directly — doing the
latter anchors `${projectDir}` to `nextflow/workflows/` instead of `nextflow/`, and
`bin/` isn't there.

**2. Nextflow's launch directory can't be on the PVC.** `rook-cephfs` doesn't support
the file locking Nextflow's resume-cache database (LevelDB) needs — `nextflow run`
fails outright with "Can't open cache DB ... needs a shared file system that supports
file locks" if `.nextflow/` (created wherever you `cd` before running `nextflow run`)
is on the PVC. Fix: launch from a local directory inside the head pod's own container
filesystem (e.g. `/root/runs/<name>`, as below) — only `workDir` (set to
`/workspace/work/ANI` in the profile) needs to be on the PVC. Trade-off: `-resume`
state is lost if the head pod is deleted/recreated, even though `workDir` contents on
the PVC survive — a fresh pod would recompute from scratch rather than pick up where a
prior pod left off.

## Setup (one-time)

```bash
export S3CFG=~/.s3cfg   # NRP User Portal -> S3 Tokens page, west pool
bash k8/stage_repo.sh
```

This applies `pvc.yaml` and `rbac.yaml`, creates the `nrp-s3-creds` Secret from your
`.s3cfg` (never committed — credentials only ever live in the k8s Secret), applies
`head-pod.yaml`, and clones/updates this repo onto the PVC at `/workspace/repo`. It
prints the exact `kubectl exec` command to launch a run when it finishes.

## Running

```bash
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
```

`-params-file` needs to be reachable inside the pod too — `kubectl cp` it into
`/root/runs/ANI/` (or onto the PVC) first if it's not already part of the repo
checkout.

## Tuning

- **`queueSize`** (30 by default here) caps concurrent pods — check
  `kubectl describe resourcequota -n ucr-stajichlab` before raising it.
- **`ani_batch_size`** / **`fastani_prefilter`**: same params as the HPC profile.
  For a first k8s pilot, prefer `skani` (default `ani_method`) over `fastani` —
  skani's per-group job count is much lower, which matters more on k8s where each
  task is a full pod-schedule round trip, not a cheap Slurm array step.
- **`cleanup = false`** initially, same reasoning as the IPS6 setup: keep work dirs
  until `-resume` is confirmed working, then flip to `true`.
- **Preemption**: `errorStrategy` retries once on exit 137/143 (pod eviction/OOM).
  ANI tasks are short, so this tolerates Nautilus's opportunistic scheduling well —
  much better than a multi-hour funannotate/IPS6 job would.
- **Image caching**: `k8s.pullPolicy = 'IfNotPresent'` means a node that already
  pulled an image for one task reuses it for later tasks on that same node — but
  that's per-node, not cluster-wide; a fresh node still pulls once. There's no
  Singularity-style shared SIF cache to configure here: task pods run plain OCI
  images via containerd (no evidence Nautilus supports running `.sif` files
  directly as a pod's container), so image caching is the container runtime's job,
  not something a namespace user configures further.

## What this doesn't solve for you

- `main.nf` currently fails to parse on this branch regardless of `--pipeline`,
  because it eagerly includes every workflow file and
  `subworkflows/local/FUNANNOTATE_GENOME_PREP.nf` references three
  `modules/funannotate/genome/*/main.nf` files that don't exist in the repo
  (a missed-`git add` from the funannotate decomposition commit, unrelated to
  k8s — recoverable from the pre-decomposition `funannotate.nf`, not yet done).
- Confirming PVC quota/storage headroom at real pipeline scale — `bfd-work-pvc` is
  requested at 200Gi; only validated with a 2-genome smoke comparison so far.
