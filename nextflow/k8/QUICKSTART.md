# compare_ANI on Nautilus — quick reference

> **Updated 2026-09-27 — runs are Jobs now; there is no head pod.** NRP prohibits
> idle interactive pods ("Running in interactive mode (`sleep infinity` command and
> manual start of computation) ... is prohibited, and user can be banned" —
> https://nrp.ai/documentation/userdocs/running/jobs/), so `head-pod.yaml` is gone.
> Each pipeline run is a Kubernetes Job whose command is `nextflow run`
> (`k8/bin/nf-job.sh`, used by `ani-run.sh`, `ani-gather.sh`, `ips6-run.sh`); it
> exits when the pipeline does. Launch dir, work dir and logs are on the PVC
> (`/workspace/runs/<name>`), so relaunching the same run resumes. Where the text
> below mentions the head pod, `kubectl exec` into it, `/root/runs`, or its 6 h
> cap, it describes the old setup.
>
> - Update the checkout: `bash k8/stage_repo.sh` (finite git Job; also PVC/RBAC/Secret).
> - Start / resume: `k8/bin/ani-run.sh --taxon GENUS:Yarrowia --compare SPECIES`
> - Follow: `k8/bin/ani-status.sh [ani-<name>]` or `kubectl logs -f -n ucr-stajichlab job/nf-ani-<name>`
> - Stop: `kubectl delete job -n ucr-stajichlab nf-ani-<name>`
> - Many taxa: `k8/bin/ani-suite.sh` runs at most `MAX_RUNS` (default 2) at once.
> - Task pods are `opportunistic` (preemptible) and avoid GPU nodes; at most 25
>   per run (NRP CPU-only guidance, https://nrp.ai/documentation/userdocs/running/cpu-only/).
> - The CephFS "launch dir can't be on the PVC" gotcha no longer holds: Nextflow
>   25.10.7 resumed fine from a launch dir on `bfd-work-pvc` (tested 2026-09-27).
> - The 6 h cap applies to pods without a controller. The run itself is a Job
>   now; task pods are still bare pods, so individual tasks keep the 6 h limit.

Terse operational cheat sheet. Full detail/rationale: `k8/README_compare_ani.md`.

## Before every session

```bash
kubectl get pods -n ucr-stajichlab --no-headers | grep '^nf-' | grep -v Running
```
Delete anything `Completed`/`Error` (safe — task pods are stateless):
```bash
kubectl delete pods -n ucr-stajichlab --field-selector=status.phase=Succeeded
```

Check the head pod is alive — **it dies on its own after 6h, no exceptions**
(cluster-enforced `activeDeadlineSeconds=21600`, verified live 2026-07-30):
```bash
kubectl get pod bfd-nextflow-head -n ucr-stajichlab
```
If it's gone / `Error` / `DeadlineExceeded`, recreate it (also re-syncs the repo
checkout — do this after every `git push` of pipeline changes too):
```bash
bash k8/stage_repo.sh
```

## Launch a run

```bash
k8/bin/ani-run.sh --taxon RANK:Value --compare RANK
```
`RANK` ∈ `PHYLUM SUBPHYLUM CLASS SUBCLASS ORDER FAMILY GENUS SPECIES`.
`--taxon` picks the subtree; `--compare` sets the rank pairwise groups are
formed at. Detached (`nohup` inside the head pod) — returns immediately.

**1. Single genus, compared at species level** (small first test):
```bash
k8/bin/ani-run.sh --taxon GENUS:Neurospora --compare SPECIES
```

**2. Broader — a whole order or family** (compare at genus level within it):
```bash
k8/bin/ani-run.sh --taxon ORDER:Pichiales    --compare GENUS
k8/bin/ani-run.sh --taxon FAMILY:Pichiaceae  --compare GENUS
```

**3. Later — a whole phylum** (this is the big one; expect hours, watch the
6h head-pod cap above — see README gotcha #3):
```bash
k8/bin/ani-run.sh --taxon PHYLUM:Ascomycota --compare ORDER
```

Batch several taxa at once (each launches independently, all inside the same
head pod):
```bash
k8/bin/ani-suite.sh --compare GENUS --taxa ORDER:Pichiales,ORDER:Saccharomycetales
```

## Check on it

```bash
k8/bin/ani-status.sh              # run logs + live head processes + task pods
k8/bin/ani-status.sh neurospora   # tail one run's log
```

## Get results (once compute looks done)

```bash
k8/bin/ani-gather.sh --taxon GENUS:Neurospora --compare SPECIES
```
Safe to re-run if `all_pairs.csv` looks incomplete right after a big batch —
S3 listing can lag a publish by a few seconds. Results land at
`s3://stajichlab/BFD/results/ANI/<method>/<compare-rank>/`.

## Pulling results and logs to your laptop

**Results are already on S3** — nothing to copy off the cluster for these,
just read them directly:
```bash
s3cmd ls s3://stajichlab/BFD/results/ANI/skani/GENUS/
# or
aws --endpoint-url https://s3-west.nrp-nautilus.io s3 sync \
  s3://stajichlab/BFD/results/ANI/skani/GENUS/Neurospora/ ./neurospora_results/
```

**Logs live on the PVC (`/workspace`)** — survive pod restarts, but need
`kubectl exec`/`kubectl cp` to reach your laptop:
```bash
# tail without copying
kubectl exec -n ucr-stajichlab bfd-nextflow-head -- \
  tail -100 /workspace/logs/cli-runs/neurospora.log

# pull the trace/report/timeline HTML+txt
kubectl cp ucr-stajichlab/bfd-nextflow-head:/workspace/logs/nextflow/ \
  ./nf-logs/ -c nextflow

# one failed task's .command.err — find it via `nextflow log` in the pod
# rather than hunting the work-dir hash by hand
kubectl exec -it -n ucr-stajichlab bfd-nextflow-head -- sh -c \
  'cd /root/runs/neurospora && nextflow log last -f hash,process,exit,workdir'
```

With `nf-job.sh` Jobs the resume cache is not in the run dir. Inside the running Job
pod `NXF_CACHE_DIR` is already set, so `nextflow log` works there. From another pod,
point it at a snapshot:
`NXF_CACHE_DIR=<run dir>/.nextflow-snapshots/<newest snap-*>/.nextflow nextflow log`.

**`.nextflow.log` (Nextflow's own engine log) lives only in the head pod's
container filesystem** — `/root/runs/<name>/.nextflow.log` — and is lost for
good if the pod dies (6h cap below) before you grab it:
```bash
kubectl cp ucr-stajichlab/bfd-nextflow-head:/root/runs/neurospora/.nextflow.log \
  ./neurospora.nextflow.log
```

Bottom line: results always land on S3 regardless of what happens to the head
pod, which is why `ani-gather.sh` is safe to re-run anytime — it's only the
logs that require a deliberate copy off the cluster.

## Known gotchas (see README for the "why")

1. `${projectDir}/bin/*.py` only resolves from a checkout *on the PVC* —
   don't launch from a laptop or ad hoc dir.
2. Nextflow's own launch dir (where you `cd` before `nextflow run`) must be
   **local to the head pod** (`/root/runs/...`), never on the PVC — CephFS
   can't do the file locking the resume-cache DB needs.
3. **The head pod is hard-capped at 6h by the cluster**, not by us. A batch
   spanning longer loses its in-flight run(s) and their `-resume` cache
   (already-published per-group results on S3 are safe). Watch
   `ani-status.sh` on long batches; recreate with `stage_repo.sh` +
   re-`ani-run.sh` if it dies.
4. `main.nf --pipeline compare_ani` currently fails to parse on this branch
   (unrelated missing funannotate module files) — always use
   `ani-run.sh`/`ani-gather.sh`, not the `main.nf` path in the main README.
