# compare_ANI on Nautilus — quick reference

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
