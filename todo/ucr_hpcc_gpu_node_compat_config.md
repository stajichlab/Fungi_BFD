# Extract UCR HPCC GPU node compatibility settings into a shared cluster config

**Status**: idea (not started)
**Priority**: medium
**Category**: infrastructure
**Added**: 2026-08-01
**Author**: Jason Stajich

## Motivation

While debugging a `RUN_IDP` (aiupred) failure in a downstream project
(`compare/Exophiala`, driven by this repo's `nextflow/` pipeline), found that
`short_gpu` spans 8 GPU generations on UCR HPCC — K80 (Kepler, sm_37), P100
(sm_60), A100 (sm_80), H100 (sm_90), Ada6000 (sm_89), and Blackwell6000
(sm_120). The `aiupred` module's PyTorch build (`2.6.0+cu124`) only ships
kernels for sm_50–sm_90, so jobs landed on K80 or Blackwell6000 nodes fail
with `CUDA error: no kernel image is available for execution on the device`.
Nextflow's built-in retry (`maxRetries=2`) isn't reliable protection — a run
against `compare/Exophiala` finished with 3 permanently failed `RUN_IDP` tasks
after exhausting retries, because retries don't change which node class gets
picked.

Fixed short-term by adding `--exclude=gpu01,gpu02,gpu03,gpu13,gpu14` to the
`idp` label's `clusterOptions` in `conf/profile_BFD.config`. Slurm node
features don't distinguish Blackwell from A100/H100/Ada6000 (all share
`gpu_latest`/`gpu_highmem`), so this had to be an explicit node exclusion
rather than a `-C` feature constraint.

## The idea

This exclusion list is a fact about **UCR HPCC's hardware topology**, not
about the BFD pipeline. It's already duplicated in spirit — both
`conf/profile_BFD.config` (`idp` label) and `conf/profile_funannotate.config`
(likely `signalp`, which also runs a GPU-dependent tool on `short_gpu`) will
need the same node-exclusion list, and it will silently drift if only one
profile gets updated the next time a node is added/retired or a tool's torch
build changes its supported compute-capability range.

Proposed structure: a `conf/hpcc_ucr.config` (or similar) defining a
reusable param, e.g. `params.gpu_exclude_legacy_blackwell =
'gpu01,gpu02,gpu03,gpu13,gpu14'`, included once from `nextflow.config`
alongside the existing profile includes. Each GPU-labeled process's
`clusterOptions` across all profile configs then interpolates that param
instead of hardcoding the node list.

## Open questions

- Does `signalp/6-gpu`'s torch/CUDA build have the same sm_50–90 constraint,
  or a different compatible range? Worth checking before assuming the same
  exclusion list applies.
- Should this live as a single flat exclude list, or a small compatibility
  matrix (per-tool minimum/maximum compute capability) that generates the
  right `--exclude` per label? Flat list is simpler and sufficient for the
  two known GPU-dependent tools today; a matrix is more correct if a third
  GPU tool with different constraints shows up.

## Not urgent

This is a refactor, not a bugfix — the immediate failure is already patched
in `profile_BFD.config`. Fine to defer until the next GPU-dependent process
needs the same treatment.
