# Singularity / Apptainer Images for Nextflow Runs

## How Nextflow Requests Singularity Images

### 1. The `container` process directive (Nextflow-native)

Add a `container` directive to each process (or a `withLabel`/`withName` block in config),
and set `singularity.enabled = true` in `nextflow.config`:

```groovy
// nextflow.config — container blocks are env-driven; nothing to edit for a new run.
// params.singularity_cache is resolved from the environment (first match wins):
//   NXF_APPTAINER_CACHEDIR / NXF_SINGULARITY_CACHEDIR / APPTAINER_CACHEDIR
//   / SINGULARITY_CACHEDIR
apptainer  { enabled = false; cacheDir = params.singularity_cache }
singularity{ enabled = false; cacheDir = params.singularity_cache }
```

The shared cache on this cluster is `/bigdata/stajichlab/shared/singularity_cache`.
It is **never hardcoded** in this repo's configs on purpose — the `run_*.sh` launchers
export it for you (e.g. `run_funannotate.sh` defaults
`NXF_SINGULARITY_CACHEDIR`/`APPTAINER_CACHEDIR` to it, keeping a backwards-compat
fallback), and interactive `nextflow run` commands must export one of the env vars
above. The BFD and funannotate workflows fail loudly at startup with the full list
if none is set. Every `*_sif` path param derives from `params.singularity_cache`,
so there is exactly one place the cache location is decided.

`apptainer` and `singularity` are enabled/disabled per profile
(`apptainer.enabled = true` in `profile_BFD.config` / `profile_funannotate.config`,
or via the `-profile apptainer` / `-profile singularity` selectors).

Then add a `container` line in the pipeline `.nf` file or in a config block:

```groovy
// In the .nf file:
process PFAM_SCAN {
    label 'pfam'
    container 'docker://quay.io/biocontainers/hmmer:3.4--hdbdd923_1'
    ...
}

// Or in conf/profile_*.config (avoids touching .nf files):
process {
    withLabel: 'pfam' {
        container = 'docker://quay.io/biocontainers/hmmer:3.4--hdbdd923_1'
    }
}
```

When `singularity.enabled = true`, Nextflow automatically:
- Pulls the image and converts it to a `.sif` file in `cacheDir`
- Runs the process task inside that container
- Mounts the work directory into the container

The `container` string can be:
- `docker://org/image:tag` — pulls from Docker Hub and converts to SIF
- `oras://ghcr.io/org/image:tag` — pulls an OCI artifact directly
- `/absolute/path/to/image.sif` — uses a pre-built local SIF

### 2. Manual Singularity inside process scripts

Used when you need a persistent **instance** (e.g. a running database service) rather than
a one-shot container execution. `funannotate.nf` already does this for MariaDB:

```bash
# apptainer is loaded inline in the process script (module loads in beforeScript
# do not survive into the task on this cluster — see .living/learnings.md).
module load apptainer
apptainer instance start --writable-tmpfs \
    ${params.mariadb_sif} mysqldb_${asmid} /usr/bin/mysqld_safe
```

`params.mariadb_sif` is one of the `*_sif` params derived from
`params.singularity_cache` (`${params.singularity_cache}/mariadb-10.3.9.sif`).

The standard `container` directive only supports one-shot execution; use manual
`singularity instance start/stop` when the process needs a running service alongside it.

### 3. Current approach in this project (apptainer, `docker://` URIs resolved at run time)

The BFD and funannotate pipelines run inside apptainer containers referenced by
`docker://` URIs (pinned tag per tool, in `nextflow.config` / `profile_BFD.config` /
`profile_funannotate.config`), with `apptainer.enabled = true` turned on per profile.
The cache is warmed ahead of time by `nextflow/bin/setup_singularity_cache.sh` (Section
Option A) so runs never pull. Note it: apptainer/singularity are also referenced as
`docker://` biocontainer URIs whose *canonical module filename* is used for the
materialized `.sif`, so a `docker://quay.io/biocontainers/foo:tag` in config and a
`foo_tag.sif` in cache must match exactly what the setup script mirrors.

A few labels still use `conda` or `beforeScript` `module load` — check
`conf/profile_BFD.config` for which. To migrate one of those labels to a container:

1. Pre-pull the SIF into the shared cache (Option A above).
2. Replace the `beforeScript`/`conda` block for that label with a `container =` line
   in the same profile label block.
3. Make sure the label is not in the profile's `executor`/`docker`-only branches, and
   keep `apptainer.enabled = true`.

---

## Building Singularity / Apptainer Images

### Option A: Pull from BioContainers (easiest — no build needed)

Almost every bioinformatics tool has a pre-built image on
[quay.io/biocontainers](https://quay.io/organization/biocontainers).
Pre-pull to the shared cache to avoid pulling at run time. The maintained way to
warm the cache is the repo's setup script (mirrors the pinned module names that
`nextflow.config`, `profile_BFD.config`, and `profile_funannotate.config` reference):

```bash
# DRY-RUN first (echoes every `apptainer pull` it would run, does nothing):
nextflow/bin/setup_singularity_cache.sh --dry-run --cache-dir /bigdata/stajichlab/shared/singularity_cache
# Then apply the pulls:
nextflow/bin/setup_singularity_cache.sh --cache-dir /bigdata/stajichlab/shared/singularity_cache
```

Manual equivalent for one image:

```bash
module load apptainer
apptainer pull --dir /bigdata/stajichlab/shared/singularity_cache \
    docker://quay.io/biocontainers/hmmer:3.4--hdbdd923_1
```

This writes `hmmer_3.4--hdbdd923_1.sif` into the cache directory.

### Option B: Build from a Singularity Definition File (`.def`)

Write a definition file, e.g. `hmmer.def`:

```singularity
Bootstrap: docker
From: quay.io/biocontainers/hmmer:3.4--hdbdd923_1

%post
    # extra setup if needed
    apt-get -y install some-extra-dep 2>/dev/null || true

%environment
    export PATH=/usr/local/bin:$PATH
```

Build it (requires root or `--fakeroot` with user namespaces enabled on the cluster):

```bash
# With fakeroot (most HPC clusters allow this):
singularity build --fakeroot hmmer.sif hmmer.def

# With root (build node or interactive job on a root-allowed node):
sudo singularity build hmmer.sif hmmer.def
```

### Option C: Build from a Dockerfile via Podman → SIF

```bash
# 1. Build as an OCI image with Podman (rootless)
podman build -t myimage:latest .

# 2. Save to a tarball
podman save myimage:latest -o myimage.tar

# 3. Convert to SIF
singularity build myimage.sif docker-archive:myimage.tar
```

Or build directly from a running Podman daemon without saving:

```bash
singularity build myimage.sif docker-daemon:myimage:latest
```

### Option D: Sandbox build (interactive development / debugging)

```bash
# Create a writable sandbox directory from a base image
singularity build --sandbox myimage_sandbox/ docker://ubuntu:22.04

# Enter it and install / test tools interactively
singularity shell --writable myimage_sandbox/

# Convert the finished sandbox to a read-only SIF
singularity build myimage.sif myimage_sandbox/
```

---

## Wiring a New Container Into This Project

Cache location is `params.singularity_cache` (env-resolved, see top of this doc) —
do not hardcode `/bigdata/...` paths in any config. Minimal steps to add a container:

1. Pre-pull the SIF into the shared cache (Option A above, ideally by adding the
   pinned canonical filename to `nextflow/bin/setup_singularity_cache.sh`).
2. Add a `*_sif` param in `nextflow.config` (or the relevant `conf/profile_*.config`)
   as `"${params.singularity_cache}/<canonical>.sif"` — i.e. derive it from the cache
   param, don't hardcode the directory.
3. In `conf/profile_*.config`, replace the `beforeScript`/`conda` block for that label
   with a `container = params.<param>_sif` line.
4. Keep `apptainer.enabled = true` in the profile (or the `-profile apptainer` selector).

The `interproscan6` profile overrides `singularity_cache` with its own separate
containers dir (`conf/profile_interproscan6.config`) rather than the shared cache —
that override is intentional and should be left alone.
