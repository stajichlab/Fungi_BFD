# Environments & Installations

This is a fungal genomics / annotation project (BFD — Fungi_BFD) run on the
**UCR HPCC SLURM cluster**. Most compute happens through Nextflow DSL2
pipelines (see `nextflow/`) and SLURM batch scripts, using environment
modules, `pixi`, and Singularity containers rather than a single project-wide
Python env.

## Primary Environment

- **Cluster**: UCR HPCC (SLURM), Rocky Linux 8.x
- **Default `python3`**: 3.9 (system miniconda, `py39_4.12.0`)
- **Newer interpreter available**: `python3.12` at `/usr/bin/python3.12`

## Mycelium tooling — IMPORTANT

The mycelium living-repo scripts (`init_repo.py`, `generate_index.py`,
`recall_lessons.py`, `install_convention.py`, etc.) require **Python ≥ 3.11**
(they use `from datetime import UTC` and PEP 604 `X | None` syntax). The
default `python3` on this cluster is 3.9 and will raise `ImportError` /
`TypeError` on import.

**Always invoke mycelium scripts with `python3.12`**, e.g.:

```bash
MYC=/rhome/jstajich/.claude/plugins/marketplaces/mycelium/skills/core/scripts
python3.12 $MYC/recall_lessons.py --living-dir .living/ --tag <tag>
python3.12 $MYC/generate_index.py --living-dir .living/ --summary-heuristic
```

Note: the mycelium scripts live in the plugin marketplace dir above, **not**
in a `skills/core/scripts/` path inside this repo (the CLAUDE.md template
references the in-repo path, which does not exist here).

### Hook caveat

The SessionStart `mycelium-health.sh` hook calls `generate_index.py` via plain
`python3` (3.9), so automatic `.living/INDEX.md` regeneration **silently
no-ops** in this environment. Regenerate the index manually with `python3.12`
when learnings/decisions change. The hooks' lightweight `python3 -c
"import sys,json"` JSON helpers work fine on 3.9.

## Pipelines & dependencies

- **Nextflow** DSL2 workflows under `nextflow/` (ANI comparison, annotation,
  rnaseq, etc.); SLURM execution via layered `nextflow.config` +
  `conf/profile_*.config`.
- **funannotate**, **EarlGrey**, **Trinity**, **PHYling**, **Kraken2/Bracken**
  and related tools are provisioned per-process via modules / pixi / Singularity.
- See `HOWTO_singularity.md`, `HOWTO_k8s.md`, and the various `do_*.sh` /
  `run_*.sh` launchers for specifics.

## System Dependencies

- SLURM (sbatch/srun), environment modules, Singularity/Apptainer.

## Singularity images

- **GeneMark** (`genemark_sif`, used by `GENEMARK_RUN`): defaults to a local
  conversion of the public `docker://teambraker/braker3:v3.1.1` image
  (`/bigdata/stajichlab/shared/singularity_cache/braker3-v3.1.1.sif`), built
  by the [Gaius-Augustus/BRAKER](https://github.com/Gaius-Augustus/BRAKER)
  team. Confirmed 2026-08-18 to bundle GeneMark-ES/ET 4.72 (matching this
  project's earlier host module version) plus AUGUSTUS's `bam2hints` and
  `join_mult_hints.pl`, all already on `PATH` inside the image. GeneMark's
  license still requires a user-obtained key (`~/.gm_key`) regardless of
  which image runs it -- the image being public only means the GeneMark
  *binary* ships in it, not that a license is granted. Requires
  `--bind /opt/linux:/opt/linux` at runtime to resolve this host's
  `~/.gm_key` symlink chain to the actual key file (see
  `GENEMARK_RUN/main.nf`).
  - A privately-built, GeneMark-only alternative also exists in this cache
    (`genemark-4.72_lic.sif`), built via
    [hyphaltip/genemark-container](https://github.com/hyphaltip/genemark-container)
    and never pushed to a public registry (kept private since GeneMark's
    license forbids public redistribution of a *purpose-built* GeneMark
    image). Useful if a smaller, GeneMark-only image is preferred over
    braker3's full ~3.2 GB toolchain.

- **Container cache location** is `params.singularity_cache`, which `nextflow.config`
  resolves from the first set of `NXF_APPTAINER_CACHEDIR` / `NXF_SINGULARITY_CACHEDIR`
  / `APPTAINER_CACHEDIR` / `SINGULARITY_CACHEDIR` (no hardcoded default — the shared
  dir is `/bigdata/stajichlab/shared/singularity_cache`). Export one before an
  interactive `nextflow run`; the `run_*.sh` launchers do this for you. BFD and
  funannotate fail loudly at workflow start if it is unset. Warm/pre-mirror the cache
  with `nextflow/bin/setup_singularity_cache.sh` (`--dry-run` first). Every `*_sif`
  path param derives from it, so nothing else has a hardcoded cache path.
- **Apptainer/singularity** must be loadable in the *launching* shell for
  `apptainer.enabled = true` runs (driver-process needs it on PATH; see
  `.living/learnings.md`).
