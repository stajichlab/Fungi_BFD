#!/usr/bin/env bash
# setup_singularity_cache.sh — warm the shared Apptainer/Singularity cache that
# the BFD Nextflow pipelines reference, so SLURM tasks never race to pull/convert
# OCI images at run time.
#
# It manages EXACTLY the images the pipeline configs actively reference:
#
#   Section A — docker:// container URIs kept in the configs (see conf/*.config).
#       Pre-pulls each one into the cache. The pull writes a .sif AND warms
#       Apptainer's OCI layer-blob cache in the same dir, so Nextflow's later
#       on-the-fly docker:// -> .sif conversion for a job reuses the cached
#       blobs instead of re-downloading (the conversion race that corrupted the
#       shared blob cache under concurrent SLURM tasks).
#
#   Section B — pre-built .sif files pinned by absolute path in the configs.
#       Ensures the exact named file exists in the cache:
#         * public images (AAFTF, funannotate, braker3): copied from
#           --source-cache if present, else pulled (unless --skip-build);
#           they must exist under their pinned name because the configs point
#           at the file directly.
#         * mariadb-10.3.9.sif: correctness fix — must be BUILT from
#           nextflow/docs/mariadb-10.3.9.def (a bare `apptainer pull
#           docker://mariadb:10.3.9` starts an empty instance whose mysqld never
#           launches; see conf/profile_funannotate.config). Copied from
#           --source-cache if present, else built from the .def (unless
#           --skip-build).
#         * license/private images (DeepTMHMM-1.0.sif academic license,
#           antismash-standalone-8.0.4.sif commercial standalone): copy-only.
#           Never pulled/built. If absent from the cache AND --source-cache,
#           warn + skip (the referencing step will fail loudly at run time).
#
# CACHE DIR — env-only, NO baked-in HPCC default (project convention, 2026-08-25).
# Resolved, first match wins:
#   1. --cache-dir DIR
#   2. $NXF_APPTAINER_CACHEDIR
#   3. $NXF_SINGULARITY_CACHEDIR
#   4. $APPTAINER_CACHEDIR
#   5. $SINGULARITY_CACHEDIR
# If none is set the script fails with exit 2. Same resolution feeds the
# pipeline configs' params.singularity_cache at run time, so warm exactly the
# dir that value points at.
#
# RUNTIME — prefers `apptainer` on PATH, falls back to `singularity`, then tries
# `module load apptainer` and `module load singularity`. Runs on the login node
# or an sbatch job; no sudo/root needed (pull needs none; build of the mariadb
# .def may need network + user namespaces — see the file help).
#
# Usage:
#   setup_singularity_cache.sh [--cache-dir DIR] [--source-cache DIR]
#                              [--repo-root DIR] [--skip-build]
#                              [--dry-run] [--check] [--verbose] [-h]
#
#   --cache-dir DIR     Cache to populate (default: env chain above).
#   --source-cache DIR  Existing populated cache to copy license/private .sif
#                       (and pinned public .sif) from. Optional.
#   --repo-root DIR     Repo checkout root (for the mariadb .def discovery).
#                       Defaults to $PWD if nextflow/main.nf exists there, else
#                       the nextflow/ dir sibling of this script.
#   --skip-build        Do not pull named public .sif or build mariadb from its
#                       .def; only copy from --source-cache / warn.
#   --dry-run           Report what would be done without changing anything.
#   --check             Verify presence of every referenced image; exit 0 if all
#                       present, 1 if anything missing. Never mutates.
#   --verbose           Log the per-image action taken.

set -euo pipefail

# ├─ parse args ────────────────────────────────────────────────────────────────
CACHE_DIR=""
SOURCE_CACHE=""
REPO_ROOT=""
SKIP_BUILD=0
DRY_RUN=0
CHECK_ONLY=0
VERBOSE=0

usage() {
    cat <<'EOF'
Usage:
  setup_singularity_cache.sh [--cache-dir DIR] [--source-cache DIR]
                             [--repo-root DIR] [--skip-build]
                             [--dry-run] [--check] [--verbose] [-h]

  --cache-dir DIR     Cache to populate (default: env chain below).
  --source-cache DIR  Existing populated cache to copy license/private .sif
                      (and pinned public .sif) from. Optional.
  --repo-root DIR     Repo checkout root (for the mariadb .def discovery).
                      Defaults to $PWD if nextflow/main.nf exists there, else
                      the nextflow/ dir sibling of this script.
  --skip-build        Do not pull named public .sif or build mariadb from its
                      .def; only copy from --source-cache / warn.
  --dry-run           Report what would be done without changing anything.
  --check             Verify presence of every referenced image; exit 0 if all
                      present, 1 if anything missing. Never mutates.
  --verbose           Log the per-image action taken.

Cache dir is env-only, NO baked-in default (first match wins):
  --cache-dir DIR, then $NXF_APPTAINER_CACHEDIR, $NXF_SINGULARITY_CACHEDIR,
  $APPTAINER_CACHEDIR, $SINGULARITY_CACHEDIR. If none set: exit 2.
EOF
    exit "${1:-0}"
}

while [ $# -gt 0 ]; do
    case "$1" in
        --cache-dir)   CACHE_DIR="$2";                 shift 2 ;;
        --source-cache) SOURCE_CACHE="$2";             shift 2 ;;
        --repo-root)   REPO_ROOT="$2";                 shift 2 ;;
        --skip-build)  SKIP_BUILD=1;                   shift   ;;
        --dry-run)     DRY_RUN=1;                      shift   ;;
        --check)       CHECK_ONLY=1;                   shift   ;;
        --verbose)     VERBOSE=1;                      shift   ;;
        -h|--help)     usage 0 ;;
        *) echo "error: unknown option: $1" >&2; usage 2 ;;
    esac
done

# ├─ resolve cache dir (env-only, no default) ──────────────────────────────────
if [ -z "${CACHE_DIR}" ]; then
    CACHE_DIR="${NXF_APPTAINER_CACHEDIR:-${NXF_SINGULARITY_CACHEDIR:-${APPTAINER_CACHEDIR:-${SINGULARITY_CACHEDIR:-}}}}"
fi
if [ -z "${CACHE_DIR}" ]; then
    echo "ERROR: no cache dir. Export one of NXF_APPTAINER_CACHEDIR / NXF_SINGULARITY_CACHEDIR / APPTAINER_CACHEDIR / SINGULARITY_CACHEDIR, or pass --cache-dir DIR." >&2
    exit 2
fi

# ├─ resolve repo root (for the mariadb .def) ─────────────────────────────────
if [ -z "${REPO_ROOT}" ]; then
    if [ -f "$(pwd)/nextflow/main.nf" ]; then
        REPO_ROOT="$(pwd)"
    elif [ -f "$(pwd)/main.nf" ]; then
        REPO_ROOT="$(pwd)"
    elif [ -n "${SLURM_SUBMIT_DIR:-}" ] && [ -f "${SLURM_SUBMIT_DIR}/nextflow/main.nf" ]; then
        REPO_ROOT="${SLURM_SUBMIT_DIR}"
    elif [ -n "${BASH_SOURCE[0]:-}" ] && [ -f "$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/nextflow/main.nf" ]; then
        REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
    else
        echo "WARNING: could not auto-detect repo root (nextflow/main.nf not found in \$PWD or near this script) — mariadb .def build won't be available. Pass --repo-root DIR." >&2
        REPO_ROOT=""
    fi
fi
MARIADB_DEF=""
if [ -n "${REPO_ROOT}" ]; then
    for cand in "${REPO_ROOT}/nextflow/docs/mariadb-10.3.9.def" "${REPO_ROOT}/docs/mariadb-10.3.9.def"; do
        [ -f "${cand}" ] && { MARIADB_DEF="${cand}"; break; }
    done
fi

# ├─ resolve apptainer/singularity runtime ─────────────────────────────────────
log2() { printf '%s\n' "$*" >&2; }
info() { log2 "[setup_singularity_cache] $*"; }

resolve_runtime() {
    local rt=""
    if command -v apptainer >/dev/null 2>&1; then rt="apptainer"; fi
    if [ -z "${rt}" ] && command -v singularity >/dev/null 2>&1; then rt="singularity"; fi
    if [ -z "${rt}" ]; then
        if command -v module >/dev/null 2>&1 || { command -v source >/dev/null 2>&1 && [ -f /etc/profile.d/modules.sh ]; }; then
            source /etc/profile.d/modules.sh 2>/dev/null || true
        fi
        if [ -z "${rt}" ] && command -v module >/dev/null 2>&1; then
            module load apptainer    2>/dev/null || true
            command -v apptainer >/dev/null 2>&1 && rt="apptainer"
        fi
        if [ -z "${rt}" ] && command -v module >/dev/null 2>&1; then
            module load singularity  2>/dev/null || true
            command -v singularity >/dev/null 2>&1 && rt="singularity"
        fi
    fi
    if [ -z "${rt}" ]; then
        echo "ERROR: no apptainer or singularity found on PATH and module load failed. module avail shows apptainer/1.4.5, singularity-ce. Load manually and retry." >&2
        exit 3
    fi
    RT="${rt}"
    log2 "[setup_singularity_cache] runtime: ${RT}"
}
RT=""
resolve_runtime
# read() on this script (usage) is the only BASH_SOURCE use; actual work below is SLURM-safe.

# ├─ inventory ─────────────────────────────────────────────────────────────────
# Section A: docker:// container directives kept in the configs.
#   format: "<config-ref-name>\t<docker:// URI>"
DOCKER_REFS=$(cat <<'EOF'
pfam_sif	docker://quay.io/biocontainers/hmmer:3.4--hdbdd923_1
diamond_sif	docker://quay.io/biocontainers/diamond:2.2.5--he361c42_0
blastp_sif	docker://quay.io/biocontainers/blast:2.16.0--h66d330f_5
eggnog_sif	docker://quay.io/biocontainers/eggnog-mapper:2.1.15--pyhdfd78af_0
trnascan_sif	docker://quay.io/biocontainers/trnascan-se:2.0.13--pl5321hab16a5f_0
infernal_sif	docker://rnacentral/infernal:latest
dbcan_cgc_sif	docker://quay.io/biocontainers/dbcan:5.2.9--pyhdfd78af_0
skani_container	docker://quay.io/biocontainers/skani:0.3.2--h79ce301_0
mash_container	docker://quay.io/biocontainers/mash:2.3--hb105d93_10
sourmash_container	docker://quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0
fastani_container	docker://quay.io/biocontainers/fastani:1.34--hb66fcc3_5
report_container	docker://python:3.12
antismash_sif	docker://quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1
wgd_container	docker://quay.io/biocontainers/wgd:2.0.38--pyhdfd78af_0
EOF
)
# wgd_container: NOT yet referenced by any config/module (no funannotate/BFD/ANI
# step calls wgd today) -- pre-staged here ahead of that wiring, 2026-08-26.
# quay.io/biocontainers/wgd has no literal "latest" tag (confirmed via the
# quay.io API: `?specificTag=latest` returns zero tags) -- `docker pull
# quay.io/biocontainers/wgd:latest` would fail outright. 2.0.38--pyhdfd78af_0
# is the actual newest version (built 2024-05-24; no release since), pinned
# here to match every other entry in this table instead of an unpinned tag.

# Section B: pre-built .sif (or sandbox dir) pinned by absolute path in the configs.
#   name|kind|config-param|source docker URI|source .def
#   kind: public  = pullable docker image (build unless --skip-build)
#         build   = must build from a .def (mariadb correctness fix)
#         license = copy-only from --source-cache (warn + skip if absent)
#         sandbox = extracted directory tree, not a .sif -- built from the
#                   docker URI (reusing an already-cached Section-A .sif for
#                   that same URI if present, else pulling fresh). Apptainer
#                   bind-mounts a sandbox directly with no FUSE/squashfuse
#                   layer, unlike a .sif -- use this kind for containers run
#                   at very high task concurrency, where many simultaneous
#                   squashfuse mounts of the same .sif can time out under
#                   contention (see skani_0.3.2--h79ce301_0.sandbox below,
#                   added 2026-08-26 after squashfuse_ll timeouts failed
#                   ~23% of SKANI_COMPARE tasks on a real compare_ani run;
#                   the .sif itself was already fully cached, so this was
#                   never a "missing image" problem -- see
#                   nextflow/conf/profile_ANI.config's skani_container
#                   comment for the full incident).
PINNED_SIFS="\
AAFTF-latest.sif|public|aaftf_sif|docker://ghcr.io/stajichlab/aaftf:latest|
funannotate-1.9.0-beta.8.sif|public|funannotate_sif|docker://ghcr.io/nextgenusfs/funannotate:1.9.0-beta.8|
braker3-v3.1.1.sif|public|genemark_sif|docker://teambraker/braker3:v3.1.1|
mariadb-10.3.9.sif|build|mariadb_sif||${MARIADB_DEF}
DeepTMHMM-1.0.sif|license|deeptmhmm_sif||
antismash-standalone-8.0.4.sif|license|antismash_standalone_sif||
skani_0.3.2--h79ce301_0.sandbox|sandbox|skani_container|docker://quay.io/biocontainers/skani:0.3.2--h79ce301_0|
mash_2.3--hb105d93_10.sandbox|sandbox|mash_container|docker://quay.io/biocontainers/mash:2.3--hb105d93_10|
sourmash_4.9.4--hdfd78af_0.sandbox|sandbox|sourmash_container|docker://quay.io/biocontainers/sourmash:4.9.4--hdfd78af_0|
"

# known-but-not-actively-wired (comment-referenced fallbacks / other projects);
# NOT managed here, listed only so a reader doesn't mistake them for covered.
NOT_MANAGED="genemark-4.72_lic.sif signalp6-fast.sif signalp6-fast-gpu-models"

# ├─ helpers ───────────────────────────────────────────────────────────────────
# Apptainer's default pull output name: last URI path component with ':' -> '_'.
sif_name_from_uri() {
    local uri="$1"
    local name
    name="${uri#docker://}"                        # strip scheme
    name="${name##*/}"                             # last path component
    name="${name%%@*}"                             # drop digest if present
    printf '%s.sif' "${name//:/_}"
}

cache_has() { [ -s "${CACHE_DIR}/$1" ]; }
# Sandbox images are directories, not regular files -- `-s` always fails on a
# directory, so this needs its own existence check (non-empty dir).
cache_has_dir() { [ -d "${CACHE_DIR}/$1" ] && [ -n "$(ls -A "${CACHE_DIR}/$1" 2>/dev/null)" ]; }

do_cp_from_source() {
    local dst="$1" src="$2"
    if [ -z "${SOURCE_CACHE}" ]; then return 1; fi
    if [ ! -s "${SOURCE_CACHE}/${src}" ]; then
        info "source-cache lacks ${src}; nothing to copy"
        return 1
    fi
    if [ "${CHECK_ONLY}" = 1 ] || [ "${DRY_RUN}" = 1 ]; then
        log2 "[would copy]  ${SOURCE_CACHE}/${src} -> ${CACHE_DIR}/${dst}"
        return 0
    fi
    mkdir -p "${CACHE_DIR}"
    cp -p "${SOURCE_CACHE}/${src}" "${CACHE_DIR}/${dst}"
    log2 "[copied]      ${SOURCE_CACHE}/${src} -> ${CACHE_DIR}/${dst}"
    return 0
}

do_pull() {
    local dst="$1" uri="$2" extra=""
    if [ "${CHECK_ONLY}" = 1 ] || [ "${DRY_RUN}" = 1 ]; then
        log2 "[would pull]  ${uri}"
        return 0
    fi
    mkdir -p "${CACHE_DIR}"
    # Pin the .sif name with the positional <target> arg before the URI
    # (portable across apptainer AND singularity; --name is not reliably
    # supported). With no target, --dir funnels the derived-name output into
    # CACHE_DIR.
    if [ -n "${dst}" ]; then
        "${RT}" pull "${CACHE_DIR}/${dst}" "${uri}"
    else
        "${RT}" pull --dir "${CACHE_DIR}" "${uri}"
    fi
}

do_build_def() {
    local dst="$1" def="$2"
    if [ "${CHECK_ONLY}" = 1 ] || [ "${DRY_RUN}" = 1 ]; then
        log2 "[would build] ${def} -> ${CACHE_DIR}/${dst}"
        return 0
    fi
    mkdir -p "${CACHE_DIR}"
    "${RT}" build "${CACHE_DIR}/${dst}" "${def}"
}

do_build_sandbox() {
    local dst="$1" uri="$2"
    if [ "${CHECK_ONLY}" = 1 ] || [ "${DRY_RUN}" = 1 ]; then
        log2 "[would build]  ${dst}/ (sandbox) from ${uri}"
        return 0
    fi
    mkdir -p "${CACHE_DIR}"
    # Prefer building from an already-cached local .sif for this exact URI
    # (fast, local, no network) if Section A (or a prior run of this script)
    # already pulled it under its default derived name; else build the
    # sandbox straight from the docker:// URI (apptainer supports that
    # directly -- slower, needs network, but works standalone too).
    local local_sif
    local_sif="$(sif_name_from_uri "${uri}")"
    if [ -s "${CACHE_DIR}/${local_sif}" ]; then
        log2 "[build]        ${dst}/ (sandbox) from local ${local_sif}"
        "${RT}" build --sandbox "${CACHE_DIR}/${dst}" "${CACHE_DIR}/${local_sif}"
    else
        log2 "[build]        ${dst}/ (sandbox) from ${uri} (no local .sif cached -- pulling)"
        "${RT}" build --sandbox "${CACHE_DIR}/${dst}" "${uri}"
    fi
}

MISSING=0

# ├─ Section A: prefetch docker:// refs ────────────────────────────────────────
if [ "${VERBOSE}" = 1 ]; then info "warm docker:// refs referenced by the configs"; fi
while IFS=$'\t' read -r ref uri; do
    [ -z "${ref}" ] && continue
    dst="$(sif_name_from_uri "${uri}")"
    if cache_has "${dst}"; then
        log2 "[present]      ${dst} (${ref})"
        continue
    fi
    if [ "${CHECK_ONLY}" = 1 ]; then
        log2 "[MISSING]     ${dst} (${ref})"
        MISSING=1
        continue
    fi
    do_pull "" "${uri}"
done <<< "${DOCKER_REFS}"

# ├─ Section B: pinned .sif ────────────────────────────────────────────────────
if [ "${VERBOSE}" = 1 ]; then info "ensure pinned .sif files referenced by absolute path"; fi
IFS=$'\n'
for row in ${PINNED_SIFS}; do
    [ -z "${row}" ] && continue
    name="${row%%|*}"; rest="${row#*|}"; kind="${rest%%|*}"; rest="${rest#*|}"
    param="${rest%%|*}"; rest="${rest#*|}"; srcuri="${rest%%|*}"; rest="${rest#*|}"
    def="${rest}"

    if [ "${kind}" = "sandbox" ]; then
        present=false; cache_has_dir "${name}" && present=true
    else
        present=false; cache_has "${name}" && present=true
    fi
    if [ "${present}" = true ]; then
        log2 "[present]      ${name} (${param})"
        continue
    fi
    if [ "${CHECK_ONLY}" = 1 ]; then
        log2 "[MISSING]     ${name} (${param})"
        MISSING=1
        continue
    fi

    case "${kind}" in
        license)
            if do_cp_from_source "${name}" "${name}"; then
                log2 "[copied]      ${name} (${param}) from --source-cache"
            else
                log2 "[WARN]        ${name} (${param}) absent and no --source-cache copy — license image cannot be pulled/built. The referencing step will fail loudly at run time."
                MISSING=1
            fi
            ;;
        public)
            if do_cp_from_source "${name}" "${name}"; then
                log2 "[copied]      ${name} (${param}) from --source-cache"
            elif [ "${SKIP_BUILD}" = 1 ]; then
                log2 "[WARN]        ${name} (${param}) absent; --skip-build set, not pulling. Run without --skip-build or point --source-cache at a populated cache."
                MISSING=1
            else
                log2 "[pull]        ${name} (${param}) from ${srcuri}"
                do_pull "${name}" "${srcuri}"
            fi
            ;;
        build)
            if do_cp_from_source "${name}" "${name}"; then
                log2 "[copied]      ${name} (${param}) from --source-cache"
            elif [ -n "${def}" ] && [ "${SKIP_BUILD}" != 1 ]; then
                log2 "[build]       ${name} (${param}) from ${def}"
                do_build_def "${name}" "${def}"
            elif [ -z "${def}" ]; then
                log2 "[WARN]        ${name} (${param}) absent and mariadb .def not found (--repo-root). Cannot build the correctness-fixed image; a bare docker:// pull would launch an EMPTY mysqld instance."
                MISSING=1
            else
                log2 "[WARN]        ${name} (${param}) absent; --skip-build set. Run without --skip-build or point --source-cache at a populated cache."
                MISSING=1
            fi
            ;;
        sandbox)
            if [ "${SKIP_BUILD}" = 1 ]; then
                log2 "[WARN]        ${name} (${param}) absent; --skip-build set, not building. Run without --skip-build (a --source-cache copy is not supported for sandbox dirs -- cp -r manually if needed)."
                MISSING=1
            else
                do_build_sandbox "${name}" "${srcuri}"
            fi
            ;;
    esac
done
unset IFS

# ├─ summary ───────────────────────────────────────────────────────────────────
if [ "${CHECK_ONLY}" = 1 ]; then
    if [ "${MISSING}" = 1 ]; then
        log2 "[setup_singularity_cache] CHECK FAILURE: one or more referenced images missing from ${CACHE_DIR}"
        exit 1
    fi
    log2 "[setup_singularity_cache] CHECK OK: all referenced images present in ${CACHE_DIR}"
    exit 0
fi

log2 "[setup_singularity_cache] cache dir: ${CACHE_DIR}"
log2 "[setup_singularity_cache] not managed here (fallbacks/other projects): ${NOT_MANAGED}"
log2 "[setup_singularity_cache] reminder: runs must export NXF_APPTAINER_CACHEDIR/NXF_SINGULARITY_CACHEDIR + APPTAINER_CACHEDIR pointing at the SAME dir (the Nextflow cacheDir and the apptainer binary blob cache are separate knobs)."
log2 "[setup_singularity_cache] reminder: mariadb must NEVER come from bare 'apptainer pull docker://mariadb:10.3.9' — build from nextflow/docs/mariadb-10.3.9.def (correctness fix)."
if [ "${MISSING}" = 1 ] && [ "${SKIP_BUILD}" = 1 ]; then
    log2 "[setup_singularity_cache] one or more images were skipped (--skip-build); re-run without it or point --source-cache at a populated cache before launching a run."
elif [ "${MISSING}" = 1 ]; then
    log2 "[setup_singularity_cache] one or more images could not be provisioned (see WARN above); the referencing step will fail loudly at run time."
fi
exit 0
