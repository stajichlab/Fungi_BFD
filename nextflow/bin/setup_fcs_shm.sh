# setup_fcs_shm.sh — stage the NCBI FCS-GX database into node-local /dev/shm.
#
# SOURCED (not executed) by the GENOME_CLEAN process before AAFTF fcs_gx_purge,
# which expects the GX database under  /dev/shm/gxdb/all{.gxi,.gxs,...}.
# Keeping the DB in /dev/shm (RAM) is required for acceptable FCS-GX performance.
#
#   source scripts/setup_fcs_shm.sh
#   AAFTF fcs_gx_purge --db /dev/shm/gxdb/all ...
#
# Configure the source location of the GX database via FCS_GX_DB_SRC (a directory
# containing all.gxi / all.gxs etc). The default below is a placeholder — set it
# to your site's FCS-GX database path (or export FCS_GX_DB_SRC before launching
# the pipeline / in the profile env).
#
# This script is idempotent: if /dev/shm/gxdb/all.gxi already exists it does
# nothing, so concurrent GENOME_CLEAN tasks on the same node share one copy.

: "${FCS_GX_DB_SRC:=/srv/projects/db/ncbi-fcs/0.5.4/gxdb}"
: "${FCS_GX_SHM_DIR:=/dev/shm/gxdb}"

if [ -f "${FCS_GX_SHM_DIR}/all.gxi" ]; then
    echo "[setup_fcs_shm] FCS-GX db already present in ${FCS_GX_SHM_DIR}; reusing" >&2
else
    if [ ! -d "${FCS_GX_DB_SRC}" ]; then
        echo "[setup_fcs_shm] ERROR: FCS_GX_DB_SRC not found: ${FCS_GX_DB_SRC}" >&2
        echo "[setup_fcs_shm] Set FCS_GX_DB_SRC to your FCS-GX gxdb directory." >&2
        return 1 2>/dev/null || exit 1
    fi
    echo "[setup_fcs_shm] Syncing FCS-GX db ${FCS_GX_DB_SRC} -> ${FCS_GX_SHM_DIR}" >&2
    mkdir -p "${FCS_GX_SHM_DIR}"
    # --inplace + large files: rsync keeps the RAM copy current without doubling space
    rsync -a --inplace "${FCS_GX_DB_SRC}/" "${FCS_GX_SHM_DIR}/"
fi

# Register cleanup so the (large) RAM copy is removed when this shell exits,
# unless FCS_GX_KEEP_SHM=1 (useful when many tasks share a node sequentially).
if [ "${FCS_GX_KEEP_SHM:-0}" != "1" ]; then
    trap 'rm -rf "'"${FCS_GX_SHM_DIR}"'" 2>/dev/null || true' EXIT
fi
