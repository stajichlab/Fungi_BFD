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

# Timing is logged to stderr (captured in .command.log by Nextflow) and
# appended as a TSV line to fcs_gx_shm_timing.tsv in the task work dir, so
# per-task sync cost can be aggregated/reported across the pipeline run.
FCS_GX_TIMING_LOG="${FCS_GX_TIMING_LOG:-${PWD}/fcs_gx_shm_timing.tsv}"
if [ ! -f "${FCS_GX_TIMING_LOG}" ]; then
    printf 'timestamp\thost\tstatus\telapsed_sec\tbytes\tmbps\n' > "${FCS_GX_TIMING_LOG}"
fi

if [ -f "${FCS_GX_SHM_DIR}/all.gxi" ]; then
    echo "[setup_fcs_shm] FCS-GX db already present in ${FCS_GX_SHM_DIR}; reusing" >&2
    printf '%s\t%s\treused\t0\t-\t-\n' "$(date -Is)" "$(hostname -s)" >> "${FCS_GX_TIMING_LOG}"
else
    if [ ! -d "${FCS_GX_DB_SRC}" ]; then
        echo "[setup_fcs_shm] ERROR: FCS_GX_DB_SRC not found: ${FCS_GX_DB_SRC}" >&2
        echo "[setup_fcs_shm] Set FCS_GX_DB_SRC to your FCS-GX gxdb directory." >&2
        return 1 2>/dev/null || exit 1
    fi
    echo "[setup_fcs_shm] Syncing FCS-GX db ${FCS_GX_DB_SRC} -> ${FCS_GX_SHM_DIR}" >&2
    mkdir -p "${FCS_GX_SHM_DIR}"
    # gxdb is a flat directory dominated by two huge files (all.gxi ~300GB,
    # all.gxs ~165GB); a single rsync stream can't saturate the IB fabric /
    # parallel filesystem, so fan out one rsync per file (-P4 matches typical
    # node core count). --whole-file skips delta-transfer checksumming, which
    # is pure overhead here since the destination is always freshly emptied.
    fcs_gx_sync_start=$(date +%s)
    find "${FCS_GX_DB_SRC}/" -maxdepth 1 -type f -printf '%f\0' | \
        xargs -0 -P 4 -I{} rsync -a -W --inplace "${FCS_GX_DB_SRC}/{}" "${FCS_GX_SHM_DIR}/{}"
    fcs_gx_sync_status=$?
    fcs_gx_sync_end=$(date +%s)
    fcs_gx_sync_elapsed=$(( fcs_gx_sync_end - fcs_gx_sync_start ))
    fcs_gx_sync_bytes=$(du -sb "${FCS_GX_SHM_DIR}" 2>/dev/null | cut -f1)
    fcs_gx_sync_mbps="-"
    if [ -n "${fcs_gx_sync_bytes}" ] && [ "${fcs_gx_sync_elapsed}" -gt 0 ]; then
        fcs_gx_sync_mbps=$(awk -v b="${fcs_gx_sync_bytes}" -v s="${fcs_gx_sync_elapsed}" \
            'BEGIN { printf "%.1f", (b / 1048576) / s }')
    fi
    if [ "${fcs_gx_sync_status}" -eq 0 ]; then
        echo "[setup_fcs_shm] Sync completed in ${fcs_gx_sync_elapsed}s (${fcs_gx_sync_bytes:-?} bytes, ${fcs_gx_sync_mbps} MB/s)" >&2
        printf '%s\t%s\tsynced\t%s\t%s\t%s\n' "$(date -Is)" "$(hostname -s)" \
            "${fcs_gx_sync_elapsed}" "${fcs_gx_sync_bytes:--}" "${fcs_gx_sync_mbps}" >> "${FCS_GX_TIMING_LOG}"
    else
        echo "[setup_fcs_shm] ERROR: rsync failed after ${fcs_gx_sync_elapsed}s (exit ${fcs_gx_sync_status})" >&2
        printf '%s\t%s\tfailed\t%s\t%s\t-\n' "$(date -Is)" "$(hostname -s)" \
            "${fcs_gx_sync_elapsed}" "${fcs_gx_sync_bytes:--}" >> "${FCS_GX_TIMING_LOG}"
        return 1 2>/dev/null || exit 1
    fi
fi

# Register cleanup so the (large) RAM copy is removed when this shell exits,
# unless FCS_GX_KEEP_SHM=1 (useful when many tasks share a node sequentially).
if [ "${FCS_GX_KEEP_SHM:-0}" != "1" ]; then
    trap 'rm -rf "'"${FCS_GX_SHM_DIR}"'" 2>/dev/null || true' EXIT
fi
