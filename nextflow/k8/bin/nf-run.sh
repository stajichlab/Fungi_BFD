#!/usr/bin/env bash
# k8/bin/nf-run.sh — runs INSIDE an nf-job.sh Job: `nextflow "$@"` with the
# resume cache (.nextflow: LevelDB task cache + run history) on pod-local disk,
# snapshotted to the PVC atomically. nf-job.sh ships it in the Job's ConfigMap.
#
# Why: a node lost mid-write leaves the cache on CephFS corrupt, and every
# retry then fails with "Can't open cache DB ... Corruption" (seen with
# nf_funannotate1 on NRP, 2026-09-28). Nextflow's LevelDB (iq80) memory-maps
# its log and MANIFEST into 1 MB zero-filled files; the dead node's mapped
# writes never reach CephFS, and they are invisible to a file copy too.
#
# How: NXF_CACHE_DIR is an emptyDir (/nxf-local) and leveldb.mmap=false.
# Every CACHE_SYNC_SECONDS the cache is copied (retrying until no file changed
# during the copy) to $RUN_DIR/.nextflow-snapshots/.tmp-<ts>, fsynced, then
# renamed to snap-<ts>; the rename is atomic, so a snap-* is always complete.
# The newest CACHE_SNAPSHOTS_KEEP are kept, one more is taken when Nextflow
# exits, and the newest is restored on start (a run dir that still has
# .nextflow from before this change is migrated once). If the newest snapshot
# will not open, delete it and relaunch to fall back to the previous one.
#
# Env: RUN_DIR (launch dir, on the PVC), NF_LOG (log file to append to),
#      CACHE_SYNC_SECONDS (300), CACHE_SNAPSHOTS_KEEP (3).
set -uo pipefail

LOCAL=/nxf-local
SNAPS="$RUN_DIR/.nextflow-snapshots"
INTERVAL="${CACHE_SYNC_SECONDS:-300}"
KEEP="${CACHE_SNAPSHOTS_KEEP:-3}"
NF_LOG="${NF_LOG:-$RUN_DIR/nf_run.log}"

log() { echo "[nf-run $(date -u +%H:%M:%S)] $*" | tee -a "$NF_LOG"; }

cd "$RUN_DIR" || exit 1
mkdir -p "$SNAPS"
rm -rf "$SNAPS"/.tmp-*

latest=$(ls -1d "$SNAPS"/snap-* 2>/dev/null | sort | tail -1)
if [ -n "$latest" ]; then
    log "restoring cache from $latest"
    cp -a "$latest/.nextflow" "$LOCAL/"
elif [ -d "$RUN_DIR/.nextflow" ]; then
    log "migrating $RUN_DIR/.nextflow to local disk"
    cp -a "$RUN_DIR/.nextflow" "$LOCAL/"
    mv "$RUN_DIR/.nextflow" "$RUN_DIR/.nextflow.migrated-$(date -u +%Y%m%dT%H%M%S)"
fi
mkdir -p "$LOCAL/.nextflow"
export NXF_CACHE_DIR="$LOCAL/.nextflow"
export NXF_OPTS="${NXF_OPTS:-} -Dleveldb.mmap=false"

state() { ls -lAR --time-style=full-iso "$LOCAL/.nextflow" 2>/dev/null | md5sum; }

snapshot() {
    local ts tmp before i
    ts=$(date -u +%Y%m%dT%H%M%S)
    tmp="$SNAPS/.tmp-$ts"
    for i in 1 2 3 4 5 6 7 8 9 10; do
        rm -rf "$LOCAL/.snap"
        before=$(state)
        cp -a "$LOCAL/.nextflow" "$LOCAL/.snap" || return 1
        [ "$before" = "$(state)" ] && break
        [ "$i" = 10 ] && { log "cache kept changing; snapshot skipped"; return 1; }
        sleep 1
    done
    mkdir -p "$tmp"
    cp -a "$LOCAL/.snap" "$tmp/.nextflow" && sync -f "$tmp" || { rm -rf "$tmp"; return 1; }
    mv "$tmp" "$SNAPS/snap-$ts" || return 1
    ls -1d "$SNAPS"/snap-* | sort | head -n -"$KEEP" | while read -r old; do rm -rf "$old"; done
}

nextflow "$@" > >(tee -a "$NF_LOG") 2>&1 &
NF_PID=$!

# `kubectl delete job` sends SIGTERM here: pass it on so Nextflow cleans up
# its task pods, then fall through to the final snapshot.
trap 'log "SIGTERM: stopping Nextflow"; kill -TERM "$NF_PID" 2>/dev/null' TERM

while kill -0 "$NF_PID" 2>/dev/null; do
    sleep "$INTERVAL" & wait $! 2>/dev/null
    kill -0 "$NF_PID" 2>/dev/null && snapshot
done
wait "$NF_PID"; rc=$?
while kill -0 "$NF_PID" 2>/dev/null; do wait "$NF_PID"; rc=$?; done
snapshot && log "final cache snapshot saved" || log "WARNING: final cache snapshot failed"
exit "$rc"
