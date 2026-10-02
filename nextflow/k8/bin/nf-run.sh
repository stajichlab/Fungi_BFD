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
# .nextflow from before this change is migrated once). A restored snapshot that
# will not open is renamed bad-* and the Job's retry falls back to the previous
# one. A new pod waits for a previous pod still stopping on the same run dir
# (heartbeat file), so `kubectl delete` + relaunch can't copy a live cache.
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

# One Nextflow per run dir. After `kubectl delete job`, the old pod can still be
# shutting down (and writing its final snapshot) when the relaunched Job's pod
# starts, so wait for its heartbeat to go stale. (The image has no flock.)
HB="$RUN_DIR/.nf-run.heartbeat"
heartbeat() { echo "$HOSTNAME $(date +%s)" > "$HB.tmp" && mv "$HB.tmp" "$HB"; }
for i in $(seq 60); do
    [ -f "$HB" ] && read -r hb_host hb_time < "$HB" || break
    [ "$hb_host" = "$HOSTNAME" ] && break
    [ $(( $(date +%s) - ${hb_time:-0} )) -gt 90 ] && break
    [ "$i" = 1 ] && log "waiting for $hb_host to finish with $RUN_DIR"
    [ "$i" = 60 ] && { log "ERROR: $hb_host still active after 30 min"; exit 1; }
    sleep 30
done
heartbeat
rm -rf "$SNAPS"/.tmp-*

latest=$(ls -1d "$SNAPS"/snap-* 2>/dev/null | sort | tail -1)
if [ -n "$latest" ]; then
    log "restoring cache from $latest"
    cp -a "$latest/.nextflow" "$LOCAL/"
elif [ -d "$RUN_DIR/.nextflow" ]; then
    # A Nextflow from before this change (no heartbeat) may still be stopping:
    # wait until its .nextflow.log has been quiet for 2 min.
    while [ $(( $(date +%s) - $(stat -c %Y .nextflow.log 2>/dev/null || echo 0) )) -lt 120 ]; do
        log "waiting for .nextflow.log to go quiet (previous Nextflow stopping?)"
        sleep 30
    done
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

last=$(date +%s)
while kill -0 "$NF_PID" 2>/dev/null; do
    heartbeat
    sleep 30 & wait $! 2>/dev/null
    if kill -0 "$NF_PID" 2>/dev/null && [ $(( $(date +%s) - last )) -ge "$INTERVAL" ]; then
        snapshot; last=$(date +%s)
    fi
done
wait "$NF_PID"; rc=$?
while kill -0 "$NF_PID" 2>/dev/null; do wait "$NF_PID"; rc=$?; done

if [ "$rc" -ne 0 ] && [ -n "$latest" ] && grep -q "Can't open cache DB" "$NF_LOG"; then
    # The restored snapshot is unreadable: set it aside (never snapshot it
    # again) so the Job's retry falls back to the previous one.
    mv "$latest" "$SNAPS/bad-$(basename "$latest")"
    log "ERROR: $(basename "$latest") would not open; moved to bad-*. Retry uses the previous snapshot."
else
    snapshot && log "final cache snapshot saved" || log "WARNING: final cache snapshot failed"
fi
rm -f "$HB"
exit "$rc"
