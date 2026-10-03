#!/usr/bin/bash
# Find FUNANNOTATE_TRAIN/UPDATE tasks hung after starting MariaDB (Fungi_BFD_runs
# DECISIONS D138). Before Fungi_BFD 8e0b33c, a mariadbd that took longer than 30s to
# come up sent the task down a retry path that ran an unbounded `wait` on the live
# mariadbd; the task then idled until its time limit. A Nextflow run started before
# 8e0b33c keeps the old task code, so this finds the hung tasks.
#
# A task is reported when its SLURM job is RUNNING, the last line of .command.out is the
# "[INFO] Starting ... mariadbd" line, and .command.out is older than MINUTES (default 15).
# With --cancel the job is cancelled; with executor.queueGlobalStatus = true Nextflow
# marks the task failed and retries it.
#
# Usage: find_hung_mariadb_tasks.sh <run folder> [MINUTES] [--cancel]
RUN=${1:?run folder}
MIN=${2:-15}
CANCEL=${3:-}
now=$(date +%s)
squeue -u "$USER" -h -t R -o "%i|%Z|%j" | grep -E "nf-TRAIN|nf-UPDATE|TRAIN|UPDATE" | while IFS='|' read -r job wd name; do
    case "$wd" in "$RUN"/*) ;; *) continue ;; esac
    out="$wd/.command.out"
    [ -f "$out" ] || continue
    last=$(grep -v '^[[:space:]]*$' "$out" | tail -n 1)
    case "$last" in "[INFO] Starting "*mariadbd*) ;; *) continue ;; esac
    age=$(( (now - $(stat -c %Y "$out")) / 60 ))
    [ "$age" -ge "$MIN" ] || continue
    task=$(grep -m1 "^### name:" "$wd/.command.run" 2>/dev/null | sed "s/^### name: //")
    echo "$(date -Iseconds) HUNG job=$job idle=${age}min task=$task workdir=$wd"
    if [ "$CANCEL" = "--cancel" ]; then
        scancel "$job" && echo "$(date -Iseconds) cancelled $job"
    fi
done
