#!/usr/bin/env bash
# k8/bin/ani-suite.sh — launch ani-run.sh for a whole list of taxa at once.
#
# Each launch returns almost immediately (ani-run.sh backgrounds the actual
# work on the cluster), so this effectively starts them all concurrently —
# each with its own nextflow head process and its own queueSize=30 pod cap.
# Watch your namespace's real capacity if you launch a large batch (see
# k8/README_compare_ani.md's Tuning section); there's no hard quota to stop
# you from over-requesting, only cluster-wide fair-share.
#
# Usage:
#   ani-suite.sh --compare SPECIES --taxa GENUS:Yarrowia,GENUS:Aspergillus
#   ani-suite.sh --compare SPECIES --taxa-file genera.txt   # one RANK:VALUE per line, # comments ok
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

COMPARE="GENUS"
METHOD="skani"
TAXA=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --compare) COMPARE="$2"; shift 2 ;;
    --method) METHOD="$2"; shift 2 ;;
    --taxa) IFS=',' read -r -a TAXA <<< "$2"; shift 2 ;;
    --taxa-file)
      while IFS= read -r line; do
        line="${line%%#*}"
        line="$(echo "$line" | xargs || true)"
        [[ -n "$line" ]] && TAXA+=("$line")
      done < "$2"
      shift 2 ;;
    *) echo "Unknown arg: $1" >&2; exit 1 ;;
  esac
done

if [[ ${#TAXA[@]} -eq 0 ]]; then
  echo "Usage: $0 --compare LEVEL (--taxa RANK:VAL,RANK:VAL,... | --taxa-file file)" >&2
  exit 1
fi

echo "Launching ${#TAXA[@]} run(s), compare=${COMPARE} method=${METHOD}:"
for taxon in "${TAXA[@]}"; do
  "${HERE}/ani-run.sh" --taxon "$taxon" --compare "$COMPARE" --method "$METHOD"
done

echo
echo "All launched. Check progress with: k8/bin/ani-status.sh"
echo "Once a taxon's compute phase looks done, gather it with: k8/bin/ani-gather.sh --taxon RANK:VALUE --compare ${COMPARE}"
