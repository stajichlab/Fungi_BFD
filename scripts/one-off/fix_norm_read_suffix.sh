#!/usr/bin/env bash
# Add /1 and /2 read-name suffixes to <species>_norm_R1/R2.fastq.gz files that lack them.
# Dry-run by default. Pass --apply to change files.
# Usage: fix_norm_read_suffix.sh [--apply] SPECIES [SPECIES ...]
# Originals are kept as <file>.pre_rename in the same directory until you delete them.
set -euo pipefail

PROJ=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
READS=$PROJ/rnaseq_reads
HDR=$PROJ/nextflow/bin/fix_fastq_header_trinity
APPLY=0
if [[ "${1:-}" == "--apply" ]]; then APPLY=1; shift; fi
[[ $# -gt 0 ]] || { echo "no species given" >&2; exit 2; }
[[ -x $HDR ]] || { echo "missing $HDR" >&2; exit 1; }
PIGZ=$(command -v pigz || true)
[[ -n $PIGZ ]] || { echo "pigz not found" >&2; exit 1; }

first_hdr() { zcat "$1" | head -1 || true; }
count_reads() { zcat "$1" | awk 'END{print NR/4}'; }

for sp in "$@"; do
  for r in 1 2; do
    f=$READS/${sp}_norm_R${r}.fastq.gz
    [[ -s $f ]] || { echo "FAIL $sp R$r: missing or empty $f" >&2; exit 1; }
    h=$(first_hdr "$f")
    if [[ "$h" == */$r ]]; then echo "SKIP $sp R$r: already has /$r ($h)"; continue; fi
    if [[ -e $f.pre_rename ]]; then echo "FAIL $sp R$r: $f.pre_rename exists" >&2; exit 1; fi
    echo "TODO $sp R$r: $h"
    [[ $APPLY -eq 1 ]] || continue
    tmp=$f.tmp_rename
    n_in=$(count_reads "$f")
    "$HDR" --read "$r" "$f" | "$PIGZ" -c > "$tmp"
    n_out=$(count_reads "$tmp")
    h2=$(first_hdr "$tmp")
    [[ "$n_in" == "$n_out" ]] || { echo "FAIL $sp R$r: reads $n_in -> $n_out" >&2; rm -f "$tmp"; exit 1; }
    [[ "$h2" == */$r ]] || { echo "FAIL $sp R$r: new header $h2" >&2; rm -f "$tmp"; exit 1; }
    mv "$f" "$f.pre_rename"
    mv "$tmp" "$f"
    echo "DONE $sp R$r: $n_in reads; $h2"
  done
done
