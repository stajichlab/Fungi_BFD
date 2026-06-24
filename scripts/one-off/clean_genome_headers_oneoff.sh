#!/usr/bin/env bash
#
# One-off: strip multi-word FASTA headers from recently-added genomes.
#
# Looks in input_clean_genomes/ for files newer than 2026-06-19. For each, reads
# the first line; if the header is a single word (e.g. ">ABCD123") it is already
# clean and skipped. If the header has a description (e.g. ">ABCD124 there are
# more words") the file is passed through scripts/clean_genome_fa.py, which keeps
# only the first whitespace-delimited token of every header. The cleaned file is
# recompressed with `pigz -n`, given the original file's mtime, and then placed.
#
# MODES (arg 1, default: dryrun):
#   dryrun         overwrite-style dry run  -- print actions, change nothing
#   dryrun-backup  backup-style dry run     -- print actions, change nothing
#   overwrite      REAL: cleaned file replaces the original in place
#   backup         REAL: original moved to cleaned_backup/, cleaned file placed
#
# Usage:  scripts/clean_genome_headers_oneoff.sh [MODE]
set -euo pipefail

DIR="input_clean_genomes"
NEWER_THAN="2026-06-19"
CLEANER="scripts/clean_genome_fa.py"
BACKUP_DIR="${DIR}/cleaned_backup"
MODE="${1:-dryrun}"

case "$MODE" in
    dryrun|dryrun-backup|overwrite|backup) ;;
    *) echo "ERROR: unknown mode '$MODE' (use: dryrun|dryrun-backup|overwrite|backup)" >&2; exit 2 ;;
esac

: "${SCRATCH:?SCRATCH is not set}"
[[ -d "$DIR" ]]      || { echo "ERROR: $DIR not found (run from project root)" >&2; exit 1; }
[[ -f "$CLEANER" ]]  || { echo "ERROR: $CLEANER not found" >&2; exit 1; }
command -v pigz >/dev/null || { echo "ERROR: pigz not on PATH" >&2; exit 1; }

is_real=0
[[ "$MODE" == overwrite || "$MODE" == backup ]] && is_real=1
[[ "$MODE" == backup ]] && { [[ $is_real -eq 1 ]] && mkdir -p "$BACKUP_DIR"; }

echo "== mode=$MODE  dir=$DIR  newer-than=$NEWER_THAN  scratch=$SCRATCH =="

n_total=0; n_clean=0; n_dirty=0

# -maxdepth 1 -type f matches your find exactly; subdirs (e.g. clean/) excluded.
while IFS= read -r file; do
    n_total=$((n_total + 1))
    base="$(basename "$file")"

    # First line = the FASTA header. NF==1 -> single word (already clean).
    header="$(zcat "$file" | head -n 1 || true)"
    nfields="$(printf '%s' "$header" | awk '{print NF}')"

    if [[ "${nfields:-0}" -le 1 ]]; then
        n_clean=$((n_clean + 1))
        # Uncomment for verbose skip logging:
        # echo "SKIP  $base  (single-word header: $header)"
        continue
    fi

    n_dirty=$((n_dirty + 1))
    tmp_plain="$SCRATCH/$(basename "$file" .gz)"   # e.g. GCA_x.fa  or  GCA_x.masked.fasta
    tmp_gz="${tmp_plain}.gz"

    echo "CLEAN $base"
    echo "      header: $header"

    if [[ $is_real -eq 0 ]]; then
        echo "      would: pigz -dc '$file' | python $CLEANER > '$tmp_plain'"
        echo "      would: pigz -n '$tmp_plain'  (-> $tmp_gz)"
        echo "      would: touch -r '$file' '$tmp_gz'"
        if [[ "$MODE" == dryrun-backup ]]; then
            echo "      would: mv '$file' '$BACKUP_DIR/$base'  (preserve original)"
            echo "      would: mv '$tmp_gz' '$file'"
        else
            echo "      would: mv '$tmp_gz' '$file'  (overwrite original)"
        fi
        continue
    fi

    # ---- real work ----
    pigz -dc "$file" | python "$CLEANER" > "$tmp_plain"
    pigz -n -f "$tmp_plain"               # produces $tmp_gz
    touch -r "$file" "$tmp_gz"            # match original mtime
    if [[ "$MODE" == backup ]]; then
        mv "$file" "$BACKUP_DIR/$base"
    fi
    mv "$tmp_gz" "$file"
    echo "      done -> $file"
done < <(find "$DIR" -maxdepth 1 -type f -newermt "$NEWER_THAN")

echo "== scanned=$n_total  already-clean=$n_clean  cleaned=$n_dirty (mode=$MODE) =="
