#!/usr/bin/bash -l
#SBATCH -p epyc -c 32 --mem 32G --time 2:00:00 -J yj2025_aai
#SBATCH -o yj2025_aai.%j.log
# Placement check for "Saccharomycotina sp. YJ-2025a/b" (FF90B97F, F6EB82EF): one-way AAI of each
# proteome against a panel of Pucciniomycotina genus representatives, every Rhodotorula species and
# a few Ascomycota/Basidiomycota controls. R. pacifica (FBBC614A) is run as a third query to calibrate.
# See YJ2025_TAXONOMY_CHECK.md.
#
# Reads BFD tables and input/pep (read only). Writes only to the submit directory ($PWD) and $SCRATCH.
# Submit from a purgeable work folder, for example:
#   mkdir -p ~/projects/BFD_dashboard/bench/yj2025_taxonomy && cd $_
#   sbatch /path/to/Fungi_BFD/analysis/yj2025_taxonomy_check/run.sh
set -euo pipefail
HERE="${HERE:-$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")}"
[ -f "$HERE/select_panel.py" ] || HERE="$SLURM_SUBMIT_DIR/scripts"   # sbatch copies the script; scripts/ beside it
RUNS="${RUNS:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs}"
PY=/usr/bin/python3.12
CPU="${SLURM_CPUS_ON_NODE:-4}"
OUT="$PWD"
WORK="${SCRATCH:-$OUT/tmp}/yj2025_aai"
mkdir -p "$WORK/db" "$WORK/q" "$OUT/hits"
module load diamond/2.1.7
diamond version

$PY "$HERE/select_panel.py" --runs "$RUNS" --out "$OUT/panel.tsv"
echo "panel: $(($(wc -l < "$OUT/panel.tsv") - 1)) genomes"

# protein counts and DIAMOND databases (one per genome)
: > "$OUT/nprot.tsv"
# cut, not read: empty TSV fields collapse under IFS=tab
tail -n +2 "$OUT/panel.tsv" | cut -f1,11 | while IFS=$'\t' read -r loc pep; do
  zcat "$pep" | sed 's/\*$//' > "$WORK/q/$loc.fa"
  printf '%s\t%s\n' "$loc" "$(grep -c '^>' "$WORK/q/$loc.fa")" >> "$OUT/nprot.tsv"
done
tail -n +2 "$OUT/panel.tsv" | cut -f1 | xargs -P "$CPU" -I{} \
  diamond makedb --quiet --in "$WORK/q/{}.fa" -d "$WORK/db/{}"

# every query against every panel genome (queries against each other too)
QUERIES=$(awk -F'\t' '$2=="query"{print $1}' "$OUT/panel.tsv")
JOBS="$WORK/jobs.txt"; : > "$JOBS"
for q in $QUERIES; do
  for ref in $(tail -n +2 "$OUT/panel.tsv" | cut -f1); do
    [ "$q" = "$ref" ] || echo "$q $ref" >> "$JOBS"
  done
done
echo "$(wc -l < "$JOBS") pairwise runs"
PAR=$(( CPU / 4 )); [ "$PAR" -ge 1 ] || PAR=1
xargs -P "$PAR" -L 1 sh -c 'diamond blastp --quiet --threads 4 --sensitive --max-target-seqs 1 --evalue 1e-5 \
  -q "'"$WORK"'/q/$0.fa" -d "'"$WORK"'/db/$1" -o "'"$OUT"'/hits/$0__$1.tsv" \
  --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore' < "$JOBS"

$PY "$HERE/aai_summary.py" --panel "$OUT/panel.tsv" --hits "$OUT/hits" --nprot "$OUT/nprot.tsv" --out "$OUT/aai.tsv"
for q in $QUERIES; do echo "== $q top 12"; awk -F'\t' -v q="$q" 'NR==1 || ($1==q && $3<=12)' "$OUT/aai.tsv" | cut -f3-13; done
tar czf "$OUT/hits.tar.gz" -C "$OUT" hits && rm -rf "$OUT/hits"
echo done "$(date -Is)"
