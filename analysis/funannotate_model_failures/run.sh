#!/usr/bin/env bash
# Reproduce the funannotate "too few training models" failure analysis.
set -euo pipefail
cd "$(dirname "$0")/../.."   # project root

python3 analysis/funannotate_model_failures/parse_predict_failures.py \
    --log-roots genome_annotation do_annotation/genome_annotation \
    --samples   samples.csv \
    --asm-stats results/genome_stats/asm_stats \
    --outdir    analysis/funannotate_model_failures \
    --phylum    Basidiomycota

# Compute fresh seqkit metrics for the failing genomes (most are absent from the
# stale asm_stats table) and apply the screen rule.
module load seqkit 2>/dev/null || true
OUT=analysis/funannotate_model_failures/failed_genome_asmstats.tsv
printf "ASMID\ttotal_bp\tcontigs\tN50\tmax_len\n" > "$OUT"
for a in $(tail -n +2 analysis/funannotate_model_failures/failed_genomes.tsv | cut -f1 | sort -u); do
  f="input_clean_genomes/${a}.fa.gz"; [ -f "$f" ] || f="input_clean_genomes/${a}.masked.fasta.gz"
  [ -f "$f" ] || continue
  seqkit stats -a -T "$f" 2>/dev/null | tail -1 | \
    awk -v A="$a" 'BEGIN{FS=OFS="\t"} {print A, $5, $4, $13, $8}' >> "$OUT"
done
echo "Done. See crosstab_summary.txt and failed_genomes_annotated.tsv"
