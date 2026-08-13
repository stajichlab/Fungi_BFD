#!/bin/bash
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 6:00:00
#SBATCH --job-name=gmk_ab
#SBATCH --output=logs/%x.%j.log
# Rerun funannotate predict for one genome, with GeneMark either at its normal
# production weight (baseline) or forced to 0 (nogenemark), against the SAME
# already-completed PASA training data used by production. Both variants use
# the same funannotate/dev-1.9 host module as production FUNANNOTATE_PREDICT
# (nextflow/modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf) so this
# is apples-to-apples with what's already in genome_annotation/.
#
# Required env (pass via --export=ALL,OUT=...,ASMID=...,...):
#   OUT           species/strain dir name (matches genome_annotation/<OUT>)
#   ASMID         assembly id (matches input_clean_genomes/<ASMID>.masked.fasta.gz)
#   SPECIES       e.g. "Penicillium citrinum"
#   STRAIN        e.g. "NRRL 1841"
#   LOCUSTAG      e.g. FB1A4A3B
#   BUSCO         e.g. dikarya
#   TRANSL_TABLE  e.g. 1
#   MODE          baseline | nogenemark
set -euo pipefail

PROJECT_ROOT=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
RUNS_ROOT=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
ANALYSIS_DIR="$PROJECT_ROOT/analysis/genemark_es_contribution"
GENOME_GZ="$RUNS_ROOT/input_clean_genomes/${ASMID}.masked.fasta.gz"
TRAIN_DIR="$PROJECT_ROOT/genome_annotation_training/${OUT}/training"
PREDICTDIR="$ANALYSIS_DIR/outputs/predict_runs/${OUT}__${MODE}"

[ -s "$GENOME_GZ" ] || { echo "ERROR: genome not found: $GENOME_GZ" >&2; exit 1; }
[ -d "$TRAIN_DIR" ] || { echo "ERROR: training dir not found: $TRAIN_DIR" >&2; exit 1; }

rm -rf "$PREDICTDIR"
mkdir -p "$PREDICTDIR"
ln -sfn "$TRAIN_DIR" "$PREDICTDIR/training"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load miniconda3
eval "$(conda shell.bash hook)"
module load funannotate/dev-1.9

# PRIVATE writable copy of the shared AUGUSTUS_CONFIG_PATH, one per job (not
# shared across OUT/MODE combinations). predict trains a new species config
# into $AUGUSTUS_CONFIG_PATH/species/<name>/, and (a) the live shared dir is
# written to concurrently by real production jobs (do_annotation_asco etc.)
# right now, so we never point at it directly, and (b) baseline/nogenemark
# for the SAME species would otherwise race on the same species/ subdir if
# they shared one copy and happened to run concurrently. 39M, cheap to copy
# per job. The leaf directory MUST be literally named "config" -- funannotate
# predict.py only assigns AUGUSTUS_BASE (needed later for bam2hints/gff2gb
# script paths) when `os.path.basename(AUGUSTUS_CONFIG_PATH) == "config"`;
# anything else silently skips that assignment and later crashes with
# UnboundLocalError: local variable 'AUGUSTUS_BASE' referenced before
# assignment (hit this on the first real submission of this job).
SHARED_AUGCFG=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/lib/augustus/3.5/config
PRIVATE_AUGCFG="$PREDICTDIR.augustus_config/config"
rm -rf "$(dirname "$PRIVATE_AUGCFG")"
mkdir -p "$(dirname "$PRIVATE_AUGCFG")"
cp -r "$SHARED_AUGCFG" "$PRIVATE_AUGCFG"
export AUGUSTUS_CONFIG_PATH="$PRIVATE_AUGCFG"
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
export TMPDIR="${SCRATCH:-/tmp}/gmk_ab.${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"

GENOME_IN="$TMPDIR/${OUT}.genome_input.fa"
pigz -dc "$GENOME_GZ" > "$GENOME_IN"

WEIGHT_ARGS=(-w codingquarry:0 glimmerhmm:0)
if [ "$MODE" = "nogenemark" ]; then
    WEIGHT_ARGS=(-w codingquarry:0 glimmerhmm:0 genemark:0)
fi

echo "[INFO] MODE=$MODE OUT=$OUT weights=${WEIGHT_ARGS[*]}"
time funannotate predict --name "$LOCUSTAG" -i "$GENOME_IN" --strain "$STRAIN" \
    -o "$PREDICTDIR" -s "$SPECIES" --cpu "${SLURM_CPUS_PER_TASK:-16}" --busco_db "$BUSCO" \
    --AUGUSTUS_CONFIG_PATH "$AUGUSTUS_CONFIG_PATH" "${WEIGHT_ARGS[@]}" \
    --min_training_models 30 --tmpdir "$TMPDIR" --SeqCenter NCBI \
    --keep_no_stops --header_length 24 --protein_evidence "$PROJECT_ROOT/lib/swissprot_fungi.faa" \
    --max_intronlen 3000 --min_intronlen 10 \
    --tbl2asn "-l paired-ends" --table "$TRANSL_TABLE" --auto-skip-genemark

echo "[INFO] Done. GBK: $PREDICTDIR/predict_results/${OUT}.gbk"
rm -rf "$TMPDIR" "$(dirname "$PRIVATE_AUGCFG")"
