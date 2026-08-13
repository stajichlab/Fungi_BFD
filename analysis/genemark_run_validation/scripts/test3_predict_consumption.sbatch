#!/bin/bash
#SBATCH -N 1 -n 1 -c 16 --mem 32gb --time 6:00:00
#SBATCH --job-name=gmrun_test3
#SBATCH --output=logs/%x.%j.log
# Test 3: real funannotate predict --genemark_gtf consumption of a GTF
# produced by the standalone GENEMARK_RUN module (test 1's output). Mirrors
# FUNANNOTATE_PREDICT/main.nf's actual command construction (same params,
# same training data reuse) so this validates predict's real handling of
# --genemark_gtf, not a hand-simplified stand-in. Pass criteria:
#   1. Command completes and produces predict_results/<out>.gbk
#   2. Log does NOT show "Running GeneMark-ES on assembly" (confirms
#      --genemark_gtf short-circuits predict's own internal GeneMark call,
#      per the design doc's open question -- this is what actually answers it)
#   3. Final gene count is in the same ballpark as the A/B test's baseline
#      (11,198 for this genome) -- confirms the externally-supplied GTF
#      produces an equivalent annotation, not a degraded one.
set -euo pipefail

PROJECT_ROOT=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
RUNS_ROOT=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
VALDIR="$PROJECT_ROOT/analysis/genemark_run_validation"
OUT=Penicillium_citrinum_NRRL_1841_gmtest3
GENOME_GZ="$RUNS_ROOT/input_clean_genomes/GCA_020284165.1_ASM2028416v1.masked.fasta.gz"
TRAIN_DIR="$PROJECT_ROOT/genome_annotation_training/Penicillium_citrinum_NRRL_1841/training"
GENEMARK_GTF="$VALDIR/outputs/Penicillium_citrinum_NRRL_1841_gmtest.genemark.gtf"
PREDICTDIR="$VALDIR/outputs/predict_runs/${OUT}"

[ -s "$GENOME_GZ" ] || { echo "ERROR: genome not found: $GENOME_GZ" >&2; exit 1; }
[ -d "$TRAIN_DIR" ] || { echo "ERROR: training dir not found: $TRAIN_DIR" >&2; exit 1; }
[ -s "$GENEMARK_GTF" ] || { echo "ERROR: genemark gtf not found (run test 1 first): $GENEMARK_GTF" >&2; exit 1; }

rm -rf "$PREDICTDIR"
mkdir -p "$PREDICTDIR"
ln -sfn "$TRAIN_DIR" "$PREDICTDIR/training"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load miniconda3
eval "$(conda shell.bash hook)"
module load funannotate/dev-1.9

# Private per-job Augustus config copy, same reasoning as the earlier A/B
# test (analysis/genemark_es_contribution/) -- never write into the live
# shared config while real production jobs may be using it concurrently.
SHARED_AUGCFG=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/lib/augustus/3.5/config
PRIVATE_AUGCFG="$PREDICTDIR.augustus_config/config"
rm -rf "$(dirname "$PRIVATE_AUGCFG")"
mkdir -p "$(dirname "$PRIVATE_AUGCFG")"
cp -r "$SHARED_AUGCFG" "$PRIVATE_AUGCFG"
export AUGUSTUS_CONFIG_PATH="$PRIVATE_AUGCFG"
export FUNANNOTATE_DB=/bigdata/stajichlab/shared/lib/funannotate_db
export TMPDIR="${SCRATCH:-/tmp}/gmrun_test3.${SLURM_JOB_ID:-$$}"
mkdir -p "$TMPDIR"

GENOME_IN="$TMPDIR/${OUT}.genome_input.fa"
pigz -dc "$GENOME_GZ" > "$GENOME_IN"

echo "[INFO] Running funannotate predict with --genemark_gtf $GENEMARK_GTF"
funannotate predict --name GMTEST3 -i "$GENOME_IN" --strain "NRRL 1841" \
    -o "$PREDICTDIR" -s "Penicillium citrinum" --cpu "${SLURM_CPUS_PER_TASK:-16}" \
    --busco_db dikarya \
    --AUGUSTUS_CONFIG_PATH "$AUGUSTUS_CONFIG_PATH" -w codingquarry:0 glimmerhmm:0 \
    --min_training_models 30 --tmpdir "$TMPDIR" --SeqCenter NCBI \
    --keep_no_stops --header_length 24 --protein_evidence "$PROJECT_ROOT/lib/swissprot_fungi.faa" \
    --max_intronlen 3000 --min_intronlen 10 \
    --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark \
    --genemark_gtf "$GENEMARK_GTF"

echo "[INFO] Done."
echo "[CHECK] Internal GeneMark-ES invocation avoided?"
if grep -q "Running GeneMark-ES on assembly" "$PREDICTDIR/logfiles/funannotate-predict.log" 2>/dev/null; then
    echo "[FAIL] predict ran its own internal GeneMark-ES despite --genemark_gtf -- --genemark_gtf did NOT short-circuit as expected" >&2
else
    echo "[PASS] predict did not run its own internal GeneMark-ES -- --genemark_gtf correctly short-circuited the internal call"
fi
echo "[CHECK] GeneMark GTF actually used in log?"
grep -i "genemark" "$PREDICTDIR/logfiles/funannotate-predict.log" | head -10 || true
echo "[CHECK] Final gene count vs A/B baseline (11198 for this genome):"
grep -m1 "Collecting final annotation files" "$PREDICTDIR/logfiles/funannotate-predict.log" || true

rm -rf "$TMPDIR" "$(dirname "$PRIVATE_AUGCFG")"
