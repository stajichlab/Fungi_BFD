#!/usr/bin/bash -l
#
# run_unified.sh — Unified pipeline: samples.csv → annotated genomes
#
# Runs four phases in sequence, each with its own params file and -resume:
#
#   Phase 0: Genome cleaning (GENOME_CLEAN / GENOME_CLEAN_BATCH, --only_clean)
#            → input_clean_genomes/<asmid>.fa.gz (filtered + purged, unmasked)
#
#   Phase 1: BFD genome stats (ASM_STATS + BUSCO + TELOMERES + BFD_MERGE)
#            → tables/asm_stats.parquet, tables/busco_genome.parquet
#            Reads input_clean_genomes/<asmid>.fa.gz from Phase 0 (raw genome,
#            no masking/annotation required — see resolveGenomeFile()).
#
#   Phase 2: compare_ANI (ANI compute + PICK_REPRESENTATIVE_STRAIN)
#            → genome_annotation/_reuse_assignments/abinitio_reuse_assignments.csv
#            Also reads input_clean_genomes/<asmid>.fa.gz from Phase 0.
#
#   Phase 3: funannotate (MASK → RNASEQ → TRAIN → PREDICT → ANNOTATE)
#            → annotated genomes in genome_annotation/
#            CLEAN is a no-op here: every genome was already cleaned in Phase 0,
#            so FUNANNOTATE_GENOME_PREP's jobs_to_clean filter skips straight to
#            MASKREPEAT_TANTAN_RUN (tantan soft-masking, the default; an earlgrey
#            alternative can be spliced in here later via --pipeline earlgrey_mask).
#
# Prerequisites:
#   - samples.csv in the project root
#   - Genome source directory populated (see params_unified_funannotate.yaml)
#   - Nextflow module loaded (this script does it for you)
#
# Usage:
#   sbatch nextflow/run_unified.sh                    # Full run
#   sbatch nextflow/run_unified.sh --n_test 5         # Test with 5 samples
#   sbatch nextflow/run_unified.sh --skip_phase0      # Skip CLEAN (already done)
#   sbatch nextflow/run_unified.sh --skip_phase1      # Skip BFD (already done)
#   sbatch nextflow/run_unified.sh --skip_phase2      # Skip ANI (already done)
#   sbatch nextflow/run_unified.sh --only_phase3      # Only run funannotate
#
# SBART header:
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00 -p batch
#SBATCH --job-name=nxf_unified
#SBATCH --output=logs/unified_launch.%j.log

set -euo pipefail

# ── Resolve project root ────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_DIR"

# ── Load Nextflow ───────────────────────────────────────────────────────────
source /etc/profile.d/modules.sh 2>/dev/null || true
module load nextflow

mkdir -p logs/nextflow

# ── Params files ────────────────────────────────────────────────────────────
# Phase 0 reuses FUN_PARAMS with --only_clean true tacked on (see run_phase call
# below) rather than a separate yaml — CLEAN and MASK live in the same
# FUNANNOTATE_GENOME_PREP subworkflow, so one params file keeps both phases'
# clean_batch_size/source/etc settings in sync.
CLEAN_PARAMS="nextflow/params_unified_clean.yaml"
BFD_PARAMS="nextflow/params_unified_bfd.yaml"
ANI_PARAMS="nextflow/params_unified_ani.yaml"
FUN_PARAMS="nextflow/params_unified_funannotate.yaml"

# ── Parse command-line args ─────────────────────────────────────────────────
SKIP_PHASE0=false
SKIP_PHASE1=false
SKIP_PHASE2=false
SKIP_PHASE3=false
EXTRA_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip_phase0)     SKIP_PHASE0=true;  shift ;;
        --skip_phase1)     SKIP_PHASE1=true;  shift ;;
        --skip_phase2)     SKIP_PHASE2=true;  shift ;;
        --skip_phase3)     SKIP_PHASE3=true;  shift ;;
        --only_phase0)     SKIP_PHASE1=true; SKIP_PHASE2=true; SKIP_PHASE3=true; shift ;;
        --only_phase1)     SKIP_PHASE0=true; SKIP_PHASE2=true; SKIP_PHASE3=true; shift ;;
        --only_phase2)     SKIP_PHASE0=true; SKIP_PHASE1=true; SKIP_PHASE3=true; shift ;;
        --only_phase3)     SKIP_PHASE0=true; SKIP_PHASE1=true; SKIP_PHASE2=true; shift ;;
        -resume|-r)        EXTRA_ARGS+=("-resume"); shift ;;
        -stub-run|-s)      EXTRA_ARGS+=("-stub-run"); shift ;;
        --help|-h)
            head -38 "$0"
            exit 0
            ;;
        *)
            EXTRA_ARGS+=("$1")
            shift
            ;;
    esac
done

# ── Helper: run a phase ─────────────────────────────────────────────────────
run_phase() {
    local phase_name="$1"
    local profile="$2"
    local params_file="$3"
    local pipeline="$4"
    shift 4
    local phase_args=("$@")   # extra CLI overrides specific to this phase (e.g. --only_clean true)

    echo ""
    echo "═══════════════════════════════════════════════════════════════════════════"
    echo "  Phase: ${phase_name}"
    echo "  Profile: ${profile}"
    echo "  Params: ${params_file}"
    echo "  Pipeline: ${pipeline}"
    echo "═══════════════════════════════════════════════════════════════════════════"
    echo ""

    NXF_OPTS="-Xms512m -Xmx4g" \
    nextflow run nextflow/main.nf \
        -c nextflow/nextflow.config \
        -profile "${profile}" \
        -params-file "${params_file}" \
        -resume \
        "${phase_args[@]}" \
        "${EXTRA_ARGS[@]}"

    local exit_code=$?
    if [ $exit_code -ne 0 ]; then
        echo ""
        echo "✗ Phase ${phase_name} failed with exit code ${exit_code}"
        echo "  Fix the issue and re-run with --skip_phase${phase_num} or --only_phase${phase_num}"
        exit $exit_code
    fi

    echo ""
    echo "✓ Phase ${phase_name} completed successfully"
}

# ── Phase 0: genome cleaning ─────────────────────────────────────────────────
if [ "$SKIP_PHASE0" = false ]; then
    echo ""
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    echo "  UNIFIED PIPELINE: samples.csv → annotated genomes"
    echo "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
    run_phase "0: Genome cleaning (GENOME_CLEAN/GENOME_CLEAN_BATCH)" "funannotate" "$CLEAN_PARAMS" "funannotate" 
else
    echo "Skipping Phase 0 (genome cleaning)"
fi

# ── Phase 1: BFD genome stats ───────────────────────────────────────────────
if [ "$SKIP_PHASE1" = false ]; then
    run_phase "1: BFD genome stats (ASM_STATS + BUSCO + MERGE)" "BFD" "$BFD_PARAMS" "BFD"
else
    echo "Skipping Phase 1 (BFD genome stats)"
fi

# ── Phase 2: compare_ANI ────────────────────────────────────────────────────
if [ "$SKIP_PHASE2" = false ]; then
    run_phase "2: compare_ANI (ANI + representative selection)" "ani" "$ANI_PARAMS" "compare_ani"
else
    echo "Skipping Phase 2 (compare_ANI)"
fi

# ── Phase 3: funannotate ────────────────────────────────────────────────────
if [ "$SKIP_PHASE3" = false ]; then
    run_phase "3: funannotate (MASK → RNASEQ → PREDICT → ANNOTATE)" "funannotate" "$FUN_PARAMS" "funannotate"
else
    echo "Skipping Phase 3 (funannotate)"
fi

echo ""
echo "═══════════════════════════════════════════════════════════════════════════"
echo "  Unified pipeline completed successfully!"
echo "═══════════════════════════════════════════════════════════════════════════"
echo ""
