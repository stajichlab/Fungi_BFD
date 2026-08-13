#!/bin/bash
# Submit the GeneMark-ES contribution A/B test: for each of 3 candidate
# genomes (real training data + already-completed baseline predict, picked
# from analysis/funannotate_predict_stage_timing/outputs/per_run_summary.csv
# for fast genemark_es_train_seconds so the pair fits in a few hours each),
# submit two funannotate predict SLURM jobs -- baseline (production weights)
# and nogenemark (-w genemark:0) -- against the SAME existing PASA training
# data. See GENEMARK_ES_CONTRIBUTION.md for the full design and how to run
# scripts/compare_results.py once all 6 jobs finish.
set -euo pipefail
cd "$(dirname "$0")"

mkdir -p outputs/predict_runs logs

# OUT|ASMID|SPECIES|STRAIN|LOCUSTAG|BUSCO|TRANSL_TABLE
GENOMES=(
  "Penicillium_citrinum_NRRL_1841|GCA_020284165.1_ASM2028416v1|Penicillium citrinum|NRRL 1841|FB1A4A3B|dikarya|1"
  "Saccharomyces_kudriavzevii_IFO10991|GCA_000257085.1_Saccharomyces_kudriavzevii_strain_IFO10991_v1.0|Saccharomyces kudriavzevii|IFO10991|F620E934|dikarya|1"
  "Kluyveromyces_marxianus_YG-4|GCA_053539435.1_ASM5353943v1|Kluyveromyces marxianus|YG-4|F25AF76C|dikarya|1"
)

for row in "${GENOMES[@]}"; do
  IFS='|' read -r OUT ASMID SPECIES STRAIN LOCUSTAG BUSCO TRANSL_TABLE <<< "$row"
  for MODE in baseline nogenemark; do
    jobname="gmk_ab_${OUT}_${MODE}"
    echo "Submitting $jobname"
    sbatch --job-name="$jobname" \
      --export=ALL,OUT="$OUT",ASMID="$ASMID",SPECIES="$SPECIES",STRAIN="$STRAIN",LOCUSTAG="$LOCUSTAG",BUSCO="$BUSCO",TRANSL_TABLE="$TRANSL_TABLE",MODE="$MODE" \
      scripts/run_predict_variant.sbatch
  done
done

echo
echo "Submitted. Check status with: squeue -u \$USER -n gmk_ab_*"
echo "Once all 6 finish, run: python3 scripts/compare_results.py"
