#!/usr/bin/bash -l
#SBATCH -p epyc -c 4 --mem 8G -t 12:00:00 -J cleanup_arms
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/predict_arms/logs/%x_%j.out
# Usage: sbatch cleanup_predict_arms.sh [--apply]   (no argument = dry run, writes the plan only)
set -euo pipefail
A=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/predict_arms
module load zstd 2>/dev/null || true
/usr/bin/python3.12 $A/cleanup_predict_arms.py $A $A/cleanup_record "$@"
