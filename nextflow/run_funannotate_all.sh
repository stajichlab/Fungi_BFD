#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 2 --mem 8gb --time 7-00:00:00
#SBATCH --job-name=nxf_funannotate
#SBATCH --output=logs/funannotate_launch.%j.log

# Thin wrapper kept for backward compatibility -- see run_funannotate.sh,
# which toggles between representatives-only and all-strains passes:
#   sbatch nextflow/run_funannotate.sh all [extra args...]
set -euo pipefail
exec "$(dirname "${BASH_SOURCE[0]}")/run_funannotate.sh" all "$@"
