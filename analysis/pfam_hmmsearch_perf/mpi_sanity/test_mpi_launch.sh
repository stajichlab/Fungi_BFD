#!/usr/bin/bash -l
#SBATCH -p highclock -N 1 -n 4 -c 1 --mem 2gb --time 10:00 --job-name mpi_sanity
#SBATCH --out logs/mpi_sanity.%j.log

source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.2_slurm-24.11.1_mpi1-compat

BIN=analysis/pfam_hmmsearch_perf/mpi_sanity/hello

echo "=== plain srun -n 4 (no explicit --mpi=) ==="
srun -n 4 $BIN

echo ""
echo "=== srun --mpi=pmix -n 4 ==="
srun --mpi=pmix -n 4 $BIN

echo ""
echo "=== srun --mpi=pmi2 -n 4 ==="
srun --mpi=pmi2 -n 4 $BIN || echo "pmi2 FAILED (exit $?)"
