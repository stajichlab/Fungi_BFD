#!/usr/bin/bash -l
#SBATCH -p highclock -N 1 -n 4 -c 1 --mem 2gb --time 10:00 --job-name cpubind_test
#SBATCH --out logs/cpubind_test.%j.log

source /etc/profile.d/modules.sh 2>/dev/null || true
module load openmpi/4.1.2_slurm-24.11.1_mpi1-compat

BIN=analysis/pfam_hmmsearch_perf/mpi_sanity/hello

echo "=== scontrol show job (allocation detail) ==="
scontrol show job $SLURM_JOB_ID

echo ""
echo "=== srun --cpu-bind=none -n 4 ==="
srun --cpu-bind=none -n 4 $BIN || echo "FAILED exit $?"

echo ""
echo "=== srun --cpu-bind=threads -n 4 ==="
srun --cpu-bind=threads -n 4 $BIN || echo "FAILED exit $?"

echo ""
echo "=== srun --cpu-bind=none --mpi=pmix -n 4 ==="
srun --cpu-bind=none --mpi=pmix -n 4 $BIN || echo "FAILED exit $?"
