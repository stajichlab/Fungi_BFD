#!/usr/bin/bash -l
#SBATCH -p highclock -N 1 -n 4 -c 1 --mem 16gb --time 4:00:00 --job-name pfam_ab_mpi
#SBATCH --out logs/pfam_ab_mpi.%j.log

set -euo pipefail
source /etc/profile.d/modules.sh 2>/dev/null || true
module load hmmer/3.4-mpi
module load db-pfam

OUTDIR="analysis/pfam_hmmsearch_perf"
PROT="input/pep/Malassezia_brasiliensis_CBS_14135.proteins.fa"

echo "=== Host: $(hostname) ==="
echo "=== PFAM_DB: $PFAM_DB ==="

echo ""
echo "=== TEST C: hmmsearch --mpi with 4 MPI tasks (pfam_tasks=4, pfam_nodes=1) ==="
/usr/bin/time -v srun -N 1 -n 4 hmmsearch --mpi --cut_ga --noali \
    --domtbl "$OUTDIR/testC_mpi4.domtblout" \
    --tblout "$OUTDIR/testC_mpi4.tblout" \
    "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testC_mpi4.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testC_mpi4.timing.txt"

echo ""
echo "=== TEST D: hmmsearch --mpi with 2 MPI tasks (pfam_tasks=2, pfam_nodes=1) ==="
/usr/bin/time -v srun -N 1 -n 2 hmmsearch --mpi --cut_ga --noali \
    --domtbl "$OUTDIR/testD_mpi2.domtblout" \
    --tblout "$OUTDIR/testD_mpi2.tblout" \
    "$PFAM_DB/Pfam-A.hmm" "$PROT" > /dev/null 2> "$OUTDIR/testD_mpi2.timing.txt"
grep -E "Elapsed|Maximum resident" "$OUTDIR/testD_mpi2.timing.txt"

echo ""
echo "=== Sanity: domain counts vs Test A baseline (single-task, cpu=4 threads) ==="
grep -vc '^#' "$OUTDIR/testC_mpi4.tblout" || true
grep -vc '^#' "$OUTDIR/testD_mpi2.tblout" || true

echo ""
echo "=== DONE ==="
