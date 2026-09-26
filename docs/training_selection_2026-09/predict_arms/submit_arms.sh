#!/usr/bin/bash
# Submit all predict arms with dependencies. Usage: bash submit_arms.sh BENCH_JOB PREP_JOB_NC PREP_JOB_AN PREP_JOB_BC
# Writes the submitted job ids to jobs.tsv.
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
BENCH=${1:-}
# (prep job ids are no longer used; pass DEP instead)
SB="env -u SLURM_JOB_ID -u SLURM_JOBID -u SLURM_ARRAY_TASK_ID sbatch --parsable"
# GENOMES (env) overrides the genome list; DEP (env) adds a dependency to every step A job; JOBS (env) names the job table
JOBS=${JOBS:-jobs.tsv}
echo -e "genome\tarm\tstep\tevid\tjob" > $JOBS
for GEN in ${GENOMES:-Neurospora_crassa_OR74A Aspergillus_nidulans_FGSC_A4 Botrytis_cinerea_B05.10}; do
  for ARM in old oldnew tx cdsS cdsB busco; do
    ja=$($SB -J A_${ARM}_${GEN:0:10} ${DEP:+--dependency=afterok:$DEP} arm.sh $GEN $ARM A)
    echo -e "$GEN\t$ARM\tA\t-\t$ja" >> $JOBS
    jb=$($SB -J B_${ARM}_${GEN:0:10} --dependency=afterok:$ja arm.sh $GEN $ARM B fixed)
    echo -e "$GEN\t$ARM\tB\tfixed\t$jb" >> $JOBS
    if [[ $ARM == tx || $ARM == cdsS || $ARM == cdsB ]]; then
      jo=$($SB -J O_${ARM}_${GEN:0:10} --dependency=afterok:$ja arm.sh $GEN $ARM B own)
      echo -e "$GEN\t$ARM\tB\town\t$jo" >> $JOBS
    fi
  done
done
cat $JOBS
