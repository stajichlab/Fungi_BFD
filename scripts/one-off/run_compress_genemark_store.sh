#!/usr/bin/bash -l
#SBATCH -N 1 -n 1 -c 8 --mem 8gb --time 12:00:00 -p stajichlab
#SBATCH --job-name=gm_store_gz
#SBATCH --output=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/analysis/genemark_store_census/compress_%j.log
# Compress <store>/<out>/<out>.genemark.gtf in place (see compress_genemark_store.py).
# Safe to resubmit: finished folders are skipped. Submit:  sbatch run_compress_genemark_store.sh
set -euo pipefail
STORE=${STORE:-/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genemark_store}
/usr/bin/python3.12 /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD/scripts/one-off/compress_genemark_store.py \
    --store "$STORE" --apply --remove-superseded --workers 8 --min-age-min 60
