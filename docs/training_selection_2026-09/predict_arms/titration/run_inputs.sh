#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]
#SBATCH -p epyc -c 4 --mem 32G -t 2:00:00 -J titrIn -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/predict_arms/titration/logs/make_inputs.out -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/predict_arms/titration/logs/make_inputs.err
export PATH=/rhome/jstajich/projects/funannotate/funannotate-live/.pixi/envs/default/bin:$PATH
cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/predict_arms/titration && python make_inputs.py
