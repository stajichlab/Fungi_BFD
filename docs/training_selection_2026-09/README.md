# Training-data selection study (2026-09)

Copied 2026-09-26 from `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/` (the working directory; genome, BAM and predict output files stay there).

- `DECISIONS.md`: shared decision log of the two Claude sessions (D01-D86), with evidence and file paths.
- `REPORT_training_selection.md`: detailed report; section 4 is the summary.
- `training_data_selection_methods.md`: methods write-up for the funannotate paper (also in funannotate-live docs/, commit 66c75ea). `methods_tables.py` regenerates tables M1-M8 from the result files.
- `predict_arms/`: predict-arm scripts (`arm.sh`, `prep.sh`, `submit_arms.sh`, `scorecard.py`, `single_exon_score.py`), results (`scorecard.tsv`, `single_exon_scores.tsv`), job tables, and the diff of each frozen code snapshot against the base commit.
- `predict_arms/titration/`: experiment A (threshold calibration) scripts, task list and pool sizes. `titration_scores.tsv` holds one row per task (it replaces the per-task `rows/` files) and `titration_summary.tsv` holds mean/SD per genome x N. Both are built by `aggregate.py`.
  - State at the 2026-09-26 copy: N. crassa, A. nidulans and Botrytis are complete. The 34 H99 tasks failed on a missing container bind (DECISIONS D94) and were rerunning. The BUSCO-forced comparator repeats (N=busco, draws 1-3, D96) were still running. Refresh this folder when both finish.
- `refseq_benchmark/`: RefSeq benchmark scripts and results (`predict_scorer.py` and `diversity.py` were written by the PASA review session).

Related commits: funannotate-live 66c75ea (selection code), 8422006 (R1 bam2gff3 fix), Fungi_BFD c5ce230 (Nextflow).
