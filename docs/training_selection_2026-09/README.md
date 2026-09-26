# Training-data selection study (2026-09)

Copied 2026-09-26 from `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/` (the working directory; genome, BAM and predict output files stay there).

- `DECISIONS.md`: shared decision log of the two Claude sessions (D01-D86), with evidence and file paths.
- `REPORT_training_selection.md`: detailed report; section 4 is the summary.
- `training_data_selection_methods.md`: methods write-up for the funannotate paper (also in funannotate-live docs/, commit 66c75ea). `methods_tables.py` regenerates tables M1-M8 from the result files.
- `predict_arms/`: predict-arm scripts (`arm.sh`, `prep.sh`, `submit_arms.sh`, `scorecard.py`, `single_exon_score.py`), results (`scorecard.tsv`, `single_exon_scores.tsv`), job tables, and the diff of each frozen code snapshot against the base commit.
- `predict_arms/titration/`: experiment A (threshold calibration) scripts, task list and pool sizes. `titration_scores.tsv` holds one row per task (it replaces the per-task `rows/` files) and `titration_summary.tsv` holds mean/SD per genome x N. Both are built by `aggregate.py`.
  - Complete (2026-09-26): 149 of 149 rows with status ok. The H99 tasks were rerun after a container bind fix (DECISIONS D94). N=busco draws 1-3 are BUSCO-forced repeats with the titration code (D96); run-to-run noise is at most 0.1 point (D100). N=busco_code_new is the earlier single run. The final threshold analysis is REVIEW's (see evidence_and_alignment_methods.md and the funannotate docs page assessment_pasa2.6_fun1.9).
- `refseq_benchmark/`: RefSeq benchmark scripts and results (`predict_scorer.py` and `diversity.py` were written by the PASA review session).

Related commits: funannotate-live 66c75ea (selection code), 8422006 (R1 bam2gff3 fix), Fungi_BFD c5ce230 (Nextflow).
