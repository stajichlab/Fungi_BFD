# predict_arms cleanup (2026-09-26)

This folder held the predict-arm and titration (experiment A) runs for the
funannotate training-selection work. On 2026-09-26 it was reduced to the records
needed to audit the results or rerun an arm. The KEEP/REMOVE rules were agreed
with the REVIEW session. They are in `cleanup_predict_arms.py` (run by
`cleanup_predict_arms.sh`, SLURM job 29114458).

## Record

- `cleanup_record/plan_keep.tsv`, `plan_remove.tsv`: every path with its size before the cleanup.
- `cleanup_record/manifest_kept.tsv`: every kept file after compression, with bytes and sha256.
- `cleanup_record/summary.txt`: totals.
- `cleanup_record_dryrun_29114339/`: the dry run before the last rule change (stock AUGUSTUS config under ab_initio_parameters).

## What was kept

- Run records: `arm_info.txt`, `run_status.tsv`, `score.tsv`, `busco_short_summary.txt`, `busco.log`, capture logs.
- funannotate `logfiles/` (predict log, EVM log, augustus logs, `predict_training_gate.tsv`).
  These runs used frozen code snapshots that predate `training_decisions.tsv`, so none have that file.
- Evidence and training inputs: arm `pasa.gff3`, titration `pasa.subset.gff3` and `n_subset.txt`, `pasa_inputs/`,
  `trinity.genome*.bam(.bai)`, `split.*_chroms.txt`, `genemark.genome*.gtf`, `transcripts.genome*.gff3`, `titration/inputs/`.
- `predict_misc/`: `weights.evm.txt`, `final_training_models.gff3`, `pasa_predictions.gff3`,
  `pasa.training.tmp.f.good.gtf` / `f.bad.gtf`, `gene_predictions.gff3`, `protein_alignments.gff3`,
  and `ab_initio_parameters/` (trained AUGUSTUS species dir, BUSCO_* training dirs, snap `.hmm`).
- `predict_results/`: final `.gff3`, `.proteins.fa`, `parameters.json`, `stats.json`, validation and error reports.
- All scripts, job tables, `scorecard.tsv`, `single_exon_scores.tsv`, titration tables and analysis outputs.

Plain-text files larger than 1 MB were compressed with `zstd -19` (suffix `.zst`).

## What was removed

- tbl2asn work dirs, EVM partitions and other `predict_misc` intermediates (hints, softmasked genome, p2g, repeats, tRNA).
- `predict_results` `.gbk`, `.scaffolds.fa`, `*-transcripts.fa` and `.tbl` (these can be rebuilt from the GFF3 and the genome).
- `genome.fa` / `genome_train.fa` copies. `prep.sh` rebuilds `genome_train.fa` from the production genome and `split.train_chroms.txt`.
- The per-arm AUGUSTUS config copies (`augustus/`), and the stock config under `ab_initio_parameters/augustus`
  (model, extrinsic, profile, parameters, species generic/anidulans). These are identical to `lib/augustus/3.5/config`
  and to funannotate_db `trained_species/anidulans`.
- `gmes/`, `busco_downloads/`.
- `code_new`, `code_new2`, `code_new3` and `code_new4` snapshots. Each is funannotate commit
  41a2fd799233b7d66b5d9923804ea98ff2ee3b20 (on GitHub, an ancestor of `target_1.9/rust_EVM_trinity_PASA`)
  plus the diff in `Fungi_BFD/docs/training_selection_2026-09/predict_arms/code_snapshots/`. Each diff was
  applied to that commit and reproduced every tracked file of its snapshot. The build-time `_version.txt`
  files are archived next to the diffs.

## Using the kept files

- `scorecard.py` and `single_exon_score.py` read `.gff3` or `.gff3.zst`.
- `refseq_benchmark/predict_scorer.py` runs gffcompare, which needs a plain GFF3. Decompress into node-local
  scratch first, for example: `zstd -dc out/predict_results/X.gff3.zst > $SCRATCH/X.gff3`.
- To rerun an arm, first rebuild `genome.fa` / `genome_train.fa` with `prep.sh` (or link the production genome),
  and apply the code snapshot diff to 41a2fd7.
