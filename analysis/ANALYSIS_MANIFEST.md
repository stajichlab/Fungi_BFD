# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### ideas/2026-06-18-genome-size-composition-framework

```yaml
name: ideas-genome-size-composition-framework
type: ideation-session
status: complete
created: 2026-06-18
last_updated: 2026-06-18
datasets: [tables/asm_stats.tsv.gz, tables/ (BFD.nf compositional tables)]
algorithms: []
parent_analysis: null
key_findings:
  - "9 personas (Broader-9) generated 18 grounded ideas for exploring fungal genome-size variation and its compositional drivers."
  - "2 low-effort/data-ready ideas: #8 order-parameter/phase-transition in TE expansion (asm_stats only), #17 Zipf-Mandelbrot codon/AA lexicon shifts."
  - "Within-taxonomy variation -> variance partitioning (#1/#13/#14); contributing indicators -> causal/scaling/diversity (#3/#4/#7/#9/#12); framework deliverable -> residualized PCA->UMAP (#12) + VAE (#11)."
  - "Common prerequisite: materialize the All_Taxa compositional merge (one BFD.nf run); universal confound: assembly quality + taxonomic non-independence."
report: analysis/ideas/2026-06-18-genome-size-composition-framework/00_index.md
tags: [ideation, genome-size, transposons, composition, comparative-genomics, framework]
```

Persona-based ideation session (idea-generator convention) over the BFD
comparative-genomics tables. Produced 18 ideas spanning evolutionary biology,
causal inference, information theory, statistical physics, ecology,
representation learning, quantitative genetics, topology, and linguistics —
indexed in `00_index.md`, grounded in `CONTEXT.md`. Aimed at understanding how
genome size varies within taxonomy, which indicators (TEs especially) contribute,
and how to build a reusable exploratory framework.

### rnaseq-trinity-qc

```yaml
name: rnaseq-trinity-qc
status: active
created: 2026-06-18
last_updated: 2026-06-18
datasets: [misc/TrinityGG_summary.tsv, misc/TrinityGG_failed.tsv, misc/trinity_lowcount_report.tsv, rnaseq_rerun_report.tsv, samples.rnaseq_sra.csv]
algorithms: []
parent_analysis: null
key_findings:
  - "239 of 1178 attempted species have a Trinity-GG FASTA; 939 are no-assembly failures (874 MISSING + 65 EMPTY)."
  - "Existing assemblies are overwhelmingly usable: 234 PASS / 3 BORDERLINE / 2 FAIL; tier headline robust to threshold choice."
  - "Failed report's 5 LOW_COUNT rows double-list successes (not separate failures); counting them as failures inflates the denominator."
  - "All 20 curated low-count species reconciled: 4 still successes (4/4 flagged sub-PASS), 16 now MISSING/EMPTY."
report: null
tags: [rnaseq, trinity, qc, annotation-evidence, funannotate, sensitivity-analysis]
```

QC classification of genome-guided Trinity transcript assemblies into usability
tiers for funannotate annotation evidence, plus a reconciled accounting of the
success/failure landscape. The dominant QC problem is *missing* assemblies
(~80% of attempted species), not poor-quality ones; among assemblies that exist,
all but five are strongly usable. A strict disjointness check surfaced that the
upstream failure report double-counts low-count successes — corrected here.

### funannotate_model_failures

```yaml
name: funannotate-model-failures
type: failure-screen + pipeline-guard
status: complete
created: 2026-06-24
last_updated: 2026-06-25
datasets: [genome_annotation/*/logfiles/funannotate-predict.log (8082), samples.csv, results/genome_stats/asm_stats/, input_clean_genomes/]
algorithms: [parse_predict_failures.py]
parent_analysis: null
key_findings:
  - "8,054/8,075 genomes completed funannotate predict; 9 failed with 'Not enough gene models N to train Augustus (30 required)'."
  - "All 9 too-few-models failures are small AND fragmented (127 kb–15.7 Mb); a small-AND-fragmented screen catches 9/9 with zero false-skips of complete small genomes (Malassezia, Ashbya, Rozella, Microsporidia verified safe)."
  - "6 INCOMPLETE_unknown are the opposite class: oversized genomes (0.5–1.3 Gb rusts; good 18–34 Mb assemblies) timing out — a resource problem, not size-detectable."
  - "results/genome_stats/asm_stats/ is stale (4,094 rows; none of the 21 failing genomes present) — pre-screen needs CALC_ASM_STATS re-run."
report: analysis/funannotate_model_failures/FUNANNOTATE_MODEL_FAILURES.md
tags: [funannotate, augustus, training-models, asm-stats, too-small, basidiomycota, suppress-txt, pipeline-guard]
```

Ground-truth screen of all persisted `funannotate-predict.log` files for the
too-small/fragmented "not enough training models" failure, joined to taxonomy
(unique `--name LOCUSTAG` key) and assembly metrics. Produced a per-genome outcome
table + crosstab validating that a small-AND-fragmented asm_stats rule cleanly
separates the unfixable too-small failures from the oversized-genome resource
failures. Drove a two-layer guard in `funannotate.nf` (pre-flight gate +
post-predict catch) and 7 additions to `suppress.txt`. Reproduce: `run.sh`.

### nextflow_memory_profile

```yaml
name: nextflow-memory-profile
type: pipeline-tuning
status: active
created: 2026-07-17
last_updated: 2026-07-23
datasets: [logs/nextflow/funannotate_trace*.txt (all roots), nextflow/conf/profile_funannotate.config]
algorithms: [scripts/profile_nextflow_memory.py]
parent_analysis: null
key_findings: [FUNANNOTATE_TRAIN attempt-1 cpus=2 is CPU-adequate (p90 ~105% of 2 cores); all 4 profiled processes are OOM_RISK by retry-escalation design, not under-provisioned]
report: analysis/nextflow_memory_profile/NEXTFLOW_MEMORY_PROFILE.md
tags: [nextflow, pipeline-tuning, memory, resource-allocation, trace-analysis, funannotate]
```

Profile Nextflow funannotate task memory usage (peak RSS) from all historical run
directories to identify under- and over-provisioned `withName` process blocks in
`nextflow/conf/profile_funannotate.config`. Compares observed peak RSS against
configured attempt-1 memory allocations to flag processes for retuning:
`OOM_RISK` (exit code 137/140 or RSS > 85% config), `OVERPROVISIONED`
(p90 RSS < 40% config), or `OK`. Deduplicates task hashes across duplicate trace
files. Reproduce: `bash analysis/nextflow_memory_profile/run.sh`.

### funannotate_train_stage_timing

```yaml
name: funannotate-train-stage-timing
type: pipeline-tuning
status: active
created: 2026-07-23
last_updated: 2026-07-23
datasets: [genome_annotation_training/*/logfiles/funannotate-train.log, genome_annotation/*/logfiles/funannotate-train.log (400/5471 random sample, seed=42)]
algorithms: [scripts/profile_train_stage_timing.py]
parent_analysis: null
key_findings: ["PASA (Launch_PASA_pipeline.pl) is median 84.8% (p90 94.4%) of a training run's total wall time; median PASA time 73.2 min of median 73.8 min total run", "Trinity-GG genome-guided assembly never runs in 336/336 PASA-dominant sampled logs (per-strain runs reuse shared Trinity via --trinity); only 4/400 runs use an older hisat2+StringTie path", "TransDecoder training-set extraction is a distant second at 6.4% of sampled wall-clock"]
report: analysis/funannotate_train_stage_timing/FUNANNOTATE_TRAIN_STAGE_TIMING.md
tags: [funannotate, pasa, trinity, throughput, log-parsing, trace-analysis, upstream-bugfix]
```

Parse `[MM/DD/YY HH:MM:SS]:`-timestamped `funannotate-train.log` lines to attribute
wall-clock time to named pipeline stages (Trinity-GG, seqclean, transcript alignment,
PASA, TransDecoder, Kallisto, best-model selection), using a sticky/stateful
"current stage" tracker so multi-hour gaps following an unmatched raw-subprocess
command-echo line still attribute correctly to the stage that launched them. Built to
decide whether the upstream Trinity-GG threading fix
([nextgenusfs/funannotate#1178](https://github.com/nextgenusfs/funannotate/issues/1178))
or PASA optimization work matters more for real throughput — answer: PASA, by a wide
margin. Reproduce: `bash analysis/funannotate_train_stage_timing/run.sh`.

### funannotate_predict_stage_timing

```yaml
name: funannotate-predict-stage-timing
type: pipeline-tuning
status: active
created: 2026-07-23
last_updated: 2026-07-23
datasets: [genome_annotation_training/*/logfiles/funannotate-predict.log, genome_annotation/*/logfiles/funannotate-predict.log (400/11320 random sample, seed=42)]
algorithms: [scripts/profile_predict_stage_timing.py]
parent_analysis: null
key_findings: ["Ab-initio training (GeneMark-ES + Augustus [PASA or BUSCO path] + SNAP) is median 39.8% (mean 38.7%) of a predict run's total wall time", "Augustus training cost is bimodal: PASA-path (RNA-seq available, 29% of sample) is negligible at median 0.3 min; BUSCO-path (no RNA-seq, 67% of sample -- the majority case) costs median 25.9 min (p90 61.2 min), comparable to GeneMark-ES", "GeneMark-ES (251.9 total sampled hours) and Augustus BUSCO-path training (143.8 total sampled hours) are comparable in aggregate cost, not GeneMark-ES-dominant as originally assumed", "A single hand-inspected log had suggested RNA-seq hints prep dominated (~72%); the 400-log aggregate shows this was an outlier -- median share is 0.0%, mean 9.5%"]
report: analysis/funannotate_predict_stage_timing/FUNANNOTATE_PREDICT_STAGE_TIMING.md
tags: [funannotate, genemark, augustus, snap, busco, throughput, log-parsing, trace-analysis, ab-initio-reuse]
```

Parse `[MM/DD/YY HH:MM:SS]:`-timestamped `funannotate-predict.log` lines to attribute
wall-clock time to named pipeline stages (RNA-seq hints prep, GeneMark-ES self-training,
Augustus training [PASA or BUSCO path], SNAP train+predict, EVM, tRNA/tbl2asn finishing),
using the same sticky/stateful stage-tracking method as
`analysis/funannotate_train_stage_timing/`. Built to test the cost premise behind
`todo/species_level_abinitio_reuse.md` (sharing trained ab-initio parameters across
ANI-qualified strains of the same species) — directly informed that plan's decision to
share all three predictors (not just GeneMark-ES) in the first rollout wave, since Augustus's
BUSCO-seeded training path turned out to be both expensive and the majority case in this
dataset. Reproduce: `bash analysis/funannotate_predict_stage_timing/run.sh`.

```yaml
name: funannotate-abinitio-reuse-gene-diff
type: qc-characterization
status: complete
created: 2026-07-25
last_updated: 2026-07-25
datasets: [genome_annotation/*/predict_results/*.gff3+proteins.fa (condition A, T-004), do_annotation/.claude/notrain_test/*/predict_results/*.gff3+proteins.fa (condition B, T-013)]
algorithms: [bedtools intersect -v]
parent_analysis: null
key_findings: ["A_only loci (300-634/strain) consistently ~10-20x larger than B_only loci (21-45/strain) across all 12 strains -- T-013's gene-count deficit vs T-004 is overwhelmingly one-directional (T-004/PASA calls genes T-013 misses), not scattered noise in both directions", "A_only count is larger for Aspergillus fumigatus (329-634) than Beauveria bassiana (291-355), consistent with F-006's species-level gene-count-delta asymmetry", "Protein sequences for both unique sets extracted per strain but not yet manually inspected -- open question is whether A_only genes are genuine PASA-rescued loci or over-calling artifacts"]
report: analysis/funannotate_abinitio_reuse_gene_diff/FUNANNOTATE_ABINITIO_REUSE_GENE_DIFF.md
tags: [funannotate, gff3, bedtools, gene-diff, ab-initio-reuse, T-004, T-013, pasa]
```

Coordinate-level diff (via `bedtools intersect -v` on `mRNA` features, both conditions
predicted off the same masked genome fasta so no liftover needed) between T-004
(ab-initio reuse + full TRAIN) and T-013 (no-TRAIN, `--transcript_evidence` + ab-initio
reuse) predict outputs for the 12 strains validated in
`.living/findings/funannotate-abinitio-reuse-validation.md` (F-006), to characterize
*which* genes account for T-013's consistent gene-count deficit rather than just its
magnitude. Produces per-strain unique-locus ID lists and protein FASTAs for direct
inspection. Reproduce: `bash analysis/funannotate_abinitio_reuse_gene_diff/run.sh`.

```yaml
name: pfam-hmmsearch-perf
type: performance-profiling
status: complete
created: 2026-08-01
last_updated: 2026-08-01
datasets: [real Pfam-A HMM DB (db-pfam module), Malassezia_brasiliensis_CBS_14135.proteins.fa (3786 proteins)]
algorithms: [hmmsearch A/B timing (shared vs scratch storage), MPI task-scaling A/B (pfam_tasks=1/2/4/8/16), minimal compiled MPI hello-world for --cpu-bind launch diagnosis]
parent_analysis: null
key_findings: ["Scratch-copying the Pfam-A DB to node-local NVMe gives no measurable wall-time benefit vs reading from shared NFS storage (~1% noise) -- the hmmsearch task is compute-bound, not I/O-bound", "hmmsearch only needs the raw Pfam-A.hmm file, not the hmmpress-ed .h3* index files (those are an hmmscan requirement)", "RUN_PFAM's existing MPI code path (srun -N/-n --mpi hmmsearch) is currently broken as written on this cluster -- fails with 'CPU binding outside of job step allocation' unless --cpu-bind=none is added to srun; dormant today only because pfam_tasks/pfam_nodes default to 1", "With --cpu-bind=none, real MPI task scaling is non-monotonic: pfam_tasks=2 is slower than no MPI (HMMER's rank 0 dispatches rather than searches, so -n 2 is really 1 worker plus overhead), pfam_tasks=4 peaks at +15% faster than baseline, pfam_tasks=8/16 degrade back toward and below baseline"]
report: analysis/pfam_hmmsearch_perf/PFAM_HMMSEARCH_PERF.md
tags: [pfam, hmmer, hmmsearch, mpi, cpu-bind, slurm, srun, scratch-partition, task-scaling, performance, T-015]
```

Real (non-stub) A/B timing investigation for T-015 (`todo/pfam_hmmer_performance.md`),
raised as an adjacent-but-independent question while touching `MERGE_PFAM` during the
T-014 genome_stats/function storage reorg. Answers both open questions with data:
scratch-copying the Pfam-A DB isn't worth implementing (compute-bound, not I/O-bound),
and enabling MPI parallelism for `hmmsearch` gives a real but modest (+15%) gain at
`pfam_tasks=4`, non-monotonic beyond that -- but first requires fixing a real,
previously-undiscovered bug in `RUN_PFAM`'s MPI launch (`--cpu-bind=none` missing).
The `--cpu-bind` finding is a general HPCC/MPI mechanism, not PFAM-specific, and is
also logged in `.living/learnings.md` and the personal `nextflow-hpcc` Claude Code
skill for reuse beyond this repo.

### telomere_finder

```yaml
name: telomere-finder
type: pipeline-feature
status: complete
created: 2026-08-03
last_updated: 2026-08-03
datasets: [genome_annotation/*/predict_results/*.scaffolds.fa, samples.csv]
algorithms: [nextflow/bin/find_telomeres.py, nextflow/bin/summarize_telomeres.py]
parent_analysis: null
key_findings:
  - "Implemented a flexible telomere finder supporting multiple fungal monomer patterns, exact and fuzzy matching, per-tract repeat count, total telomeric length, and ~500 bp inward flanking sequence."
  - "Integrated into BFD.nf: FIND_TELOMERES per-genome step + MERGE_TELOMERES aggregate table (tables/telomeres.parquet)."
  - "Validated on Neurospora crassa OR74A: 2 terminal CCCTAA tracts, 342 bp total, 57 repeats across 2 scaffolds."
report: analysis/telomere_finder/TELOMERE_FINDER.md
tags: [telomeres, genome-statistics, nextflow, bfd-pipeline, fungi]
```

Pipeline feature adding telomere detection to the BFD genome-statistics stage.
Replaces the single-pattern reference script with a configurable multi-pattern
finder (regex + fuzzy modes, IUPAC support, canonical and --both-ends
orientation), per-genome TSV outputs with sequences, and a merged per-genome
summary table reporting total telomeric length and repeat counts. Tested on a
Neurospora crassa OR74A assembly via local Nextflow execution. Reproduce:
`bash analysis/telomere_finder/run.sh`.

### genemark_es_contribution

```yaml
name: genemark-es-contribution
type: pipeline-validation
status: complete
created: 2026-08-12
last_updated: 2026-08-12
datasets: [genome_annotation_training/*/training, genome_annotation/*/logfiles/funannotate-predict.log, samples.csv, analysis/funannotate_predict_stage_timing/outputs/per_run_summary.csv]
algorithms: []
parent_analysis: null
key_findings:
  - "GeneMark-ES's contribution to the final gene set ranges from dominant (68% of genes, fragmented 9.7Mb/1145-contig draft) to substantial (14%) to negligible/net-neutral (0.4%), n=3 genomes (F-008)."
  - "Dominant case is the smallest/most fragmented assembly, consistent with sparser PASA/RNA-seq training evidence leaving GeneMark as the largest ab-initio evidence source."
  - "Conclusion: GeneMark is not safely droppable for the rust-container migration; supports building a standalone GENEMARK_RUN Nextflow task rather than relying on --auto-skip-genemark degradation."
report: analysis/genemark_es_contribution/GENEMARK_ES_CONTRIBUTION.md
tags: [funannotate, genemark, evm, predict, ab-test, pipeline-validation]
```

Deterministic paired A/B rerun of `funannotate predict` (baseline vs.
`-w genemark:0`) on 3 already-trained genomes, against the same existing PASA
training data, to quantify GeneMark-ES's contribution to the final EVM
consensus gene set. Motivated by evaluating whether to split GeneMark into a
standalone Nextflow task (host module, feeding `--genemark_gtf` into a
container-based predict) as part of adopting the rust-optimized funannotate
container. Reproduce: `bash analysis/genemark_es_contribution/run.sh`, then
`python3 analysis/genemark_es_contribution/scripts/compare_results.py`.

### genemark_run_validation

```yaml
name: genemark-run-validation
type: pipeline-validation
status: complete
created: 2026-08-12
last_updated: 2026-08-12
datasets: [analysis/genemark_es_contribution/outputs/predict_runs/Penicillium_citrinum_NRRL_1841__baseline/predict_misc/ab_initio_parameters, input_clean_genomes]
algorithms: []
parent_analysis: genemark_es_contribution
key_findings:
  - "GENEMARK_RUN fresh --ES: 11,116 gene models, matching the A/B test baseline's internal GeneMark call (11,112) almost exactly."
  - "GENEMARK_RUN fast --predict_with reuse: completed in ~7 min at 2 cores vs fresh-ES's ~15.5 min at 8 cores -- confirms genuine training-free reuse, not seeded retrain."
  - "Real funannotate predict --genemark_gtf consumption: confirmed bypasses predict's internal GeneMark call entirely (resolves design doc's open question); final gene count 11,202 vs baseline 11,198, normal noise."
report: analysis/genemark_run_validation/GENEMARK_RUN_VALIDATION.md
tags: [funannotate, genemark, evm, predict, container-migration, validation, gmes_petap]
```

Real (non-stub) end-to-end validation of the `GENEMARK_RUN` Nextflow module
(nextflow/docs/GENEMARK_RUN_DESIGN.md) via a throwaway standalone workflow
that includes the real module file directly. All 3 tests (fresh ES, fast
reuse, predict consumption) pass. Reproduce:
`nextflow run analysis/genemark_run_validation/scripts/test_genemark_run.nf ...`
(see report for exact invocations).
