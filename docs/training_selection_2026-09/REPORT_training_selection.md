# PASA training-model selection: decisions and results

Session: do-annotation-triagepasarerun (Claude), started 2026-09-25.
The user approved using SLURM freely and asked for all decisions and results to be recorded here for a morning check-in.
PASA-side work (R1, R2, R13, the gffcompare scorer, diversity metrics) is done by the peer session `pasapipeline-df`.
Its record is `~/projects/funannotate/PASApipeline/CODE_REVIEW_20260925.md`.

## 1. Code state (funannotate-live, uncommitted)

| Change | Where | Test |
|---|---|---|
| RNA-seq concordance gate (`--min_rnaseq_map_rate`, exit 3) | `train.py`, `library.py` | `tests/test_training_gates.py` |
| PASA complete-ORF gate (`--min_pasa_complete_models`, BUSCO fallback) | `predict.py`, `library.py` | `tests/test_training_gates.py` |
| R3: train only on complete ORFs | `library.selectTrainingModels` | `tests/test_training_selection.py` |
| R5: one model per locus, ranked complete > CDS exons > CDS length > TPM | `train.getBestModel` | `tests/test_training_selection.py` |
| F12: transitive overlap removal | `library.selectTrainingModels` | `tests/test_training_selection.py` |

All test modules pass. `test_translation_table` prints the CLI help and exits 0; it behaved the same before these changes.

## 2. Decisions log

### D1 (2026-09-25): judge completeness from the protein, not from `cds_transcript`
- **Finding:** `gff2dict` stores `cds_transcript` in genome orientation for minus-strand genes. My first R3/R5 version tested that field, so every minus-strand model looked incomplete. In Meyerozyma that rejected 628 complete minus-strand models.
- **Fix:** `library.is_complete_model()` tests the translated protein instead: M start, `*` at the end, no internal `*`, codon_start 1.
  - It now finds 1,235 complete models in Meyerozyma (607 plus-strand, 628 minus-strand), which matches the independent coordinate-based counter exactly.
  - A minus-strand regression test was added; it failed on the buggy code and passes now.

### D2 (2026-09-25): the locus rule for `getBestModel` is decided by a RefSeq benchmark, not by argument
- **Old code:** a one-sided overlap on transcript spans, strand-blind. In Meyerozyma it discarded about 920 distinct complete-ORF loci.
- **Current code (tx_strand):** transcript spans, same strand. It keeps those loci, but still merges CDS-disjoint neighbors: 1,283 Meyerozyma loci, which loses 2,894 separate coding regions.
- **CDS-span rules** would keep those regions, but may also keep extra TransDecoder ORFs from one assembly. Those models reach EVM as PASA evidence at weight 6.
- **Benchmark:** `refseq_benchmark/benchmark.py` (job 29107015) compares old, tx_strand, cds_strand and cds_blind against RefSeq CDS chains for N. crassa OR74A, A. nidulans FGSC A4 and Botrytis cinerea B05.10.

### D3 (2026-09-25): predict-arm experiment design (user approved)
Question: does the training set from each selection variant give more accurate final predictions?
- **Same software stack in every arm.** All arms run in `funannotate-1.9.0-rc.1.sif`.
  - The image's funannotate code is identical to git HEAD `41a2fd7` (no diff under `funannotate/` between the image's `6f3eaad` and HEAD).
  - "New code" arms bind-mount a frozen copy of the working tree over the image's `site-packages/funannotate`, so only the Python code differs between arms.
- **Chromosome holdout** (peer's `predict_scorer.py split`): chromosomes of at least 1 Mb, sorted by length, alternated into train and holdout; the mitochondrion is excluded.
- **Step A (training):** `funannotate predict` on a genome made of the train chromosomes only, with that arm's PASA models restricted to the train chromosomes. It produces `parameters.json` (Augustus + SNAP).
- **Step B (prediction):** `funannotate predict` on the full genome with `-p parameters.json`. Every arm gets identical evidence: the same PASA GFF3 (the old `getBestModel` output), transcript alignments, GeneMark GTF and protein evidence. Step B always uses the image's own (old) code; training is not exercised in step B because parameters are pretrained.
- **Scoring:** gffcompare on held-out chromosomes only (peer's `predict_scorer.py score`, CDS level). Also BUSCO (fungi_odb10, protein mode) on the full predicted proteome.
- **GeneMark:** self-trained once per genome (ES mode, braker3 image, the same as the pipeline's `GENEMARK_RUN`) on the full genome. It is the same for all arms. ES training uses no labels, so this is not holdout leakage.
- **Transcript BAM:** the rc.1 train BAMs were not kept. The Trinity-GG FASTA is re-aligned with minimap2 splice mode for each genome, and separately to the train-chromosome genome for step A. Identical across arms.
- Arms, EVM settings and any later changes are recorded below.

### D3a (2026-09-25): predict-arm run details
- Arm script: `predict_arms/arm.sh`; submitter: `predict_arms/submit_arms.sh`; job ids: `predict_arms/jobs.tsv` (45 jobs).
- Predict flags mirror the pipeline's FUNANNOTATE_PREDICT: `--busco_db dikarya -w codingquarry:0 glimmerhmm:0 genemark:1 --min_training_models 30 --keep_no_stops --header_length 24 --protein_evidence lib/swissprot_fungi.faa --max_intronlen 3000 --min_intronlen 10 --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark --genemark_gtf ...`, plus `--rna_bam`, `--pasa_gff` and `--transcript_alignments`.
- Each arm gets its own copy of the Augustus config directory, so concurrent arms training the same species name cannot overwrite each other.
- **New-code PASA arms (oldnew, tx, cdsS, cdsB) run with `--min_pasa_complete_models 0`.** Reason: on half a genome, an arm could fall below the 500-model gate and quietly switch to BUSCO training, which would mix two experiments. The busco arm forces BUSCO training with a threshold of 1e9.
- The old arm runs the image's own code (git HEAD 41a2fd7 equivalent), so it has no gate and uses the old selection.
- Step B always runs the image's own code, because parameters are pretrained with `-p`.
- Scores: each step B job writes `predict_arms/<genome>/<arm>.B.<evid>/score.tsv` (peer scorer) and `busco_short_summary.txt`.

### D3b (2026-09-25 22:45): harness fix and resubmission
- First submission: step A for old and oldnew on N. crassa crashed at Augustus training with `UnboundLocalError: AUGUSTUS_BASE` (both old and new code).
- **Cause:** `predict.py` sets `AUGUSTUS_BASE` only when the base name of `--AUGUSTUS_CONFIG_PATH` is exactly `config`. My per-arm copy was named `augustus_config`.
- **Fix:** each arm now copies the config to `<arm>/augustus/config`. All 45 jobs were cancelled and resubmitted (new ids in `predict_arms/jobs.tsv`).
- **Side note (funannotate robustness, not fixed):** a config directory with any other name crashes Augustus training instead of giving a clear error.

## 3. Results

### 3.1 RefSeq benchmark (job 29107015; `refseq_benchmark/benchmark.tsv`; DECISIONS.md D23)
Reads: N. crassa = divergent strain HJDF (about 95-98% identity, REVIEW D19); A. nidulans and Botrytis = same strain (D20). Reported per genome only.

**Training set:** new selection (R3 + R5 + F12) against the old pipeline.

| Genome | Old: models, exact % | New: models, exact % | Exact RefSeq genes, old → new | Redundant, old → new |
|---|---|---|---|---|
| N. crassa | 2,889, 38.9% | 1,515, 78.7% | 1,124 → 1,192 | 46 → 0 |
| A. nidulans | 5,007, 50.1% | 4,123, 62.1% | 2,508 → 2,559 | 37 → 4 |
| Botrytis | 5,637, 75.5% | 4,866, 89.4% | 4,255 → 4,348 | 59 → 0 |

- The new selection removes mostly wrong models and keeps every correct one.
- The three locus rules give an identical training set.

**EVM evidence (getBestModel output):**
- tx_strand is better than old on every genome.
- cds_blind dominates cds_strand.
- tx_strand against cds_blind is a trade-off. On Botrytis, cds_blind gives +268 exact genes but also +263 models with no RefSeq overlap and +100 redundant (precision 64.8% → 62.5%).
- The decision waits for the predict arms that use each variant's own evidence (B.own).

**Side finding:** every training set (old and new code) has 0% single-exon genes, because of the multi-CDS requirement. RefSeq has about 21% single-CDS mRNAs in N. crassa. This is review item R6. **It needs a user decision.**

**Correction (see DECISIONS.md, the entry after D27):** the benchmark ran selection without the filterGeneMark keeper set, so its training-set gains are an upper bound. In real predict step A on A. nidulans, the old and new selection give 1,103 training models each (12 CDS lines differ), because filterGeneMark (1,377 keepers) already removes most incomplete models. R3 changes the result only when there are fewer than 200 keepers: 95 of 400 sampled production PASA-trained genomes (24%).

(Predict-arm results follow when jobs finish.)

### 3.2 Status update (2026-09-26 early morning)
- **User decision (D35):** no commit until the predict arms settle the locus rule (tx_strand vs cds_blind). Then apply the review fixes and ask the user again.
- **Fable 5.1 review of gates + R3/R5** (`REVIEW_R3R5_fable.md`, D41).
  - Fixed test-first in the working tree:
    - the selectTrainingModels exit when R3 leaves 0 models;
    - `is_complete_model` checks CDS length % 3;
    - the predict gate uses `lib.pasa_gate_applies` (any predictor in pasa mode; skipped on resume or with `--augustus_gff`).
  - Being decided by data: finding 1 (getBestModel ranks complete before structure, so a short complete single-exon ORF beats a long 5'-partial multi-exon one). Job 29107681 → `refseq_benchmark/rank_benchmark.tsv`, comparing complete_first / struct_first / guarded on the rc1, R1 and R1R2 PASA sources.
- **Real-world scope of R3 (D29):** in real predict, filterGeneMark already cleans the training set when there are ≥200 keepers (A. nidulans old = new, 1,103 models). R3 matters on low-keeper genomes: 24% of sampled production PASA-trained genomes.
- **Low-keeper production genomes (D30):** A. niger CBS 101883 (992 genes), A. pullulans (2,817 genes) and others collapsed despite reads mapping ≥92%. Read identity (`qc_training_sweep_20260925/identity_check/identity.out`): A. niger 100.0% and C. neoformans H99 100.0% (same strain); A. pullulans 98.0% and P. antarcticum 92.7% (divergent). So same-strain collapses need another explanation; REVIEW is looking at their PASA logs. H99 and S. commune were added to all benchmarks and arms.
- **R6 option (b)** (single-exon training genes, D38/D42/D45): implemented opt-in as `predict --training_single_exon`.
  - The share estimate uses GeneMark-ES (error +1.3 to +6.5 points vs RefSeq) instead of protein alignments (−2.5 to −13.6).
  - **The source choice needs a user decision (D45).**
  - Arms se vs tx2 on all five genomes (`predict_arms/jobs_se.tsv`).
- Full funannotate test suite: all modules pass (gates 21, selection 25+). `test_translation_table` prints help (behavior predates these changes).
- **F1 production scope (DECISIONS D47):** 5,599 production PASA sets (dated 2026-06-29 to 09-24) carry F1 frame errors, plus 3,052 genomes inherit parameters from one. Rerun candidate list v2 has 9,374 genomes (`qc_training_sweep_20260925/rerun_candidates_v2.tsv`). Sets from before 06-29 show 0 frame errors.

### 3.3 Interim predict-arm results (holdout chromosomes, gffcompare CDS level; `predict_arms/scorecard.tsv`)
Single genomes, no replicates: differences of 1-2 points show direction, not statistical significance.

**Locus rule (step B with each variant's own PASA evidence):**

| Genome | tx_strand: locus Sn/Pr, intron-chain Pr | cds_strand | cds_blind |
|---|---|---|---|
| N. crassa | 63.8/71.8, 72.1 | 63.5/70.9, 71.6 | 63.5/70.8, 71.6 |
| A. nidulans | 54.7/55.7, 54.2 | 54.7/55.1, 54.0 | pending (SIGILL rerun) |

The CDS-span rules add evidence models with no sensitivity gain and about 0.6-1.0 points lower locus precision. **Interim: keep tx_strand.**

**Training selection (step B, fixed evidence = old getBestModel output):**

| Genome | old | oldnew | tx | busco (forced) |
|---|---|---|---|---|
| N. crassa locus Sn/Pr | 63.3/71.6 | 63.5/71.8 | 64.3/72.6 | **66.0/74.1** |
| A. nidulans locus Sn/Pr | 54.6/55.7 | 54.7/55.8 | 54.6/55.8 | pending |
| Botrytis locus Sn/Pr | pending (rerun) | pending (rerun) | **83.3/82.2** | 82.0/82.0 |

- The new selection helps N. crassa by about 1 point and is neutral on A. nidulans, as expected from D29 (filterGeneMark already cleans A. nidulans).
- BUSCO training beats PASA training on N. crassa (divergent reads, 96.7% identity) but loses on Botrytis (same strain). This is consistent with the user's D37 identity gate (divergent reads → adapt or avoid PASA training).
- **Locus rule settled on 3 genomes (DECISIONS D54): keep tx_strand.** With each rule's own evidence, the CDS-span rules add about 100 mostly single-exon predictions, gain ≤0.3 locus Sn and lose 0.6-1.2 locus Pr (Botrytis tx 82.9/81.4 vs cds 83.2/80.2-80.3). Low-keeper genomes pending.
- **getBestModel ranking changed to guarded complete-first (D55)**, from the RefSeq ranking benchmark: exact matches within 0-20 of complete-first, and 117-235 more correct intron chains per genome. Confirmation arms txG are running.
- **"BUSCO beats PASA on divergent reads" is being re-tested with fixed PASA (D53)**: arms txR1, txR1R2 and txID90 (N. crassa) and txR1 (Botrytis). The interim result used the broken minimap2 parser.
- **Commit status:** the locus rule is answered (D54). Before the commit, the user will be asked, with the txG result and the Fable fixes (D41) in hand.
- **Correction to §3.3 (BUSCO vs PASA training on same-strain genomes):** A. nidulans busco-forced is locus 54.8/56.8 and intron chain 50.1/55.4, against PASA-trained tx 54.6/55.8 and 51.0/54.1. BUSCO has about 1 point higher precision; on Botrytis, PASA was ahead (83.3/82.2 vs 82.0/82.0). **So there is no consistent winner on same-strain genomes.** The differences are about 1 point, on single genomes. The only clear BUSCO lead (N. crassa, 66.0/74.1 vs 64.3/72.6) used old-parser PASA and is being re-tested with fixed PASA (D53).
- Botrytis, fixed evidence: old 83.3/82.2, oldnew 83.1/82.0, tx 83.3/82.2. Selection is neutral where filterGeneMark keeps enough models (D29).

## 4. Final summary for the user (2026-09-26 morning)
All predict arms have finished. The full per-arm table is `predict_arms/scorecard.tsv`; single/multi-exon exact matches are in `predict_arms/single_exon_scores.tsv`. Decisions are in `DECISIONS.md` D07-D79. Every number below is a holdout-chromosome gffcompare or exact-CDS score against RefSeq, from one genome each (no replicates).

**Settled by data:**
1. **Locus rule: tx_strand** (D54, D57). CDS-span rules add about 100 predictions, gain ≤0.3 Sn and lose 0.6-2.6 Pr on all 4 genomes where both were scored.
2. **Ranking: complete-first stays the default** (D69). The guarded rule is an option (`complete_min_frac=0.8`); it was prediction-neutral.
3. **R3 (complete ORFs only) must ship with the predict gate** (D78, D79). Alone it hurts a genuine low-keeper genome (S. commune 29.2 → 24.9 locus Sn). With the gate (93 complete < 500), that genome trains from BUSCO and reaches 37.4/49.2 (vs 29.2/44.2 on the old production path).
4. **F1 is the main cause of production collapses** (D47, D75, D76). H99 recovers fully with rc.1 PASA. Low-keeper genomes: 30% among F1-affected vs 10% among clean. Rerun list v2: 9,374 genomes.
5. **abu_dhabi nodes c01-c30 must be excluded** (D49). SIGILL in bam2hints/kallisto.

**Waiting for the user:**
- **D64/D77, R6 option (b) (single-exon training) as the default?**
  - 4 genomes: single-exon Sn +5 to +8; multi-exon Pr up in all; multi-exon Sn −0.9 to 0.0; net exact genes up.
  - On S. commune it admits nothing (no protein-supported single-exon models).
  - REVIEW recommends enabling it.
- **D73, action for the divergent category (90-99% identity):**
  - N. crassa: fixed PASA (R1+R2) + single-exon training 67.1/75.0 ≥ BUSCO + fixed PASA evidence 66.8/74.4 > fixed PASA alone 65.8/73.6.
  - S. commune (divergent, but low-keeper): BUSCO wins by about 8.
  - Suggestion: fixed PASA + (b) when the gate passes; BUSCO via the gate when there are too few complete models.
- **D35, commit.** Scope:
  - funannotate-live: gates, R3/R5 + tx_strand, complete-first with the optional guard, review fixes (D41), R6 (b) opt-in, R2/F4 flag passthrough (D63), tests, CHANGELOG;
  - Fungi_BFD/nextflow: gate and flag params, train log capture, exit-3 handling, predict_misc keep list.
  - Then REVIEW's R1 PR, then the PASA pin to v2.6.1-rc.2, then the image rebuild (D59).
