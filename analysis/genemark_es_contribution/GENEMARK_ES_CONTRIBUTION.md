# GeneMark-ES contribution to funannotate predict

## Motivation

While validating the `ghcr.io/nextgenusfs/funannotate` rust container for
Trinity/PASA (GeneMark is absent from that image by design — license
restrictions forbid redistribution), the question came up: how much does
GeneMark actually change the final EVM-consensus gene set today? Production
`FUNANNOTATE_PREDICT` (`nextflow/modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf`)
runs GeneMark in **ES mode** (self-training, no RNA-seq hints — confirmed by
reading `funannotate-live/funannotate/predict.py`: `--genemark_mode` defaults
to `ES` and production never passes `--rna_bam`/`--genemark_mode ET`). This
motivates a possible future split of GeneMark into a standalone Nextflow task
(`GENEMARK_RUN`, host module — GeneMark can't be containerized) feeding
`--genemark_gtf` into a container-based predict. Before investing in that,
this analysis measures whether GeneMark is worth keeping at all, and by how
much.

## Method

Deterministic paired rerun, not a sampled measurement: for 3 already-trained
genomes, rerun `funannotate predict` twice against the *same* existing PASA
training data (`genome_annotation_training/<out>/training/`) — once with
production weights (baseline) and once with `-w genemark:0` added
(nogenemark) — and diff the results. No bootstrap/permutation testing
applies here; the reportable quantity is an exact set difference between two
runs of the same deterministic pipeline on the same inputs, not an estimate
with a confidence interval.

**Genomes** (picked from `analysis/funannotate_predict_stage_timing/outputs/per_run_summary.csv`
for `augustus_evidence_mode=pasa` — i.e. real RNA-seq training data — and low
`genemark_es_train_seconds`, so each predict rerun finishes in well under an
hour):

| Genome | ASMID | genemark_es_train_seconds (baseline) | total_wall_seconds (baseline) |
|---|---|---|---|
| Penicillium_citrinum_NRRL_1841 | GCA_020284165.1_ASM2028416v1 | 611 | 2418 |
| Saccharomyces_kudriavzevii_IFO10991 | GCA_000257085.1_..._v1.0 | 676 | 2570 |
| Kluyveromyces_marxianus_YG-4 | GCA_053539435.1_ASM5353943v1 | 698 | 2070 |

**Fidelity to production**: same `funannotate/dev-1.9` host module, same
`--protein_evidence lib/swissprot_fungi.faa`, `--SeqCenter NCBI`,
`--header_length 24`, `-w codingquarry:0 glimmerhmm:0` (production already
zeroes these two), `--auto-skip-genemark`, `--min_training_models 30`,
intron length bounds — all copied from the actual command line recorded in
each genome's existing `logfiles/funannotate-predict.log`. One deliberate
deviation: `--AUGUSTUS_CONFIG_PATH` points at a **private per-job copy**
(seeded from the shared `Fungi_BFD_runs/lib/augustus/3.5/config`, ~39M),
not the live shared directory — that directory is written to by real
concurrent production jobs (`do_annotation_asco` etc. were actively running
at analysis time), and `funannotate predict` trains a new species config
into `$AUGUSTUS_CONFIG_PATH/species/<name>/`, so sharing it here risked
either colliding with production or having baseline/nogenemark for the same
species race each other.

## Running it

```bash
cd analysis/genemark_es_contribution
./run.sh                       # submits 6 SLURM jobs (3 genomes x 2 modes)
squeue -u $USER -n 'gmk_ab_*'  # watch progress

# once all 6 finish (predict_results/*.gbk present in each outputs/predict_runs/*/):
module load funannotate/dev-1.9   # brings bedtools onto PATH (used by compare_results.py)
python3 scripts/compare_results.py
```

`compare_results.py` parses each run's `Summary of gene models: {...}` log
line (raw per-source model counts fed into EVM) and does a genomic-coordinate
`bedtools intersect` between baseline and nogenemark final gene sets — a
count-only diff would miss cases where EVM substitutes a GeneMark model for
what would otherwise be an Augustus/snap call at the *same* locus, since
totals could match while individual gene identities differ.

## Status

Complete. First submission (jobs 27419586-91) failed all 6 with
`UnboundLocalError: local variable 'AUGUSTUS_BASE' referenced before
assignment` after GeneMark-ES itself had already completed — root cause:
`funannotate predict` only assigns `AUGUSTUS_BASE` when
`os.path.basename($AUGUSTUS_CONFIG_PATH) == "config"` exactly; the private
per-job Augustus config copy was named `<predictdir>.augustus_config`, not
`.../config`. Fixed (`scripts/run_predict_variant.sbatch` now copies into
`<predictdir>.augustus_config/config`) and resubmitted (jobs 27420274-79);
all 6 COMPLETED. Logged as a learning in `.living/learnings.md` (2026-08-12).

## Results

| Genome | Baseline final genes | Nogenemark final genes | Only in baseline (GeneMark-attributable) | Only in nogenemark |
|---|---|---|---|---|
| Saccharomyces_kudriavzevii_IFO10991 (9.7Mb, 1,145 contigs) | 882 | 328 | 603 (68%) | 47 |
| Kluyveromyces_marxianus_YG-4 | 2,632 | 2,256 | 381 (14%) | 3 |
| Penicillium_citrinum_NRRL_1841 | 11,198 | 11,363 | 44 (0.4%) | 152 |

Full per-genome numbers (raw EVM-input model counts, EVM-consensus counts,
final counts) in `outputs/genemark_contribution_summary.csv`.

**GeneMark's contribution is highly genome-dependent** — from dominant
(*S. kudriavzevii*: two-thirds of the final annotation exists only because of
GeneMark) to negligible/net-neutral (*Penicillium*). The dominant case is the
smallest, most fragmented assembly of the three, consistent with sparser
PASA/RNA-seq training evidence leaving GeneMark's self-trained ab-initio
model as the largest single evidence source. The baseline rerun for
*S. kudriavzevii* reproduced 882 final genes vs. the original production
run's 878 — validates the A/B setup's fidelity to production.

**Conclusion**: GeneMark is not safely droppable for the container migration.
Build the standalone `GENEMARK_RUN` Nextflow task (host module, feeding
`--genemark_gtf` into a container-based predict) rather than accepting
`--auto-skip-genemark`'s graceful degradation as good enough — that
degradation silently costs the most on exactly the genomes (small,
fragmented, RNA-seq-poor) where GeneMark turned out to matter most.

Full writeup: `.living/findings/funannotate-genemark-contribution.md` (F-008).
n=3 — the *magnitude* pattern (fragmented → GeneMark-dominant) is plausible
but not established; worth revisiting with more genomes before treating it
as a rule.
