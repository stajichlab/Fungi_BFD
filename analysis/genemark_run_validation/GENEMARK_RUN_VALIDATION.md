# GENEMARK_RUN real end-to-end validation

## Motivation

`nextflow/docs/GENEMARK_RUN_DESIGN.md`'s Fable-reviewed design for a
standalone `GENEMARK_RUN` Nextflow process was implemented and validated at
the DAG/channel-wiring level (`-stub-run`), but not at the actual-execution
level — no real `gmes_petap.pl` had run, and predict's real handling of
`--genemark_gtf` was an open question read from source, not confirmed. This
analysis closes that gap with real (non-stub) executions against a real
genome.

## Method

Standalone throwaway Nextflow workflow
(`scripts/test_genemark_run.nf`) that `include`s the real
`nextflow/modules/funannotate/predict/GENEMARK_RUN/main.nf` directly (not a
duplicate/hand-copy — same file the production subworkflow uses), run three
times against `Penicillium_citrinum_NRRL_1841` (the same genome used in
`analysis/genemark_es_contribution/`'s F-008 A/B test, so results are
directly comparable to an already-established baseline):

1. Fresh `--ES` self-training (`shared_mod=''`)
2. Fast `--predict_with <mod>` reuse, using a real `.genemark.mod` already on
   disk from the A/B test's `predict_misc/ab_initio_parameters/`
3. Real `funannotate predict --genemark_gtf <gtf from step 1>`
   (`scripts/test3_predict_consumption.sbatch`), mirroring
   `FUNANNOTATE_PREDICT/main.nf`'s actual command construction against the
   same genome's real PASA training data

## Results

| Test | Result | Wall time | Pass criteria |
|---|---|---|---|
| Fresh `--ES` | 11,116 gene models | ~15.5 min (8 cores) | Comparable to A/B baseline's `GeneMark: 11112` — **pass** |
| Fast `--predict_with` | 100,966-line GTF (near-identical to fresh-ES's 101,000) | ~7 min (2 cores) | Genuinely faster than fresh-ES despite fewer cores, confirming true training-free reuse — **pass** |
| `predict --genemark_gtf` consumption | 11,202 final gene models; predict's internal GeneMark-ES call never invoked (confirmed via log: GeneMark was "functional" but unused; `Summary of gene models` shows `GeneMark 1 11116`, the exact count from test 1) | ~23 min | Matches A/B baseline's 11,198 (4-gene difference, normal tie-breaking noise) — **pass** |

**Design doc open question resolved**: `--genemark_gtf` does bypass predict's
own `GENEMARK_PATH` preflight/internal call entirely, confirmed directly
(not inferred from source reading).

**One test-harness-only bug found and fixed**: Nextflow's CLI parser treats
`--flag ''` (empty string) as a bare boolean flag, corrupting the intended
empty value. Fixed by omitting the flag rather than passing an empty string;
confirmed this cannot affect the real pipeline (which passes these values as
Groovy string literals inside `tuple()`, never through CLI parsing).

## ET mode (T-022) — standalone recipe validation (`et_eval/`, `et_eval2/`)

Separate from the 3 tests above (which cover the wired `GENEMARK_RUN` module
in ES mode only): evaluated whether GeneMark-ET is viable at all for this
pipeline, using real standalone commands (not yet wired into the module).

**First attempt failed**: `bam2hints --source=b2h` on
`training/transcript.alignments.bam`, filtered/rescored to match
`RunGeneMarkET()`'s exact convention → `gmes_petap.pl --ET` died in
`bp_seq_select.pl: hash is empty`.

**Root cause** (confirmed against `Gaius-Augustus/BRAKER`'s real pipeline,
not guessed): the hints had no strand assigned (`.` in column 7). BRAKER
always runs raw intron hints through `filterIntronsFindStrand.pl` before
GeneMark ever sees them — checks genome sequence at each intron boundary for
canonical splice sites (GT-AG/GC-AG/AT-AC), assigns strand, and drops
non-canonical introns. Branch-point signal is inherently strand-specific;
unstranded hints gave `bp_seq_select.pl` nothing usable.

**Fixed recipe** (`et_eval2/`), verified with a real `gmes_petap.pl --ET`
run: `bam2hints --intronsonly` → `filterIntronsFindStrand.pl` (vendored at
`nextflow/bin/vendor/filterIntronsFindStrand.pl`, Artistic License, from
BRAKER) → `join_mult_hints.pl` (AUGUSTUS script, already on PATH) →
`gmes_petap.pl --ET`. 20,603 of 21,400 hints (96%) passed stranding.
**Result: 10,776 gene models**, same order of magnitude as the ES run's
11,116. **T-022 status: unblocked, recipe validated, not yet wired into
`GENEMARK_RUN/main.nf` or the Nextflow subworkflow.**

Full trace, including the BRAKER source investigation and the ruled-out
`--bp_region_length` hypothesis: `nextflow/docs/GENEMARK_RUN_DESIGN.md`'s
"ET mode" section.

## Also found while validating: a latent EVM-weight bug in `FUNANNOTATE_PREDICT`

`predict.py:567` unconditionally zeroes GeneMark's EVM weight
(`StartWeights["genemark"] = 0`) whenever `gmes_petap.pl` isn't found on the
host running `predict` — with no check for whether `--genemark_gtf` was
supplied as an alternative evidence source. Not triggered by Test 3 above
(that host had `gmes_petap.pl` present), but would silently zero out a
correctly-supplied `--genemark_gtf`'s contribution once `FUNANNOTATE_PREDICT`
actually moves onto the rust container (no `gmes_petap.pl` there at all).
Fixed in `FUNANNOTATE_PREDICT/main.nf`: explicit `-w genemark:1` whenever
`genemark_gtf` is non-empty, merged into the *same* `-w` argument group as
the existing `codingquarry:0 glimmerhmm:0` (funannotate's `-w`/`--weights`
argparse option has no `action='append'`, so a second `-w` flag silently
replaces the first — confirmed directly against argparse — this was almost a
second bug introduced while fixing the first one).

## Status

Complete for ES mode (wired, tested, shipped) and for validating the ET
recipe (proven correct standalone, not yet wired). Full detail and
provenance in `.living/learnings.md` (2026-08-12 entries) and
`nextflow/docs/GENEMARK_RUN_DESIGN.md`.

## Reproduce

```bash
module load nextflow
cd analysis/genemark_run_validation

# Test 1: fresh ES
nextflow run scripts/test_genemark_run.nf -w work -c scripts/test_local.config \
    --genome_fa <path to a masked genome.fa.gz> --out <name> \
    --asmid <asmid> --species "..." --strain "..." --transl_table 1 \
    --force_independent false

# Test 2: fast reuse (needs a real shared .mod)
nextflow run scripts/test_genemark_run.nf -w work2 -c scripts/test_local_light.config \
    --genome_fa <same> --out <name>_reuse --asmid <asmid> \
    --species "..." --strain "..." --transl_table 1 \
    --shared_mod <path to a real .genemark.mod> --force_independent false

# Test 3: predict consumption (needs test 1's GTF at outputs/<out>.genemark.gtf)
sbatch scripts/test3_predict_consumption.sbatch
```
