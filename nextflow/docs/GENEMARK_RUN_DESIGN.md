# Design: standalone `GENEMARK_RUN` Nextflow process

Status: IMPLEMENTED — ES and ET both wired and real end-to-end validated
(2026-08-12). Reviewed by Fable before
implementation — 8 concrete gaps found, all incorporated (full review kept in
`.living/decisions.md`, 2026-08-12 entry, for provenance). DAG wiring
validated via `-profile funannotate,test,test_funannotate -stub-run`: all 5
processes (`GENEMARK_RUN`, `GENEMARK_RUN_SIB`, `FUNANNOTATE_PREDICT`,
`FUNANNOTATE_PREDICT_SIB`, `BACKFILL_ABINITIO_PARAMS`) registered with no
duplicate-name errors, `withName:GENEMARK_RUN` config selector matched both
the plain and `_SIB`-aliased invocation, and all reached "Starting process >"
with zero exceptions in `.nextflow.log` -- confirms the full channel-wiring
chain (joins/maps/`sharedGenemarkModFor()`) is syntactically and structurally
sound. No actual `GENEMARK_RUN` task instance executed in this stub run
because the 2-genome test fixture never clears `GENOME_CLEAN_BATCH`
upstream -- a pre-existing gap, confirmed unrelated to this change via a
prior session's identical warning in `nextflow/session-ses_05fc.md`.

**Real end-to-end validation also complete** (`analysis/genemark_run_validation/`,
2026-08-12): real `gmes_petap.pl --ES` fresh training, real `--predict_with`
fast reuse, and real `funannotate predict --genemark_gtf` consumption all
tested against a real genome and all pass — see "Real end-to-end validation"
section below for full results. This module is now validated at both the
DAG-wiring level and the actual-execution level.

## Motivation

The rust-optimized `ghcr.io/nextgenusfs/funannotate` container (see
`.living/learnings.md`, 2026-08-12) has Trinity/PASA/EVM all working, but
**GeneMark is deliberately absent** (its license forbids redistribution in a
public image). To eventually move `FUNANNOTATE_TRAIN` and `FUNANNOTATE_PREDICT`
onto the container, GeneMark has to come out of `funannotate predict` and
become its own step, running on the host module (where the licensed GeneMark
install already lives), producing a GTF that gets handed to the
container-based predict via `--genemark_gtf`.

Before committing to this, `.living/findings/funannotate-genemark-contribution.md`
(F-008) validated that GeneMark is worth keeping — its contribution to the
final gene set ranges from negligible to a genome's *majority* of predicted
genes (n=3), so `--auto-skip-genemark`'s graceful degradation is not an
acceptable substitute for the fragmented/RNA-seq-poor genomes where GeneMark
turns out to matter most.

## What `funannotate predict` does with GeneMark today (as-is, verified by reading `predict.py`/`library.py`)

- Mode is always **ES** (self-training, no RNA-seq hints) — `--genemark_mode`
  defaults to `ES`, production never passes `--rna_bam`/`ET`. `args.rna_bam`
  *is* auto-populated from `training/funannotate_train.coordSorted.bam` when
  present (`predict.py:584-591`), but that only feeds Augustus hints/StringTie,
  never GeneMark's own mode selection.
- **Ab-initio reuse today does not skip retraining.** When a sibling strain
  reuses a species' shared parameters (`-p parameters.json`, populated by
  `BACKFILL_ABINITIO_PARAMS`), `RunGeneMarkES()` is called with
  `ini=<species>.genemark.mod` — but that only **seeds** the run
  (`gmes_petap.pl --ES --ini_mod <mod> ...`); GeneMark still runs a full
  self-training pass. This matches what we measured in
  `analysis/funannotate_predict_stage_timing/`: `genemark_es_train_seconds`
  is consistently 600-800s even for genomes whose species has a
  pre-backfilled shared store.
- `gmes_petap.pl` itself supports a genuinely training-free mode,
  **`--predict_with <mod>`**, mutually exclusive with `--ES`/`--ET`/`--EP`
  (verified: `gmes_petap.pl:2066-2077`, `count_modes != 1` is a hard error) —
  funannotate's wrapper never uses it. This is the fast path a standalone
  `GENEMARK_RUN` can exploit that today's embedded call cannot: true
  prediction-only reuse of an existing `.mod`, not seeded-retrain.
- Output consumed by `--genemark_gtf` is GeneMark's **raw native GTF**
  (`genemark.gtf`, produced directly by `gmes_petap.pl` in its working dir) —
  predict.py runs `genemark_gtf2gff3.pl` on it internally
  (`predict.py:537,1735-1736`). We hand predict the raw GTF, no conversion
  needed on our side.
- The trained model funannotate already extracts and backfills today is
  `predict_misc/ab_initio_parameters/<lower_out>.genemark.mod`, copied by
  `backfill_abinitio_params.py` into
  `gene_prediction_shared_abinitio/<species_tag>/<species_tag>.genemark.mod`.
  **This is the exact artifact `GENEMARK_RUN` needs to produce and consume** —
  the shared-store layout does not need to change.

## Proposed process: `GENEMARK_RUN`

One process, two modes selected by whether a usable shared `.mod` exists —
mirrors the existing `ABINITIO_REUSE_FLAG` / `shared_params_json` pattern in
`FUNANNOTATE_PREDICT/main.nf` exactly, not a new paradigm.

```
process GENEMARK_RUN {
    tag "$out"
    cpus   16
    memory '32 GB'
    time   '4h'

    input:
    tuple val(out), val(asmid), val(species), val(strain),
          val(genome_fa), val(transl_table), val(force_independent),
          val(shared_mod)          // '' if none available / force_independent=true

    output:
    tuple val(out), val(asmid), val(species), val(strain),
          val(genome_fa), val(transl_table),
          path("${out}.genemark.gtf"), emit: gtf
    // Keyed tuple, NOT a bare path -- Fable review point 1: a bare
    // `path(...), optional: true` output desyncs positional alignment with
    // `gtf` the moment a batch mixes fresh-train and fast-reuse rows (mod is
    // only emitted on the fresh-train path). Keying by (out, species) lets
    // downstream joins stay correct even when some rows in a channel never
    // emit this output at all.
    tuple val(out), val(species), path("${out}.genemark.mod"), optional: true, emit: mod

    script:
    """
    module load genemarkESET   # or GENEMARK_PATH set directly; host-only, never containerized
    ...
    if [ -n "${shared_mod}" ] && [ "${force_independent}" != "true" ]; then
        # fast path: true prediction-only reuse, no training
        gmes_petap.pl --predict_with "${shared_mod}" --sequence genome.fa \\
            --cores ${task.cpus} --fungus
    else
        # fresh ES self-training (representative, independent, or forced)
        gmes_petap.pl --ES --sequence genome.fa --cores ${task.cpus} \\
            --max_intron ${params.max_intronlen} --soft_mask 2000 --fungus
        cp output/gmhmm.mod ${out}.genemark.mod
    fi
    cp genemark.gtf ${out}.genemark.gtf
    """
}
```

`--gcode` handling for transl_table 6/26 and the `_genemark_supports_gcode()`
version check carry over unchanged from `RunGeneMarkES`/`RunGeneMarkET` — copy
that logic in, don't reinvent it.

## Wiring into `FUNANNOTATE_PREDICTION.nf`

Precise insertion points below (revised after Fable review — the first draft
was ambiguous about exactly which channel `GENEMARK_RUN` reads from, and that
ambiguity hides real bugs: running GENEMARK_RUN on unfiltered `branched.*`
would repeat GeneMark self-training on every -resume for genomes that are
already up to date, defeating the whole point of the migration).

**Representative + independent rows** (both predict-scope branches share this
code): `GENEMARK_RUN` must consume `rep_todo`/`indep_todo` — the
*already-filtered* channels (post `gbkResult`/`staleRnaseq`/`staleGenome`),
not `branched.representative`/`branched.independent` directly. `shared_mod`
is always `''` for these (representative is what gets shared; independents
were never eligible) → always the fresh-ES-training path. Concretely:

```groovy
def rep_genemark_in = rep_todo.map { out, a, sp, st, lt, bl, hl, tt, gfa, _spj ->
    tuple(out, a, sp, st, gfa, tt, 'false', '')   // force_independent, shared_mod
}
GENEMARK_RUN(rep_genemark_in)   // and again for indep_todo, same shape

def rep_with_gtf = rep_todo.join(GENEMARK_RUN.out.gtf.map { out, ...rest -> tuple(out, rest[-1]) })
    .map { out, a, sp, st, lt, bl, hl, tt, gfa, spj, gtf -> tuple(out, a, sp, st, lt, bl, hl, tt, gfa, spj, gtf) }
```

`join()` by `out` here is safe (unlike the species-keyed join the subworkflow's
comments warn about) — `out` is 1:1 per genome, never many-to-one, in both
`rep_todo` and `indep_todo` individually. Do the join **before** `rep_todo.mix(indep_todo)`,
not after — post-mix there's no `shared_mod`/mode signal left to tell the two
populations apart if anything downstream needed to.

**Eligible siblings**: `GENEMARK_RUN` must consume `sibling_predict_todo`
(after `readyRows`/`gated`/`availableSpeciesSet` have resolved availability),
**not** `branched.eligible_sibling`. Resolving `shared_mod` any earlier means
duplicating the availability-gating logic (`readyRows`/`blockedRows`/
`allowFallback`) a second time — exactly the "second parallel
reuse-eligibility mechanism" the first draft said to avoid but didn't
actually avoid. Add `sharedGenemarkModFor(species)` next to
`sharedParamsJsonFor()` in `utils.nf` (same existence/non-empty check,
pointed at `<species_tag>.genemark.mod`), called at the same point
`sharedParamsJsonFor(sp)` already gets called inside `readyRows`'s `.map`.

**`--predict_scope representative_only` mode**: `rep_todo` construction is
shared code before the `if/else` split, so wiring `GENEMARK_RUN` once there
covers both modes for free. Eligible-sibling wiring is explicitly a **no-op**
in this mode — siblings are never predicted here at all
(`branched.eligible_sibling` has no consumer in this branch today), so don't
build sibling-availability plumbing that only the `all`-scope branch needs.

**Module input change**: `FUNANNOTATE_PREDICT`/`FUNANNOTATE_PREDICT_SIB` gain
a new tuple element (`genemark_gtf`, alongside the existing
`shared_params_json`). Pass `--genemark_gtf "${genemark_gtf}"` when non-empty;
predict's own internal GeneMark call becomes fully unreachable once every
caller supplies this (still keep `--auto-skip-genemark` as a hard fallback
for `run_genemark=false`/host-module-unavailable cases — belt and suspenders,
not a behavior change for that path).

**Backfilling the `.mod`** (Fable review points 5-6 — this is the part the
first draft under-specified): keep backfill fully centralized in
`backfill_abinitio_params.py`; do **not** have `GENEMARK_RUN` write to the
shared store directly — the store's atomic staging-dir swap (already built
for augustus/snap/parameters.json) is what makes `staleSharedParams()`'s
single-file (`parameters.json`) mtime check a reliable freshness signal for
the *whole* bundle. A genemark-only side-channel write would let a sibling
observe a refreshed `.genemark.mod` without `parameters.json`'s mtime moving,
silently serving a stale/mismatched combination with no repredict trigger.
Concretely:

- Join `GENEMARK_RUN.out.mod` (keyed by `out`) onto `FUNANNOTATE_PREDICT.out.metadata`
  by `out` **before** building `freshBackfillInput`'s `(sp, out)` pairs —
  `out` is 1:1 here too, same safety argument as the rep/indep join above.
- Manifest format changes from today's 2-field `"${sp}\t${out}"` to 3-field
  `"${sp}\t${out}\t${genemark_mod_path}"` (empty third field when this
  representative's GENEMARK_RUN took the fast-reuse path and produced no
  fresh `.mod` — nothing new to backfill in that case, species store already
  has the current one).
- `backfill_abinitio_params.py`: `backfill_species_store()` needs an explicit
  `genemark_mod: Optional[Path] = None` parameter instead of deriving
  `genemark_src` from `rep_dir/predict_misc/ab_initio_parameters/` (that
  directory won't contain a genemark.mod anymore once GeneMark moves out of
  predict — augustus/snap still live there unchanged, only genemark's source
  moves). CLI manifest-line parsing goes from 2-field to 3-field.

## New params (nextflow_schema.json + profile_funannotate.config)

- `run_genemark` (bool, default `true`) — master on/off switch; `false`
  reproduces today's fully-skipped behavior via `--auto-skip-genemark`
  (kept as the fallback, not removed).
- `genemark_mode` (`ES`|`ET`, default `ES`) — **ES only wired/tested in this
  first pass**, per your instruction. ET's requirements are already scoped
  (see below) but not built yet.
- `genemark_force_independent_strains` — comma-separated `out` list (or a
  file, mirroring `abinitio_reuse_csv`'s convention), analogous to
  `FUNANNOTATE_PREDICTION.nf`'s existing `forceIndependentSet` parameter
  that's already threaded through as a `take:` input — same mechanism,
  scoped to GeneMark specifically rather than the whole ab-initio bundle,
  since a user might want AUGUSTUS/SNAP reuse to stay on for a strain while
  forcing GeneMark to retrain independently for it (e.g. validating whether
  reuse is appropriate for a borderline-ANI strain). Nothing in the current
  subworkflow computes this today (Fable review point 8) — needs (a) a new
  `forceIndependentGenemarkSet` `take:` input mirroring `forceIndependentSet`'s
  existing shape (loaded the same offline way, e.g. a plain `Set<String>` of
  `out` values), and (b) a per-row `force_independent` boolean computed
  alongside `is_rep`/`eligible` in the `classified` `.map{}` block, threaded
  through into the `rep_genemark_in`/sibling `GENEMARK_RUN` input tuples.

## ET mode — evaluated 2026-08-12 (T-022), corrected from the first pass above

**The first pass above was wrong about which BAM feeds the hints.**
`RunGeneMarkET()`'s filter is `"\tintron\t" in line and "\tb2h\t" in line` —
tracing where a real `\tb2h\t`-tagged line actually comes from (not just
grepping for "bam2hints" and assuming) shows it is **not** the raw RNA-seq
read alignment (`args.rna_bam` / `funannotate_train.coordSorted.bam`,
processed via the external `bam2hints` binary with its *default* `--source=E`
— tagged `E`, which `RunGeneMarkET`'s filter would silently never match).
It's **transcript-alignment evidence**: `bam2ExonsHints()`
(`funannotate/library.py:1952`, funannotate's own Python BAM→hints
converter, not the AUGUSTUS binary) explicitly sets `btag = "b2h"` when
called on minimap2-aligned transcript evidence
(`predict.py:1476`). AUGUSTUS's own `blat2hints.pl` (the BLAT-alignment path)
uses the same `b2h` convention. Both paths are about spliced *transcript*
alignments, not raw short-read RNA-seq alignments.

**Good news this correction produces**: the exact input needed —
`training/transcript.alignments.bam` — is **already produced and retained**
by `FUNANNOTATE_TRAIN` (confirmed on disk for
`Penicillium_citrinum_NRRL_1841`, alongside `funannotate_train.coordSorted.bam`,
the raw-reads BAM this design's first pass mistakenly pointed at). No new
alignment work is needed; ET mode is actually *simpler* than first scoped,
not harder.

**Verified end-to-end** (`analysis/genemark_run_validation/et_eval/`, real
run, not inferred from source): the external AUGUSTUS `bam2hints` binary,
run directly against `transcript.alignments.bam` with an explicit
`--source=b2h` override (its default is `E`), produces a hints file whose
every line is independently confirmed to match `RunGeneMarkET()`'s filter
exactly (`b2h` in column 2, `intron` in column 3):

```bash
bam2hints --intronsonly --source=b2h \
    --in=<training_dir>/transcript.alignments.bam \
    --out=raw_hints.gff
awk -F'\t' 'BEGIN{OFS="\t"} $3=="intron" && $2=="b2h" {$6="500"; print}' \
    raw_hints.gff > genemark.intron-hints.gff   # score->500, exactly matching
                                                 # RunGeneMarkET()'s own rewrite
gmes_petap.pl --ET genemark.intron-hints.gff --sequence genome.fa \
    --max_intron 3000 --soft_mask 2000 --cores N --fungus
```

21,400 intron hints produced for `Penicillium_citrinum_NRRL_1841` (correctly
b2h/intron-tagged, score rewritten to 500 matching `RunGeneMarkET()` exactly)
— **but the real `gmes_petap.pl --ET` run on them FAILS**:
```
warning, no data in specified range .../histogram.pl
error, hash is empty: .../bp_seq_select.pl
error on call: .../bp_seq_select.pl --seq_in bp_region.seq --seq_out gibbs.seq
  --max_seq 4000 --bp_region_length 50 --min_bp_region_length 40
```
**Root cause, confirmed** (not guessed): this genome's intron length
distribution is `min=32 median=65 mean=81.6 max=2934` bp (measured directly
from the 21,400 hints). GeneMark's branch-point-region extraction searches a
40-50bp window *within* each intron to find spliceosomal branch-point signal
— for a median-65bp intron, that window consumes nearly the entire intron,
leaving `bp_seq_select.pl` with too little usable sequence to build a
branch-point training set. Reading `gmes_petap.pl`'s source
(`bp_seq_select.pl` calls at lines 483/539/878) shows this branch-point step
is shared between `--ET` and `--fungus`'s internal `ES_C` sub-mode — and the
earlier **successful** `--ES --fungus` run (11,116 genes) went through an
analogous step fine, because ES's candidate introns come from its own
broader genome-wide self-discovered search, not from a caller-supplied hints
file restricted to externally-aligned regions. The failure is specific to
feeding *external* hints into the branch-point step for a fungal genome with
characteristically short introns, not something wrong with the hints file's
format or content.

**FIXED (2026-08-12), root cause was not intron length — it was missing
strand assignment.** Cross-checked against BRAKER (`Gaius-Augustus/BRAKER`
on GitHub, fetched `braker.pl` + `filterIntronsFindStrand.pl` directly — the
far more battle-tested GeneMark-ET pipeline), which never feeds raw
`bam2hints` output straight to GeneMark. Every intron hint first goes
through `filterIntronsFindStrand.pl`: it checks the genome sequence at each
intron's boundary against canonical splice-site dinucleotides (GT-AG/GC-AG/
AT-AC), assigns the correct strand, and **silently drops any intron without
a canonical splice site**. My hints had `.` (unstranded) in the strand
column — branch-point signal is inherently strand-specific, so `bp_seq_select.pl`
had nothing orientable to work with. The `--bp_region_length`/intron-length
hypothesis above was a plausible-sounding but wrong lead — BRAKER uses the
same GeneMark defaults and never touches those flags either.

Corrected recipe, verified with a real `gmes_petap.pl --ET` run
(`analysis/genemark_run_validation/et_eval2/`, matching BRAKER's
`get_genemark_hints()` exactly — `braker.pl` lines 4888-4995):

```bash
bam2hints --intronsonly --in=<training_dir>/transcript.alignments.bam --out=raw_hints.gff
perl filterIntronsFindStrand.pl genome.fa raw_hints.gff --score > stranded.gff
sort -n -k4,4 stranded.gff | sort -s -n -k5,5 | sort -s -k3,3 | sort -s -k1,1 \
    | join_mult_hints.pl > genemark.intron-hints.gff   # AUGUSTUS script, already on PATH
gmes_petap.pl --ET genemark.intron-hints.gff --sequence genome.fa \
    --max_intron 3000 --soft_mask 2000 --cores N --fungus
```

20,603 of 21,400 hints (96%) passed stranding/canonicalization for
`Penicillium_citrinum_NRRL_1841`. The `--ET` run **completed successfully**:
**10,776 gene models** — same order of magnitude as the ES run's 11,116,
plausible given RNA-seq-informed splice boundaries producing a somewhat
different (not necessarily worse) model. `filterIntronsFindStrand.pl` is not
bundled with this project's funannotate/Augustus install — fetched directly
from BRAKER's repo for this test; needs to be vendored (`nextflow/bin/`) if
ET mode gets built.

One design point clarified along the way (not a correction, a confirmation):
BRAKER derives hints from raw RNA-seq read alignment (no prior assembly
needed, since BRAKER never runs Trinity), while funannotate's `b2h` hints
come from **transcript**-alignment evidence — a real difference in general,
but for this project specifically `FUNANNOTATE_TRAIN` already runs Trinity-GG
assembly + transcript alignment as part of its existing PASA-training
pipeline, so that cost is paid regardless of whether `GENEMARK_RUN` uses it.
Using `training/transcript.alignments.bam` (funannotate's existing
convention) rather than re-deriving BRAKER-style raw-read hints costs
nothing extra here and reuses infrastructure this pipeline already has.

**Revised design for `GENEMARK_RUN`'s ET branch**: gated on
`training/transcript.alignments.bam` existing and being non-empty (mirrors
`staleRnaseq()`'s existing data-presence check pattern, just pointed at a
different file). One `genemark_mode: ET` branch inside the same
`GENEMARK_RUN` process, not a separate process.

**WIRED (2026-08-12)**: `GENEMARK_RUN/main.nf` now has a real
`mode`/`training_bam`-gated ET branch (falls back to ES when `training_bam`
is empty — no RNA-seq for that genome). `genemark_mode` (previously declared
but unreferenced anywhere — confirmed by grep before this change) and a new
`trainingTranscriptBamFor(out)` helper (`utils.nf`, mirrors
`sharedGenemarkModFor()`'s shape) are now threaded through all three
`GENEMARK_RUN`/`GENEMARK_RUN_SIB` call sites in `FUNANNOTATE_PREDICTION.nf`.
Real end-to-end smoke-tested against the actual wired module (not the
by-hand recipe) via `nextflow/genemark_run_smoke.nf` — **10,780 gene
models**, matching the standalone recipe's 10,776 — plus a fast-reuse
regression check confirming the new 10-element input tuple didn't break the
ES/`--predict_with` path. Full results:
`analysis/genemark_run_validation/GENEMARK_RUN_VALIDATION.md`.

One real bug hit and fixed while building the smoke test itself, not the
module: the first test script attempt lived under
`nextflow/tests/manual/`, and `${workflow.projectDir}` (which `GENEMARK_RUN`
uses to locate the vendored `filterIntronsFindStrand.pl`) resolves to
*whichever script was originally launched*'s own containing directory, not
a fixed repo root — so a test script anywhere other than `nextflow/` itself
(sibling to `main.nf`) resolves that path wrong. Fixed by moving the test
script to `nextflow/genemark_run_smoke.nf`. Production is unaffected (
`main.nf` already lives at `nextflow/`, confirmed via the params-summary
`projectDir` printout in an earlier stub-run), but this was a real trap
worth documenting for any future standalone module smoke test.

**No-RNA-seq genomes, confirmed against the actual wired code**:
`FUNANNOTATE_RNASEQ.nf:249` explicitly mixes a `predict_no_rnaseq` branch
into `predict_input_ch` — genomes with no RNA-seq skip `FUNANNOTATE_TRAIN`
entirely (no `training/` dir gets populated at all) but still flow into
`FUNANNOTATE_PREDICTION.nf` like any other genome, since predict itself
works fine without PASA training data. `GENEMARK_RUN`'s **ES** branch has no
dependency on this at all (`gmes_petap.pl --ES` is pure genome self-training,
no `training/` data touched) — it triggers unconditionally for every genome
reaching `rep_todo`/`indep_todo`/`sibling_predict_todo`, RNA-seq or not, and
this is correct today. The **ET** branch (now wired) needs an explicit
no-RNA-seq fallback to ES, since a no-RNA-seq genome has no
`transcript.alignments.bam` for ET to read — the
`training/transcript.alignments.bam` existence gate mentioned above
naturally provides this (falls through to ES when absent, implemented in
`GENEMARK_RUN/main.nf`'s script body), but it's worth stating as an explicit
design requirement, not an accident of the gate's
mechanics, when ET actually gets built.

### Future consideration: protein hints (BRAKER's EP/ETP mode)

`gmes_petap.pl` has distinct `--EP`/`--ETP` modes (protein-only / transcript
+protein combined) — BRAKER's `get_genemark_hints()` builds a matching
`gm_hints_prot` file by filtering the master hints for `src=P`. This project
already computes comparable protein-alignment evidence: funannotate's own
`exonerate2hints()` (`library.py:4571`, run against `--protein_evidence`,
SwissProt in this pipeline) produces `hintsP`, but tags it `src=XNT`, not
`src=P` — same category of tag mismatch as the intron-hints case above, so
combining it isn't a free drop-in (would need relabeling `XNT`→`P` or a
small translator; the `CDSpart`/`intron` feature shapes look compatible at a
glance but weren't verified). Not attempted in this pass — noted as a real,
plausible extension once ET mode itself is wired into `GENEMARK_RUN`, not
before.

**BRAKER4** (`Gaius-Augustus/BRAKER4`, Snakemake-based) exists as a newer
alternative with different mode support. Not adopted here — the ES/ET recipe
validated in this pass (calling `gmes_petap.pl` directly, matching the exact
hints-preparation steps BRAKER's Perl pipeline performs) is sufficient for
what `GENEMARK_RUN` needs and keeps this project's own Nextflow orchestration
in control, rather than taking on a second full pipeline dependency.

## Open question — RESOLVED (see "Real end-to-end validation" below)

`FUNANNOTATE_PREDICT`'s `-w genemark:0` / `--auto-skip-genemark` interaction
with an externally-supplied `--genemark_gtf` needed a quick empirical check —
does `--genemark_gtf` bypass predict's own `GENEMARK_PATH`/`which_path`
preflight entirely? **Confirmed yes**, via Test 3 below: `GENEMARK_PATH` was
available and functional on the test host, but predict's internal
`RunGeneMarkES()` was never invoked — `--genemark_gtf` took priority.
Still worth a final real confirmation specifically against the container
(where `GENEMARK_PATH` is genuinely *unavailable*, not just unused) once
`FUNANNOTATE_PREDICT` itself moves onto it, but the code-path behavior this
question was actually about is now settled.

## Known gap: `species_reuse_clusters.py`'s own backfill is not wired to GENEMARK_RUN

`nextflow/bin/species_reuse_clusters.py` has its own **separate, undeduplicated**
`backfill_species_store()` — different signature from the one in
`backfill_abinitio_params.py` (which `pick_representative_strain.py` and
`BACKFILL_ABINITIO_PARAMS` both correctly import and share), always derives
`genemark.mod` from `predict_misc/ab_initio_parameters/`, no `genemark_mod`
override param. It's only referenced from `nextflow/legacy/funannotate.nf`
(the inactive legacy pipeline) — an offline maintenance script, run manually,
not part of the active DAG this design wires into, so left untouched here.
**Consequence**: for a representative predicted through the new
`GENEMARK_RUN`-wired path, `genemark.mod` no longer lands in
`predict_misc/ab_initio_parameters/` at all (see main wiring section above),
so if someone runs `species_reuse_clusters.py` against such a representative,
it will silently backfill augustus+snap only, missing genemark — same
"missing component" degradation the script already has for any component
that isn't present, not a crash, but a real behavior gap worth fixing before
this script sees regular use again post-migration.

## Real end-to-end validation (2026-08-12, `analysis/genemark_run_validation/`)

Ran the actual `GENEMARK_RUN` module (not a duplicate/hand-copy — same file
included directly) via a minimal standalone test workflow
(`analysis/genemark_run_validation/scripts/test_genemark_run.nf`), real
`gmes_petap.pl`, real genome (`Penicillium_citrinum_NRRL_1841`, same one used
in the F-008 A/B test):

- **Fresh `--ES` path**: 11,116 gene models in the produced GTF — matches
  almost exactly the A/B test baseline's `GeneMark: 11112` from predict's own
  internal call (`analysis/genemark_es_contribution/`). Confirms
  `GENEMARK_RUN`'s fresh-training path is functionally equivalent to what
  predict does internally today. ~15.5 min wall time at 8 cores (vs. the A/B
  test's ~10.2 min `genemark_es_train_seconds` for the same genome — same
  ballpark, difference plausibly session/node contention, not a code issue).
- **Fast `--predict_with <mod>` reuse path**, using a real `.mod` already on
  disk from the A/B test's `predict_misc/ab_initio_parameters/`: 100,966 gene
  models (near-identical to fresh-ES's 101,000-line GTF total, as expected —
  same genome, same trained model, different code path to get there).
  Completed in ~7 min at only **2** cores, vs. fresh-ES's ~15.5 min at 8
  cores — real, measured confirmation that `--predict_with` genuinely skips
  the expensive training loop rather than just seeding a faster convergence.
- **One test-harness-only bug found and fixed** (not a production bug):
  Nextflow's CLI parser treats `--flag ''` (an explicitly empty string value)
  as a bare boolean flag (`shared_mod` became literal `Boolean true`, not
  `''`), corrupting the rendered shell script. This only affects this
  standalone test's own CLI-argument interface — the real
  `FUNANNOTATE_PREDICTION.nf` wiring passes these values as plain Groovy
  string literals inside `tuple()` calls, never through CLI parsing, so it
  is not exposed to this footgun. Fixed in the test harness by omitting the
  flag entirely (relying on the script's own `params.shared_mod = ''`
  default) rather than ever passing an empty string via `--shared_mod ''`.

**Test 3** (real `funannotate predict --genemark_gtf` consumption of the
fresh-ES GTF, `scripts/test3_predict_consumption.sbatch`, job 27423333,
COMPLETED in 23 min): **the open question above is now empirically
answered.** Predict's log shows `GeneMark path: .../genemarkESET/4.72_lic`,
`GeneMark appears to be functional? True` (GENEMARK_PATH was available on
this host) — but the internal `RunGeneMarkES()` call was never invoked;
`--genemark_gtf` took priority and short-circuited it, exactly as
`predict.py`'s source read predicted. Confirmed directly (not just
absence-of-log-line): the "Summary of gene models" breakdown shows
`GeneMark 1 11116` — the exact count from `GENEMARK_RUN`'s fresh-ES GTF,
proving the externally-supplied file was the one actually consumed.
**Final gene count: 11,202** vs. the A/B baseline's 11,198 (predict's
internal GeneMark call) — a 4-gene difference, well inside normal
EVM/tbl2asn tie-breaking noise. **An externally-supplied `--genemark_gtf`
produces an equivalent, non-degraded annotation.** All three tests (fresh
ES, fast reuse, predict consumption) pass.

## Remaining next steps

1. ~~ET mode~~ — wired and real end-to-end validated 2026-08-12 (T-022 closed).
2. `species_reuse_clusters.py`'s undeduplicated backfill (Known gap above).
3. A `GENEMARK_RUN_SIB` invocation has not yet been exercised through the
   real `FUNANNOTATE_PREDICTION.nf` subworkflow itself (only the underlying
   module, standalone) — worth a real `--predict_scope all` run once two
   real strains of the same species with a backfilled shared store are
   available to test against.
4. Protein hints / BRAKER's EP-ETP mode (scoped above, not built).
5. Test `ghcr.io/nextgenusfs/funannotate:v1.9.0-beta.4` (released after this
   design's container-side testing) — not yet re-checked whether it resolves
   any of the packaging bugs found against beta.2/beta.3
   (`.living/learnings.md`, 2026-08-12 container smoke-test entries).
