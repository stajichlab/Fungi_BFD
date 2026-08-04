# CLEAN_MASK Refactor Plan — extract genome cleaning + repeat masking into its own pipeline

Status: **proposal, not yet implemented**
Author: drafted 2026-08-03; reviewed by Fable 2026-08-03 (findings folded in below)
Related: `REFACTOR_NEXTFLOW_PLAN.md` (the DSL2 modularisation this builds on)

---

## 1. Why

Today "prepare a genome for annotation" is split across two unrelated entry points and one
buried subworkflow:

| Where | What it does | Entry |
|---|---|---|
| `subworkflows/local/FUNANNOTATE_GENOME_PREP.nf` | SETUP_TAXONDB → GENOME_CLEAN(_BATCH) → MASKREPEAT_TANTAN_RUN | `--pipeline funannotate` (+ `--only_clean`, `--run_repeatmasker`) |
| `workflows/earlgrey_mask.nf` | SELECT_REPS → EARLGREY_BUILD_LIB → REPEATMASK_STRAIN → DELIVER_MASK | `--pipeline earlgrey_mask` |

Consequences:

- **Masking strategy is hard-coded per entry point.** `funannotate` can only tantan-mask;
  `earlgrey_mask` can only EarlGrey. Choosing a masker means choosing a *pipeline*.
- **`--only_clean true --run_repeatmasker true` is an unobvious incantation** for "just prepare
  genomes" — it works by making `FUNANNOTATE_GENOME_PREP` emit `Channel.empty()` so every
  downstream stage silently receives nothing. Three separate launcher scripts
  (`run_clean_genome.sh`, `run_clean_mask.sh`, `run_earlgrey.sh`) exist to paper over this,
  and two of them point at the same `earlgrey`/`funannotate` profiles with different params.
- **Two masking implementations that never meet.** `MASKREPEAT_TANTAN_RUN` writes
  `<asmid>.masked.fasta.gz` via `storeDir`; `DELIVER_MASK` writes the same filename via
  `publishDir` from a different pipeline. Nothing prevents them clobbering each other and
  nothing records *which* masker produced a given file.
- **`earlgrey_mask` bypasses schema validation entirely** — it never calls
  `validateParameters()`, and none of its 10 params (`cutoff_mb`, `genome_suffix`,
  `asm_stats`, `masked_dir`, `earlgrey_workdir`, `earlgrey_version`,
  `repeatmasker_version`, `repeat_taxon`, `suppress_list`, `outdir`) appear in
  `nextflow_schema.json`. Per repo convention (see `.living/`), every param must be schema-declared.

### The refactor is cheap because the interface is already a filesystem contract

The seam already exists and is already relied on by four independent consumers:

```
input_clean_genomes/<ASMID>.fa.gz            # cleaned, unmasked
input_clean_genomes/<ASMID>.masked.fasta.gz  # cleaned + soft-masked
```

- `modules/common/utils.nf::resolveGenomeFile()` — BFD genome-stats channels
- `subworkflows/local/ANI_SAMPLES.nf` — ANI genome universe
- `FUNANNOTATE_GENOME_PREP`'s `--run_repeatmasker false` branch — picks up `.masked.fasta`
  if a prior run left one
- `earlgrey_mask.nf`'s `DELIVER_MASK` — writes into it

So extracting clean+mask is **not** re-plumbing channels between pipelines. It is moving
producers of an already-published directory into one place, and letting funannotate become a
pure consumer.

---

## 2. Target shape

### 2.1 New pipeline

```
nextflow run nextflow/main.nf --pipeline clean_mask -profile clean_mask \
    --masker earlgrey --mask_group species
```

`workflows/clean_mask.nf`, dispatched from `main.nf` alongside the existing eight.

### 2.2 Masker taxonomy — two shapes, not seven

The seven requested options do **not** each need a bespoke branch. They fall into two shapes,
and the second shape is precisely what `earlgrey_mask.nf` already implements:

| Shape | Maskers | Stages |
|---|---|---|
| **A — per-genome, self-contained** | `none`, `tantan`, `repeatmasker` (stock Dfam/`-species fungi`) | `APPLY` only |
| **B — build a library, then apply it** | `repeatmasker_lib` (user-supplied `--repeat_library`), `repeatmodeler`, `edta`, `earlgrey` | `BUILD_LIB` → `APPLY` |

Shape B's `APPLY` step is the *same process* for all four: `RepeatMasker -lib <fa> -xsmall`,
i.e. today's `REPEATMASK_STRAIN`. Only the library *producer* differs. So the new code is:

```
subworkflows/local/GENOME_CLEAN.nf     # SETUP_TAXONDB + GENOME_CLEAN(_BATCH)  [moved, ~unchanged]
subworkflows/local/GENOME_MASK.nf      # dispatch on params.masker             [new, ~120 lines]
```

and four new leaf modules, of which two are thin:

| Module | Effort | Notes |
|---|---|---|
| `MASK_REPEATMASKER` | thin | stock library, `-species ${params.repeat_taxon}` |
| `BUILD_LIB_PROVIDED` | trivial | pass `--repeat_library` through as a channel value; no process |
| `BUILD_LIB_REPEATMODELER` | medium | `BuildDatabase` + `RepeatModeler -LTRStruct`; emits `*-families.fa` — same output contract as EarlGrey's `-families.fa` |
| `BUILD_LIB_EDTA` | medium | `EDTA.pl --anno 0`; emits `*.EDTA.TElib.fa`; needs a scratch-dir dance like EarlGrey's `earlgrey_workdir` |

`EARLGREY_BUILD_LIB` becomes the fourth builder with no logic change — only a rename/relocation.

**Resource labels, not a shared shape-B label.** "Same shape" describes control flow, not cost.
EarlGrey's tuned pipeline, stock RepeatMasker, RepeatModeler `-LTRStruct`, and EDTA have wildly
different runtimes (minutes to multi-day) on the same genome. Each `BUILD_LIB_*` module needs
its own `withLabel`/`withName` resource block in `conf/profile_clean_mask.config` — do not let
them inherit one generic "mask" label, or a queued EDTA job will starve/timeout alongside a
30-second `BUILD_LIB_PROVIDED` pass-through.

### 2.3 Grouping — the one genuinely new abstraction

EarlGrey masking is species-grouped (build once on a representative, apply to conspecifics);
tantan masking is per-genome. Generalise to `--mask_group`:

| `--mask_group` | Behaviour | Reuses |
|---|---|---|
| `none` (default) | every genome gets its own library | — |
| `species` | `SELECT_REPS` picks a representative per species (`--cutoff_mb` gate); members masked with the rep's library | existing `SELECT_REPS` + `bin/select_repeat_representatives.py` verbatim |
| `ani_cluster` | *future* — reuse `abinitio_reuse_csv` clusters rather than raw SPECIES | out of scope for this refactor; leave the enum value undeclared until implemented |

Shape A maskers ignore `--mask_group` (warn if set). This is the only place where the two
shapes need to know about each other.

### 2.4 Provenance — closing a real gap

Every masked genome gets a sidecar written by `DELIVER_MASK`:

```
input_clean_genomes/<ASMID>.masked.json
  {"asmid":..., "masker":"earlgrey", "masker_version":"7.2.6",
   "library_source":"<rep_asmid>", "mask_group":"species", "date":...}
```

This is what makes multiple maskers coexisting safe: today, `<ASMID>.masked.fasta.gz` is
anonymous, so a tantan run and an EarlGrey run silently overwrite each other with no record.
`DELIVER_MASK` should refuse to overwrite a file produced by a *different* masker unless
`--force_remask true`. **This is the single highest-value item in the plan** and is worth doing
even if nothing else here ships.

### 2.5 What funannotate becomes

`FUNANNOTATE_GENOME_PREP` is deleted. `workflows/funannotate.nf` gains a
`RESOLVE_PREPARED_GENOME` step (plain Groovy, no process) that maps each job to
`input_clean_genomes/<asmid>.masked.fasta.gz` → `<asmid>.fa.gz` → error.

Back-compat, so no user's muscle memory breaks in the same commit as the code move:

- `--only_clean true` → log a deprecation warning and **delegate**: `funannotate.nf` invokes
  `CLEAN_MASK()` and returns. Removed one release later.
- `--run_repeatmasker` → deprecated alias for `--masker tantan` / `--masker none`.
- `--pipeline earlgrey_mask` → alias for `--pipeline clean_mask --masker earlgrey
  --mask_group species`, warning on use.
- `run_clean_genome.sh`, `run_earlgrey.sh` → keep as thin wrappers that call
  `run_clean_mask.sh --masker {none,earlgrey}`; delete `run_clean_mask.sh`'s current body.

---

## 3. File-by-file work

### Moves (mechanical — `git mv` + include-path updates)

```
modules/funannotate/genome/GENOME_CLEAN/         → modules/genome/clean/GENOME_CLEAN/
modules/funannotate/genome/GENOME_CLEAN_BATCH/   → modules/genome/clean/GENOME_CLEAN_BATCH/
modules/funannotate/setup/SETUP_TAXONDB/         → modules/genome/clean/SETUP_TAXONDB/
modules/funannotate/genome/MASKREPEAT_TANTAN_RUN/→ modules/genome/mask/MASK_TANTAN/
modules/earlgrey/EARLGREY_BUILD_LIB/             → modules/genome/mask/BUILD_LIB_EARLGREY/
modules/earlgrey/REPEATMASK_STRAIN/              → modules/genome/mask/APPLY_REPEATMASKER/
modules/earlgrey/SELECT_REPS/                    → modules/genome/mask/SELECT_REPS/
modules/earlgrey/DELIVER_MASK/                   → modules/genome/mask/DELIVER_MASK/
```

`SETUP_TAXONDB` moving out of `modules/funannotate/setup/` is worth a second look — it is
taxonkit DB setup, used only by cleaning today, but the name suggests broader use. Verify with
`grep -rn SETUP_TAXONDB` before moving; if anything else includes it, leave it in
`modules/common/`.

### New files

```
workflows/clean_mask.nf                          (~150 lines)
subworkflows/local/GENOME_CLEAN.nf               (~70, extracted verbatim)
subworkflows/local/GENOME_MASK.nf                (~120, new dispatch)
modules/genome/mask/MASK_REPEATMASKER/main.nf    (~40)
modules/genome/mask/BUILD_LIB_REPEATMODELER/main.nf (~50)
modules/genome/mask/BUILD_LIB_EDTA/main.nf       (~50)
conf/profile_clean_mask.config                   (merge of profile_earlgrey.config + funannotate's clean labels)
conf/params_clean_mask.yaml                      (rewrite)
run_clean_mask.sh                                (rewrite)
```

### Edited

```
main.nf                        + clean_mask dispatch, + earlgrey_mask deprecation alias
workflows/funannotate.nf       - GENOME_PREP, + genome resolution, + deprecation shims
workflows/earlgrey_mask.nf     DELETE (replaced by alias)
nextflow_schema.json           + clean_mask_options group (~14 params, 10 of them
                               pre-existing but currently undeclared)
assets/schema_input.json       unchanged
conf/profile_earlgrey.config   DELETE (folded into profile_clean_mask.config)
run_lint.sh                    + clean_mask -preview case
```

### Deleted / retired

`subworkflows/local/FUNANNOTATE_GENOME_PREP.nf`, `workflows/earlgrey_mask.nf`,
`conf/profile_earlgrey.config`.

---

## 4. Cost estimate

Assumes one person, includes stub-run and lint verification but **not** waiting on
multi-day EarlGrey/EDTA real runs.

| Phase | Work | Est. |
|---|---|---|
| −1 | **Pre-flight spike (new, blocking):** `-resume` cache-hit check against the real `launchDir` after a trial module relocation, *before* the estimate below is trusted — see Fable review note under §5. Also run the runtime benchmark in §4.1. | 0.5 d |
| 0 | Schema: declare the 10 orphaned earlgrey params + 4 new ones; `validateParameters()` in the new workflow | 0.25 d |
| 1 | Module relocation + include rewrites; `run_lint.sh` + `-stub-run` green, zero behaviour change | 0.25 d |
| 2 | Extract `GENOME_CLEAN.nf`; write `GENOME_MASK.nf` dispatch for the 3 *existing* maskers (`none`/`tantan`/`earlgrey`); `workflows/clean_mask.nf`; profile merge | 1.0 d |
| 3 | Deprecation shims in `funannotate.nf` + `main.nf`; launcher-script rewrite; stub-run parity vs. current `--only_clean` and `--pipeline earlgrey_mask` DAGs (DAG-diffing itself, not just shim-writing — budgeted up from 0.5d per Fable review) | 0.75 d |
| 4 | `DELIVER_MASK` provenance sidecar + overwrite guard | 0.25 d |
| 5 | `MASK_REPEATMASKER` (stock) + `BUILD_LIB_PROVIDED` | 0.25 d |
| 6 | `BUILD_LIB_REPEATMODELER` + `BUILD_LIB_EDTA` | **unestimated — spike first, see below** |
| 7 | Docs (`doc/`, `ENVIRONMENTS_INSTALLATIONS.md`), `.living/decisions.md`, manifests | 0.25 d |
| | **Total (phases −1 to 5, 7)** | **≈ 3.5 d** |

Phases 0–4 (**≈ 2.75 d** including the −1 spike) deliver the whole structural win — separate
pipeline, schema coverage, provenance, no new tooling. Phases 5–6 are additive and can ship
later without re-touching anything from 0–4. **Recommend shipping 0–4 as one PR and 5–6 as a
second.**

**Phase 6 is deliberately left unestimated, not "1.0 d."** EDTA is notoriously fussy about
contig naming (rejects sequence IDs > 13 chars unless `--force`) and RepeatModeler's
`-LTRStruct` is a multi-day job on a 200 Mb genome. A prior version of this plan priced phase 6
at 1.0 d assuming "one real single-genome test each" would be clean; per the Fable review, that
assumption is optimistic — the EDTA contig-naming issue alone is likely to cost debugging
cycles, not just runtime. Re-scope phase 6 as its own spike once §4.1's benchmark numbers exist,
and price it then.

### 4.1 Runtime benchmark (new — do before trusting phase 5–6 estimates or default-masker debate)

Run tantan, stock RepeatMasker, and RepeatModeler (default settings, no `-LTRStruct` yet) on
one representative fungal genome (~30–40 Mb, mid-range for BFD) and record wall-clock time and
peak memory for each. This directly informs:

- Open question §6.5 (retire tantan as default) — right now the "~100× cheaper" figure is an
  assertion, not a measurement.
- Open question §6.4 (is per-genome RepeatModeler affordable at ~8k genomes) — extrapolate
  from the single-genome number before deciding `--mask_group species`-only is required.
- Phase 6 re-estimation once EDTA/RepeatModeler costs are known rather than assumed.

Do not block phases 0–4 on this — it can run in parallel on a spare node — but do not finalize
phase 5–6 scope or touch the tantan-default question in §6.5 until the numbers are in hand.

---

## 5. Risks and how each is contained

| Risk | Severity | Containment |
|---|---|---|
| Losing `-resume`/`storeDir` cache on ~8k already-cleaned genomes | **High** | Nothing here changes `storeDir "${launchDir}/input_clean_genomes"` or the output filenames. Cleaning is `storeDir`-cached by *path*, not by task hash, so a relocated module with an identical output name is still a cache hit. **Updated per Fable review: verify *before* phase 1, as a phase −1 pre-flight spike (§4)** — not "during phase 1" as originally written. If this claim is wrong, the plan's core premise ("cheap because it's a filesystem contract") collapses at 8k-genome scale, so it must be settled before the rest of the estimate is trusted, not discovered mid-implementation. |
| In-flight `--pipeline earlgrey_mask` jobs at merge time recompute expensive BUILD_LIB work | **High (new, per Fable review)** | Module relocation changes Nextflow's task hash (it includes process source path), even though `storeDir` output filenames stay identical. So a job mid-flight on `--pipeline earlgrey_mask` when this merges could have its *expensive* Shape-B `BUILD_LIB_*` step recomputed on `-resume`, not just the cheap cleaning step the original risk table worried about. **Containment: confirm there are no in-flight `earlgrey_mask` runs before merging phases 0–4, or hold the merge for a natural gap in the run queue.** This was previously unaddressed. |
| Different maskers clobbering `<ASMID>.masked.fasta.gz` | **High** | §2.4 provenance sidecar + overwrite guard. This risk exists *today*, unguarded; the refactor is what makes it fixable. |
| `funannotate.nf` invoking `CLEAN_MASK()` cross-workflow (the `--only_clean` delegation path) is treated as "an alias" but is a bigger structural change | Medium (new, per Fable review) | DSL2 doesn't make invoking one workflow from another trivial (channel/param scoping across workflow boundaries needs to be worked out explicitly). Before phase 3, write a concrete code sketch of the `CLEAN_MASK()` invocation from `funannotate.nf`, not just the bullet point in §2.5. |
| No enforcement of masker consistency within a downstream comparison set | Low–Medium (new, per Fable review) | Provenance sidecars (§2.4) *record* which masker produced a genome, but nothing warns or blocks when a BUSCO/comparison stage mixes tantan- and EarlGrey-masked genomes in the same analysis. Out of scope for phases 0–4; flag as a phase-7-or-later follow-up if it starts to bite. |
| funannotate now fails on an unprepared genome instead of preparing it | Medium | Intentional (that is the decoupling), but must fail *loudly with the exact `run_clean_mask.sh` command to run*, not silently drop the sample the way the current `log.warn "No cleaned genome"` filter does. |
| `beforeScript` module loads in `profile_earlgrey.config` | Medium | Already a known repo hazard (CLAUDE.md). The existing earlgrey processes load modules **both** in `beforeScript` and inline in `script:` — the inline load is what actually works. When merging profiles, **drop the `beforeScript` blocks** rather than carrying them over, so nobody later "cleans up" the redundant-looking inline load. |
| Groovy `$` interpolation in the new mask modules' `script:` blocks | Medium | CLAUDE.md rule. `BUILD_LIB_EDTA`/`BUILD_LIB_REPEATMODELER` both need `$(...)`-heavy shell; write the logic in Groovy where possible, escape `\$` otherwise, and stub-run before submitting. |
| `--only_clean` delegation changes DAG shape for existing users mid-run | Low | Deprecation path keeps the same output paths; worst case a user's in-flight `-resume` recomputes nothing because `storeDir` still hits. |
| `mask_group: ani_cluster` scope creep | Low | Explicitly out of scope; do not declare the enum value. |

---

## 6. Open questions for review

1. **Is `SETUP_TAXONDB` funannotate-specific or general?** Determines whether it moves to
   `modules/genome/clean/` or `modules/common/`.
2. **Should `clean` and `mask` be one pipeline or two?** This plan says one (`clean_mask`) with
   `--masker none` covering clean-only, mirroring today's `--only_clean`. The alternative —
   `--pipeline clean` and `--pipeline mask` — is more orthogonal but duplicates the samples.csv
   filtering preamble and makes the common "do both" case two launches.
3. **Does `--masker` belong on `funannotate` at all after this?** Plan says no (funannotate
   consumes only). But a user annotating a handful of new genomes may reasonably want one
   command. Counter-proposal: `funannotate` gains `--auto_prepare true` which invokes
   `CLEAN_MASK()` for missing genomes only.
4. **Is per-genome RepeatModeler affordable at BFD scale (~8k genomes)?** If not,
   `repeatmodeler` should be documented as `--mask_group species`-only.
5. **Retire tantan as the default?** EarlGrey/RepeatModeler soft-masking is materially better
   for fungal gene prediction, but tantan is *asserted* to be ~100× cheaper and is what every
   already-annotated genome in the DB used. Changing the default silently makes old and new
   annotations non-comparable. Recommend: keep `tantan` as default, document the tradeoff —
   **but replace the "~100×" assertion with the §4.1 benchmark numbers before finalizing this
   decision** (per Fable review: the current figure is unmeasured).
6. **Should mixed-masker comparison sets be flagged?** (new, per Fable review) Provenance
   sidecars (§2.4) record which masker produced each genome, but nothing today warns or blocks
   a downstream comparison (e.g. BUSCO) from silently mixing tantan- and EarlGrey-masked
   genomes. Not required for phases 0–4; worth deciding whether it becomes a phase-7 doc note
   or an actual enforcement check later.

---

## 7. Verification gates (per phase)

1. `bash nextflow/run_lint.sh` clean (no new errors; warnings advisory).
2. `-stub-run` task-count parity: `--pipeline clean_mask --masker tantan` must produce the same
   process set as today's `--pipeline funannotate --only_clean true`, and
   `--masker earlgrey --mask_group species` the same as today's `--pipeline earlgrey_mask`.
   (This is the same clean-stub-comparison technique `REFACTOR_NEXTFLOW_PLAN.md` §used for
   `BFD.nf` vs `main.nf` — 44 tasks each.)
3. Dry `-resume` against the real launch dir showing zero re-cleaning.
4. One real end-to-end run on a single small ASMID per masker before declaring that masker
   supported.
