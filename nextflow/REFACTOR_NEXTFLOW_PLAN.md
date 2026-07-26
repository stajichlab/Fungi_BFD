# Refactoring Plan: Nextflow Module & Workflow Structure

**Branch:** `refactor_modules`
**Goal:** Migrate all pipelines to an nf-core-inspired `modules/` + `subworkflows/` + `workflows/`
layout, with conventions that survive Nextflow's strict language spec.
**Runtime:** Nextflow 26.04 installed; `nextflow.config` pins `nextflowVersion = '>=25.10'`.

> **Read §2 before writing any module.** It is the normative part of this document.
> Everything else is context, sequencing, and verification procedure.

---

## 1. Status

| Phase | Scope | State |
|-------|-------|-------|
| 1a | `modules/common/utils.nf` (shared functions) | ✅ Done |
| 1b | `modules/common/` process modules | ⚠️ Partial — `SETUP_SYMLINKS` done; `SAMPLES_READ` is dead+incorrect, must be deleted |
| 1c | `modules/BFD/` — 31 process modules | ✅ Done, verified byte-identical |
| 1d | `workflows/BFD.nf` + `main.nf` | ⚠️ Done but not decomposed into subworkflows |
| 1e | Subworkflow decomposition (`BFD_FUNCTIONAL`, `BFD_GENOME_STATS`, `BFD_MERGE`) | ❌ Not started |
| 1f | Cut `run_*.sh` over to `main.nf`, delete legacy `BFD.nf` | ❌ Not started |
| 2 | funannotate | ❌ Blocked on §2 conventions being frozen |
| 3 | ANI / comparative / earlgrey / phyling | ❌ Not started |

**Phase 1 is NOT complete.** Process extraction is done and verified; workflow-level
decomposition and the production cutover are not. Do not start Phase 2 until 1e and 1f land.

### Verified equivalence (2026-07-25)

The extraction was checked mechanically, not by eye:

- All 32 process blocks diffed against the original `BFD.nf` — byte-identical apart from
  the intentional `tablesDir()` change and (regrettably) stripped comments.
- Clean-output stub comparison of `BFD.nf` vs `main.nf`: **44 tasks each, identical process
  set and counts, 81 published output files with identical paths and identical decompressed
  content.** Command in §7.
- `nextflow lint`: 35 errors → 3. The 3 remaining (`SampleUtils`) are pre-existing and also
  present in the original `BFD.nf`; the refactor no longer adds strict-syntax debt.

---

## 2. Conventions — normative

These are binding for every module and workflow added from here on, including Phases 2 and 3.
Where an existing file violates one, it is listed as remediation work in §5.

### 2.1 No Groovy classes under `lib/`

**Rule:** shared logic goes in an includable `.nf` file exposing plain functions. Never a
class on the `lib/` classpath.

**Why:** `lib/*.groovy` auto-classpath loading is not part of Nextflow's strict language spec.
`nextflow lint` reports ``` `Utils` is not defined ``` for every reference, which breaks the
language server and will break outright as the strict parser becomes the default. Classes also
cannot see the implicit `params`, `log` and `workflow` variables, so every helper has to thread
them through as extra positional arguments — the first cut of this refactor turned
`clearIfStale(prot, outs)` into `Utils.clearIfStale(prot, outs, workflow.stubRun, log)` at 17
call sites, with two transposable arguments and no compile-time checking.

```nextflow
// modules/common/utils.nf
def tablesDir() { ... }                              // reads params directly — fine
def clearIfStale(inputFile, List storedOutputs) {    // reads workflow.stubRun and log — fine
    if (workflow.stubRun) return
    ...
}
```

```nextflow
// consumer
include { clearIfStale } from '../modules/common/utils.nf'
```

**Status:** done for `Utils`. `lib/SampleUtils.groovy` is still a class because
`funannotate.nf`, `compare_ANI.nf`, `comparative_genomics.nf`, `query_ANI.nf`, `phyling.nf` and
`development/*.nf` all depend on it. Migrate it **once**, in Phase 2, moving every consumer in
the same commit. Do not create a second copy of `makeSampleTag`/`cleanStrain` in the meantime —
a duplicated sample-tag implementation is exactly the bug described in §2.7.

### 2.2 Directory name == process name

**Rule:** `modules/<group>/<NAME>/main.nf` defines exactly `process <NAME>`. One process per
module. No prefixes that appear in the process but not the path.

Currently violated by 9 modules: `modules/BFD/PFAM/main.nf` defines `RUN_PFAM`, and likewise for
CAZY, MEROPS, SIGNALP, TMHMM, TARGETP, IDP, WOLFPSORT, PREDGPI. The other 22 modules already
comply. Rename the directories to `RUN_PFAM/` etc. (cheaper than renaming the processes, which
would invalidate `withName:` selectors and `-resume` caches). Fix before Phase 2 quadruples the
module count.

### 2.3 Group modules by concern, not in one flat namespace

`modules/BFD/` currently holds 31 flat siblings mixing three concerns. Nest them so the
subworkflows in §2.5 have an obvious thing to import:

```
modules/BFD/functional/   RUN_PFAM, RUN_CAZY, … (9 external-tool wrappers)
modules/BFD/genestats/    BATCH_AA_FREQ, CALC_GENE_STATS, BUSCO_GENOME, … (7)
modules/BFD/merge/        MERGE_PFAM, MERGE_SAMPLES, … (15)
```

### 2.4 Sample identity travels in a `meta` map

**Rule for new and migrated code:** processes take `tuple val(meta), path(...)`, where

```nextflow
meta = [
    id      : basename,     // filesystem-safe SPECIES_STRAIN tag; the primary key
    locustag: locustag,
    species : species,
    strain  : strain,
    asmid   : asmid,        // where relevant
]
```

**Why this is not optional.** Every module today takes
`tuple val(locustag), val(basename), val(species), val(strain), path(proteins)`. Positional
tuples mean: adding a field edits every module signature and every call site; nothing checks
that argument 3 is species rather than strain across 31 modules; and no nf-core module or
subworkflow can ever be dropped in, because they all expect `tuple val(meta), path(x)`.

The refactor as implemented froze the positional convention into 31 files, so the migration
now costs 31 call sites instead of the 9 inline blocks it would have cost before. **It gets
worse, not better, with every module added — so do this before Phase 2, not after.**

If you decide against `meta`, delete this subsection and write the reason in §8 so it reads as
a decision rather than an oversight. Do not leave it ambiguous.

### 2.5 `workflows/*.nf` orchestrates; it does not contain the pipeline

**Rule:** a file in `workflows/` wires named subworkflows together and does nothing else. Any
`workflows/*.nf` over ~150 lines means the subworkflows are missing.

`workflows/BFD.nf` is 522 lines and inlines every RUN, stats and MERGE call in one
`workflow BFD {}` block — structurally the same monolith as the original, moved one directory
over, and now arguably harder to read because seeing what a step does means opening one of 31
files. Extract:

```
subworkflows/local/BFD_FUNCTIONAL.nf     take: proteins_ch          → 9 RUN_* + 9 MERGE_*
subworkflows/local/BFD_GENOME_STATS.nf   take: the 7 stats channels → BATCH_*/CALC_*/BUSCO_*
subworkflows/local/BFD_MERGE.nf          take: stats outputs        → 6 stats MERGE_* + MERGE_SAMPLES
subworkflows/local/INPUT_SETUP.nf        ✅ already extracted
```

The `merge_all` / `--taxon` branching (`workflows/BFD.nf:265-366` and `:437-521`) is the most
duplicated logic in the file and belongs inside `BFD_MERGE`, expressed once.

### 2.6 Modules should not read `params.*` directly (target state)

Output locations belong in config, not in the module, so a module is portable and its layout is
overridable per invocation:

```nextflow
// conf/modules.config          ← to be created
process {
    withName: 'MERGE_.*'   { publishDir = [path: { tablesDir() }, mode: 'copy'] }
    withLabel: 'functional' { storeDir  = { "${params.outdir}/${task.process.tokenize(':')[-1].toLowerCase()}" } }
}
```

Every module currently hardcodes `storeDir "${params.outdir}/pfam_hmmscan"` or
`publishDir path: { tablesDir() }`. This was a deliberate Phase-1 scoping call and is
**carried debt, not resolved debt**. Schedule it as Phase 1.5 alongside §2.2/§2.3/§2.4 — those
are all "touch every module once" changes and should be done in a single pass.

### 2.7 One implementation of sample-tag derivation

`SampleUtils.makeSampleTag`/`cleanStrain` is the single source of truth (its own header notes
the Python equivalents in `collect_busco_stats.py` and `fix_low_trinity.py` that must stay in
sync). Any code deriving a basename from SPECIES+STRAIN calls it.

**`modules/common/SAMPLES_READ/main.nf` violates this and must be deleted.** It is never
included anywhere, and its inline-Python derivation (lines 19-23) silently diverges: no quote
stripping, no asterisk rules, no colon handling. It is advertised in this plan as the shared
canonical parser for Phase 2, which makes it a trap rather than dead code — the next person to
wire it in gets sample tags that disagree with the rest of the pipeline and with the DuckDB
build. Delete it, or reimplement it as a thin call into the canonical helper.

### 2.8 Scope process labels per pipeline

`withLabel: 'merge'` and `'setup'` are each defined in five config files
(`profile_BFD.config`, `profile_BFD_k8.config`, `profile_comparative.config`,
`profile_phyling.config`, `test.config`). Harmless while one profile is selected at a time —
but profile composition is already in use (`run_lint.sh` uses `-profile funannotate,test`), and
a single `main.nf` entry point makes composition the norm.

**Rule:** prefix generic labels with the pipeline: `bfd_merge`, `funannotate_merge`,
`bfd_setup`. Do this before `modules/funannotate/` exists, not reactively "if conflicts arise".

### 2.9 Preserve "why" comments verbatim

The extraction stripped explanatory comments from every module. Some were load-bearing:

- `CALC_GENE_STATS` — why `bedtools` is re-loaded inside the script even though the `genestats`
  `beforeScript` already loads it (Nextflow runs `.command.sh` as `bash -l`, which rebuilds
  PATH and drops it). Without this note, someone deletes the "redundant" line and breaks tRNA
  codon enrichment.
- `SETUP_SYMLINKS` — why the output is a real file rather than `val(true)` (`val(true)` is
  always identical, so `-resume` skips the script when the species set changes).
- `MERGE_IDP` — why inputs are staged into separate `idp/` and `sum/` subdirs (so the globs
  cannot cross-match).
- `toManifest` — the staleness-via-content-hashing rationale, and why `path 'inputs/*'` staging
  was abandoned at ~8k genomes ("Argument list too long").

**Rule:** extraction is verbatim, comments included. A comment explaining a non-obvious
constraint is the highest-value thing in the file.

### 2.10 Named entry workflows

`main.nf` currently claims to route `-profile` to a workflow and does not — it unconditionally
calls `BFD()`. Once a second pipeline is migrated this breaks. Use explicit entries:

```nextflow
// main.nf
workflow BFD          { include… ; BFD_WF() }
workflow FUNANNOTATE  { … }
```
invoked as `nextflow run nextflow/main.nf -entry BFD -profile BFD`. Fix `main.nf:6`'s comment
either way — a comment describing behaviour the code does not have is worse than none.

### 2.11 No dual maintenance

**Rule:** the commit that lands a migrated pipeline also repoints its `run_*.sh` launchers and
deletes the legacy `.nf`. Never both at once, not even for a day.

`run_functional.sh:26`, `run_functional_all.sh:26`, `run_test.sh:24,32` and `run_lint.sh:13`
still invoke `nextflow/BFD.nf`. Production therefore still runs the monolith and the migrated
path has never executed a real workload. Two copies of 32 processes with nothing enforcing
agreement is the highest-probability source of a silent production bug in this whole refactor.

---

## 3. Target architecture

```
nextflow/
├── main.nf                       Entry point — named entry workflows only (§2.10)
├── nextflow.config
├── nextflow_schema.json
├── assets/schema_input.json
│
├── lib/
│   └── SampleUtils.groovy        LEGACY — migrate to a function module in Phase 2 (§2.1)
│
├── modules/
│   ├── common/
│   │   ├── utils.nf              Shared functions: tablesDir, clearIfStale
│   │   └── SETUP_SYMLINKS/main.nf
│   └── BFD/
│       ├── functional/RUN_*/main.nf     (9)
│       ├── genestats/{BATCH,CALC,BUSCO}_*/main.nf  (7)
│       └── merge/MERGE_*/main.nf        (15)
│
├── subworkflows/local/
│   ├── INPUT_SETUP.nf
│   ├── BFD_FUNCTIONAL.nf
│   ├── BFD_GENOME_STATS.nf
│   └── BFD_MERGE.nf
│
├── workflows/
│   └── BFD.nf                    ≤150 lines: channel construction + subworkflow calls
│
├── conf/
│   ├── modules.config            NEW — publishDir/storeDir/ext.args (§2.6)
│   └── profile_*.config          Pipeline-scoped labels (§2.8)
└── tests/
```

---

## 4. Sequencing

Each step ends with the §7 verification. Do not start the next until it passes.

**Phase 1e — finish the decomposition.** Extract the three BFD subworkflows (§2.5). BFD is the
smallest, least stateful pipeline and the right place to practise the pattern.

**Phase 1f — cut over.** Repoint `run_*.sh` at `main.nf`, delete `BFD.nf`, delete
`modules/common/SAMPLES_READ/` (§2.7), fix `main.nf`'s routing comment (§2.10).

**Phase 1.5 — freeze conventions, apply retroactively.** One pass over the 31 existing BFD
modules applying §2.2 (dir names), §2.3 (grouping), §2.4 (`meta` map), §2.6 (`modules.config`),
§2.8 (label scoping), §2.9 (restore comments). Applying the conventions to already-working code
is a cheap test of whether the conventions are right — and BFD is where a mistake is
recoverable.

**Phase 2 — funannotate.** Only after 1.5. Migrate `SampleUtils` to a function module with all
consumers in one commit (§2.1). funannotate's stateful multi-stage processes
(`FUNANNOTATE_TRAIN`/`PREDICT`/`UPDATE`, with retry and cache logic) are the wrong place to
discover convention gaps for the first time.

**Phase 3 — ANI, comparative, earlgrey, phyling.** Mechanical by then.

### Why this order

The original plan sequenced Phases 2 and 3 as straightforward repeats of Phase 1. But Phase 1 as
implemented skipped the decomposition and left the `meta`/naming/config questions unanswered, so
"repeat the pattern" for a 2,439-line pipeline means replicating the same debt at 2-3x scale.
Freeze the conventions on the small pipeline first.

---

## 5. Outstanding remediation

| # | Item | Ref | Priority |
|---|------|-----|:--------:|
| 1 | Delete `modules/common/SAMPLES_READ/` — dead and incorrect | §2.7 | **High** |
| 2 | Repoint `run_*.sh` to `main.nf`; delete legacy `BFD.nf` | §2.11 | **High** |
| 3 | Extract `BFD_FUNCTIONAL` / `BFD_GENOME_STATS` / `BFD_MERGE` | §2.5 | **High** |
| 4 | Decide and apply the `meta` map | §2.4 | **High** |
| 5 | Migrate `SampleUtils` off `lib/` (3 remaining lint errors) | §2.1 | Medium |
| 6 | Rename 9 `PFAM/`→`RUN_PFAM/` dirs; regroup into functional/genestats/merge | §2.2, §2.3 | Medium |
| 7 | Create `conf/modules.config`; remove `params.*` from modules | §2.6 | Medium |
| 8 | Prefix `merge`/`setup` labels per pipeline | §2.8 | Medium |
| 9 | Restore stripped "why" comments in all 32 modules | §2.9 | Medium |
| 10 | Fix `main.nf` routing comment / add named entries | §2.10 | Low |
| 11 | Cover the `merge_all=true` glob branch in tests | §7 | Low |
| 12 | `params.asmid` appears unused in BFD — confirm and remove | — | Low |
| 13 | **Race:** `MERGE_AA_FREQ`/`MERGE_CODON_FREQ` fire nondeterministically under `merge_all=false` (`.filter { it.exists() }` races `publishDir`). Pre-existing; reachable in production via `--taxon`. Gate on the published output rather than testing the disk | §7.3 | **High** |

---

## 6. Config compatibility

`withLabel:` / `withName:` selectors in `conf/` match module `label` / process names unchanged,
so no config edits were needed for the extraction. Two caveats:

- Selectors match the **process name**, not the module path — so §2.2's directory renames are
  safe but renaming a *process* invalidates `withName:` selectors and `-resume` caches.
- In a subworkflow, process names become qualified (`BFD:MERGE_PFAM` in the log). `withName:`
  still matches on the simple name, but `withName: 'BFD:MERGE_.*'` is available if a
  pipeline-scoped override is needed later.

---

## 7. Verification protocol

Run **all four** after every step in §4.

### 7.1 Strict-syntax lint

```bash
nextflow lint nextflow/main.nf nextflow/workflows/ nextflow/subworkflows/ nextflow/modules/
```
Errors must not increase. Current baseline: **3 errors** (`SampleUtils`, item 5 above), 36 files
clean. Warnings (`Channel`→`channel`, implicit closure params) are advisory; clear them
opportunistically.

### 7.2 Preview

```bash
nextflow run nextflow/main.nf -c nextflow/nextflow.config -profile test -preview
```
Proves the DAG assembles and all `include` paths resolve.

### 7.3 Golden-output diff — the real safety net

A `-stub-run` alone only proves the file parses. Two things previously made it prove nothing at
all, both now fixed in `conf/test.config`:

- `pep_dir` was assigned twice, the second overriding it to a directory that does not exist.
  Zero samples passed the `prot.exists()` filter and **zero processes ran** while the run still
  reported `[SUCCESS]`.
- `nextflow.config` defaults `skip_merge = true`, which silently dropped every `MERGE_*` step —
  the half of the DAG most likely to break during refactoring.

Fixtures now exist for all five input types under `tests/data/input/{pep,cds,gff3,dna,trna}/`.

Because the RUN processes are `storeDir`-cached and `tests/output/` contains committed outputs,
a default-path run short-circuits the entire RUN half (15 tasks instead of 44). **Always
redirect the output dirs to a clean scratch location:**

```bash
OUT=$(mktemp -d)
for v in old:BFD.nf new:main.nf; do
  n=${v%%:*}; f=${v##*:}
  nextflow -q run nextflow/$f -c nextflow/nextflow.config -profile test -stub-run \
    -work-dir $OUT/wk_$n \
    --outdir $OUT/$n/function --genome_stats_outdir $OUT/$n/genome_stats \
    --stats_outdir $OUT/$n/seqstats --tables $OUT/$n/tables > $OUT/$n.log 2>&1
done
# same processes, same counts?  (subworkflows qualify names as BFD:MERGE_PFAM — strip the prefix)
pat() { grep -oE "\[PROCESS [a-f0-9/]+\] [A-Za-z_0-9:]+" "$1" | awk '{print $3}' | sed 's/^BFD://' | sort | uniq -c; }
echo "tasks: old=$(pat $OUT/old.log | awk '{s+=$1} END {print s}') new=$(pat $OUT/new.log | awk '{s+=$1} END {print s}')"
diff <(pat $OUT/old.log) <(pat $OUT/new.log)

# same output files?
diff <(cd $OUT/old && find . -type f | sort) <(cd $OUT/new && find . -type f | sort)

# same contents?  (compare decompressed — gzip embeds an mtime, so bytes differ legitimately)
while read -r f; do
  a=$OUT/old/${f#./}; b=$OUT/new/${f#./}
  case "$f" in
    *.gz) diff -q <(gzip -cd "$a") <(gzip -cd "$b") >/dev/null || echo "CONTENT DIFF: $f" ;;
    *)    diff -q "$a" "$b" >/dev/null || echo "CONTENT DIFF: $f" ;;
  esac
done < <(cd $OUT/old && find . -type f | sort)
```

Baseline established 2026-07-25: **43 of 44 tasks match exactly, with identical output paths
and identical decompressed content.** Keep the legacy `.nf` until its replacement reproduces
this, then delete it (§2.11).

**The comparison is not exactly reproducible, and the cause is a pipeline bug, not the
refactor.** `MERGE_AA_FREQ` and `MERGE_CODON_FREQ` fire or don't, run to run, for the *same*
file — six clean runs of the unmodified `BFD.nf` gave AA-only (44 tasks), CODON-only (44), and
both (45). Cause: `workflows/BFD.nf:492` and `:504` use `.filter { it.exists() }` to test the
disk at channel-evaluation time, racing `BATCH_*_FREQ`'s `publishDir`. The `aa_sync` barrier
gates on the batch *task* completing, not on its output being published. Treat a 44/45 task
count and a one-file difference in `tables/All_Taxa/{aa_freq,codon_freq}.csv.gz` as expected
noise until item 13 is fixed; any *other* difference is a real regression.

**Known gap:** the profile sets `merge_all = false`, so only the current-run merge branch is
covered. The `merge_all = true` glob branch — the production default — is untested (item 11).

### 7.4 Output-schema validation

```bash
python3 nextflow/tests/validate_outputs.py
```

> Note: `tests/output/` holds committed run artifacts that mask `storeDir` behaviour (see 7.3).
> Consider moving them to `tests/expected/` and treating them as golden files, or deleting them
> — but that is a change to the existing test harness, out of scope for this refactor.

---

## 8. Decision log

| Date | Decision | Rationale |
|------|----------|-----------|
| 2026-07-25 | Shared helpers are functions in `modules/common/utils.nf`, not classes in `lib/` | Strict-syntax compatibility; classes cannot see `params`/`log`/`workflow` (§2.1) |
| 2026-07-25 | `conf/test.config` runs on `executor = 'local'` | Stub tasks are trivial; keeps the test profile off the shared SLURM queue instead of queueing ~30 jobs to echo placeholders |
| 2026-07-25 | Test profile enables merges (`skip_merge = false`, `merge_all = false`) | Merge steps are the most refactor-fragile half of the DAG; `merge_all = false` avoids depending on leftover files in `tests/output/` |
| — | **`meta` map: adopt or reject?** | **Open — blocks Phase 1.5 (§2.4)** |
| — | **`storeDir`/`publishDir` to `conf/modules.config`?** | **Open — recommended, §2.6** |

---

## 9. Risks

| Risk | Likelihood | Mitigation |
|------|:----------:|-----------|
| Legacy `.nf` and migrated path diverge silently | **High** | §2.11 — cut over and delete in the same commit. Currently live: all `run_*.sh` still point at `BFD.nf` |
| `-stub-run` passes while exercising nothing | **Occurred** | §7.3 — assert the task count, not just the exit status |
| Positional-tuple convention hardens across 3 pipelines | **High** | §2.4 — settle `meta` before Phase 2 |
| `lib/` classes break under strict syntax | Medium | §2.1 — done for `Utils`; `SampleUtils` outstanding |
| Label collisions once profiles compose | Medium | §2.8 — scope labels before `modules/funannotate/` exists |
| Load-bearing comments lost in extraction | **Occurred** | §2.9 — restore; extract verbatim from here on |
| `-resume` confusion across the layout change | Low | RUN processes are `storeDir`-cached, so existing outputs are honoured regardless of task hash — self-healing, but verify rather than assume |
| Process-name collisions across pipelines | Low | Names are already unique; module path disambiguates |
