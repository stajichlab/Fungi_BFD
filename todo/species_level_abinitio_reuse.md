# Plan: Species-level ab initio parameter reuse for funannotate PREDICT

Status: DESIGN FINAL — not yet implemented. Grounded in current code, reviewed twice by
Fable (draft and post-decision passes), backed by a real stage-timing measurement, and
fully interrogated/decided with the user (§6). Ready for implementation.

## 1. Problem, precisely

For species represented by dozens–hundreds of strains, `funannotate predict` currently
re-trains SNAP, GeneMark-ES (self-training), and AUGUSTUS from scratch for **every strain**,
even though these ab initio models should converge to nearly identical gene-structure
parameters for strains of the same species.

The SNAP/GeneMark/Augustus training happens inside `FUNANNOTATE_PREDICT`
(funannotate.nf:1300–1468), not inside `FUNANNOTATE_TRAIN`. `FUNANNOTATE_TRAIN`
(funannotate.nf:1117–1289) only does Trinity-GG + PASA alignment, and Trinity is *already*
shared at species level (see §2). PASA dominates `FUNANNOTATE_TRAIN` wall time (median
84.8%, `analysis/funannotate_train_stage_timing/`) and is **not** a candidate for
cross-strain sharing — PASA aligns each strain's own transcripts to its own genome
coordinates, so it must run per-strain regardless.

Confirmed mechanism for reuse: `funannotate predict -p <parameters.json>` accepts a JSON
pointing to pre-trained SNAP hmm / GeneMark mod / AUGUSTUS species dir and skips
re-training them:

```json
{
  "augustus": [{"path": ".../predict_misc/ab_initio_parameters/augustus/species/penicillium_freii_ms0197"}],
  "genemark": [{"path": ".../predict_misc/ab_initio_parameters/penicillium_freii_ms0197.genemark.mod"}],
  "snap":     [{"path": ".../predict_misc/ab_initio_parameters/penicillium_freii_ms0197.snap.hmm"}]
}
```
(`-p/--parameters` is the only flag that supplies all three; `--augustus_species` and
`--genemark_mod` exist standalone but there is no standalone SNAP flag.)

**GlimmerHMM is not applicable to this plan.** Verified: the pipeline's PREDICT command
already zeroes its EVM weight (`-w codingquarry:0 glimmerhmm:0`, funannotate.nf:1420), the
log confirms `"Skipping GlimmerHMM prediction as weight set to 0"`, and a real
`parameters.json`'s `glimmerhmm` entry is empty (`[{}]`). No trained GlimmerHMM model is
ever produced in this pipeline as configured, so there's nothing to share.

Caveat: the JSON's `path` fields point into the process's ephemeral Nextflow work-dir, so
the JSON as originally produced is not directly reusable — paths must be regenerated to
point at a durable, species-level store (§4.2).

### Cost premise — measured, not assumed

A single hand-inspected log initially suggested RNA-seq hints prep dominated PREDICT wall
time. An aggregate stage-timing analysis across a **real sample of 400 completed PREDICT
runs** (`analysis/funannotate_predict_stage_timing/`, `scripts/profile_predict_stage_timing.py`,
mirroring the existing `funannotate_train_stage_timing/` methodology) shows that impression
was an outlier, not the typical case:

| Component | Share of PREDICT wall time |
|---|---|
| RNA-seq hints prep | median 0.0%, mean 9.5% (most runs have none logged) |
| **Ab-initio training (GeneMark-ES + Augustus + SNAP)** | **median 39.8%, mean 38.7%** |

Critically, Augustus training cost is **bimodal depending on evidence availability**, split
out from the same 400-run sample:

| Augustus training path | Share of sample | Median training time | Median % of that run's wall time |
|---|---|---|---|
| PASA (RNA-seq available) | 116/400 (29%) | 0.3 min | 0.2% — negligible |
| **BUSCO-only (no RNA-seq)** | **267/400 (67%)** | **25.9 min (p90 61.2 min)** | **18.8%** |

The BUSCO-only path — a full genome-mode BUSCO run (`funannotate-BUSCO2.py`) plus Augustus
training against its output — is not a minor edge case; it's the **majority** case in this
dataset (67% of sampled strains) and it's expensive. GeneMark-ES (251.9 total sampled hours
across 1135 occurrences) and Augustus's BUSCO-path training (143.8 total sampled hours
across 789 occurrences) are comparable in aggregate cost, not GeneMark-ES-dominant as
originally assumed. SNAP training uses the same BUSCO-derived training gene models when no
PASA evidence exists, so it inherits this dependency too, even though SNAP's own training
step is individually cheap (14.8 total sampled hours, 1.4% of wall time).

**Conclusion: sharing ab-initio parameters is worth doing for all three predictors, not
just GeneMark-ES** — reuse avoids running the expensive BUSCO-seeding step *at all* for
non-representative, non-RNA-seq strains, which is where the real savings are (§6 decision 4).

## 2. Existing precedent to build on

`RNASEQ_PREPARE` (funannotate.nf:1005–1111, wired at :2092–2171) already implements the
shape of what's needed, just for Trinity instead of ab initio params: assemblies grouped by
`species_tag`, first-row-in-samples.csv picked as representative, Trinity-GG run once and
`storeDir`-cached, every other strain reuses it via `--trinity <shared>` while still running
its own PASA.

This plan reuses that grouping/caching shape for ab initio parameters, with one addition
RNASEQ_PREPARE does not have: an ANI gate. Trinity-GG sharing is evidence-agnostic and
low-risk even for a mislabeled/divergent strain; sharing trained AUGUSTUS/SNAP/GeneMark
*parameters* is not — a model trained on a diverged genome can bias gene calls — so reuse
must be conditioned on measured ANI, not just a shared `SPECIES` string.

**Decided: the ab-initio representative is selected independently and does not need to
match RNASEQ_PREPARE's Trinity representative** (§6 decision 2).

## 3. ANI framework: what's already there, what's missing

`compare_ANI.nf` computes all-vs-all skani ANI within taxonomic groups (default
`--compare GENUS`), unions genomes into clusters via `bin/report_ani.py` at
`--cluster-threshold 95.0`/`--outlier-threshold 90.0`, and writes results into `ani.db`
(SQLite) with `query_species`/`ref_species`/`query_strain`/`ref_strain` columns populated
per pair.

1. **SPECIES-level `--compare` rank — decided to add it** (§6 decision 6). `validRanks` in
   `compare_ANI.nf`/`query_ANI.nf` currently stops at GENUS. For genera not yet covered by
   any GENUS-level run (confirmed live gap: **Beauveria has zero rows in `ani.db`, no output
   directory at all**, despite 109 `Beauveria_bassiana` BUSCO summaries already on disk), a
   SPECIES-scoped run avoids computing the full genus all-vs-all just to get one species'
   pairs. For genera already covered (confirmed live: **Aspergillus fumigatus has 69,751
   within-species pairs already computed**), keep reading `ani.db` filtered to
   `query_species = ref_species` — no recompute needed.
2. **No machine-readable cluster-membership output** — `report_ani.py`'s union-find only
   emits a text report. `species_reuse_clusters.py` (§4.1) is the new machine-readable
   output.
3. **`ani.db` species/strain columns — verified reliable, not a live risk.** A 2026-07-19
   bug (`.living/learnings.md`) could have left these columns blank; queried the live
   production `ani.db` directly — 263,560 total rows, 0 blank `query_species`. Treat as a
   *fragile invariant*: `species_reuse_clusters.py` should assert non-blank species columns
   on load and fail loudly rather than silently misclassifying every strain if this ever
   regresses.

## 4. Proposed pipeline shape

### 4.1 New script: `species_reuse_clusters.py` (`nextflow/bin/`)

Does double duty: (a) computes 99%-ANI reuse clusters + representative for every species,
and (b) **backfills** the shared-parameter store from already-completed legacy PREDICT runs
(no re-run of existing strains required).

Inputs: `ani.db` (SPECIES-scoped runs where available, GENUS-scoped filtered to
same-species pairs elsewhere) + BUSCO completeness from `results/genome_stats/BUSCO_genome/
*.BUSCO_summary.fungi_odb12.txt` (5,930 files already on disk; per-genome, not gated on the
stale merged `asm_stats.tsv.gz` per T-001) + `samples.csv`.

**Representative selection is unconditional** (runs for every species regardless of ANI
coverage — BUSCO/N50/tiebreak don't require ANI data). **Reuse eligibility is ANI-gated**
per strain. These are deliberately decoupled: a representative can be picked and its
ab-initio store built (§4.2) even for a species not yet ANI-covered (representative
selection alone), but individual strains only become `reuse_eligible` once their ANI pair
to that representative is confirmed ≥99%. Concretely, for Beauveria bassiana: representative
selection and store backfill can happen today (BUSCO already covers it), but
`reuse_eligible` assignments for its other strains wait on the pending SPECIES-level ANI
run (§3.1).

For every species with N ≥ 2 strains passing the ANI gate (**no minimum cluster size** — §6
decision 3):
- Build the within-species ANI subgraph (`query_species == ref_species`).
- Pick the representative: **highest BUSCO completeness, N50 as tiebreaker, alphabetical
  `out` as final deterministic tiebreaker** (§6 decision 2).
- A strain qualifies for reuse iff ANI to representative ≥ **99.0%** (§6 decision 1),
  applied uniformly to **all three predictors — AUGUSTUS, SNAP, and GeneMark-ES** (§6
  decision 4, revised from an earlier GeneMark-ES-only proposal once the stage-timing data
  showed Augustus/SNAP's BUSCO-seeding dependency is real and shared — see §1).
- **Fail-closed on missing ANI data**: skani/fastANI/mash omit low-identity pairs entirely
  rather than reporting a low value, so "no pair found" must be treated identically to
  "below threshold" (opt out), never as eligible.
- Strains that fail the gate fall back to independent training — never block, just opt out.
- **No retroactive invalidation.** Existing/legacy strains are never forced to re-predict.
  Applies prospectively to newly-scheduled PREDICT runs; post-hoc comparison against the
  legacy independently-trained baseline is the validation strategy (§5).

Output: `abinitio_reuse_assignments.csv` — columns `species, out, is_representative,
representative_out, ani_to_representative, reuse_eligible`. Doubles as the
**provenance/QC record** (§6 decision 7) — no separate tracking file needed.

### 4.2 Central shared-parameter store

`${params.target}/_shared_abinitio/<species_tag>/`:
```
├── parameters.json                       # rewritten to point at the paths below
├── provenance.json                       # source strain (out), BUSCO score, N50,
│                                          #   ANI cluster info, date captured
├── augustus/species/<species_tag>/       # copied from representative's ab_initio_parameters/
├── <species_tag>.genemark.mod
└── <species_tag>.snap.hmm
```
**Copied, not symlinked**, from the representative's own `predict_misc/ab_initio_parameters/`
— a symlink would break if the representative's `predict_misc/` is later wiped and rebuilt
by the existing staleness logic (`rm -rf predict_results predict_misc`, funannotate.nf:1358).
Component-partial: a representative may legitimately lack a `genemark.mod` if
`--auto-skip-genemark` (already passed at funannotate.nf:1424) zeroed GeneMark's weight for
a fragmented assembly (verified against `predict.py:1733-1758`: fires when longest scaffold
< 50kb) — JSON generation must build from whichever components actually exist.

`species_reuse_clusters.py` populates this store both for the forward pipeline and
retroactively for species whose representative already has a completed legacy PREDICT —
e.g. Aspergillus fumigatus can have its store built today, independent of any Nextflow
wiring changes.

### 4.3 Nextflow wiring changes in `funannotate.nf`

- Load `abinitio_reuse_assignments.csv`, join into the predict channel by `out` (mirrors the
  existing `species_tag` join pattern at funannotate.nf:2092–2137).
- Branch `predict_ch` into `representative` (today's full PREDICT, unchanged),
  `reuse_eligible` (adds `-p <shared json>` to the command at funannotate.nf:1418–1424;
  everything else unchanged), `independent` (today's behavior, unchanged).
- **Two distinct join paths are required for `reuse_eligible` — this was a real gap caught
  in Fable's second review, not just a risk.** The originally-described single mechanism
  ("wait on the representative's `.predict.done` marker") only works when the representative
  is predicted *fresh within the same run*. For the two pilot species specifically
  (Aspergillus fumigatus, Beauveria bassiana), the representative was **already predicted in
  a past run** — there is no live channel item to join against. Two paths:
  1. **In-run**: representative predicted fresh this run → join on its `.predict.done`
     marker per species_tag (as originally described).
  2. **Already-backfilled**: `_shared_abinitio/<species>/parameters.json` already exists on
     disk (from `species_reuse_clusters.py`'s backfill, §4.1) → gate on file-existence via a
     `.filter{}` reading disk state, same idiom already used for `gbkResult`/`staleRnaseq`
     (funannotate.nf:2184).
- **Runtime re-verification, not just an offline CSV.** `abinitio_reuse_assignments.csv` is
  computed offline, before the representative's PREDICT necessarily exists — so for a
  genuinely new species, `reuse_eligible=true` at CSV-generation time is a prediction, not a
  guaranteed fact (the representative's PREDICT could still be skipped, see next bullet).
  The `reuse_eligible` process must re-check at runtime that `_shared_abinitio/<species>/
  parameters.json` (or at least one component) actually exists before adding `-p`, falling
  back to independent training in-process if not.
- **Staleness must account for the representative changing**: a `reuse_eligible` strain's
  GBK must be considered stale if the representative's shared `parameters.json` is newer
  than its own GBK, in addition to the existing rnaseq/trinity staleness check
  (funannotate.nf:1337–1359).
- **Representative failure/skip fallback**: if the representative's own PREDICT is skipped
  by the too-small/fragmented preflight guard (funannotate.nf:1389–1416) or the "not enough
  training models" check (funannotate.nf:1431–1442), `reuse_eligible` strains for that
  species fall back to independent training — covered by the same runtime file-existence
  check above (no `parameters.json` → no `-p` → independent path).

### 4.4 First-N validation gate (integrated, not purely post-hoc)

For a species newly enabled for reuse, PREDICT is scheduled for the representative + first
**N** reuse-eligible strains only; their proteins are fed into the existing `BUSCO_PEP`
process (`BFD.nf:961`, `storeDir "${params.genome_stats_outdir}/BUSCO_protein"`,
`storeDir`-deduplicated so this never double-computes against the pipeline's normal
BUSCO-protein QC step). Only after those N pass a completeness threshold does the pipeline
schedule the remaining strains in that species' cluster. This makes §6 decision-C's
"first N checked before continuing" self-enforcing rather than a manual habit — with no
gate, all strains would predict before anyone noticed a bad representative, defeating the
purpose. Beyond the first-N checkpoint, ongoing monitoring stays a **decoupled, periodic
post-hoc comparison** (§5) — not a live per-strain gate for every subsequent strain.

### 4.5 Rollout controls

- `params.share_abinitio_params` (default `false`) — opt-in.
- `params.force_independent_species` (list) — escape hatch for named species regardless of
  ANI clustering.
- First-N integrated BUSCO-protein gate (§4.4) is the automated backstop; beyond that, no
  additional hard pre-rollout gate (§6 decision 8) — validate against the two pilot species
  as the feature is used, expand species-by-species as confidence builds.

## 5. Post-hoc validation approach

For each pilot species (Aspergillus fumigatus first — already ANI+BUSCO covered; Beauveria
bassiana once its SPECIES-level ANI run completes), compare shared-parameter strains against
their existing independently-trained legacy annotations:
- BUSCO completeness of resulting proteome
- gene/exon/intron count distributions
- wall-clock/CPU-hour delta

**No repeat-content/GC-skew gate built in** (§6 decision 5) — whole-genome ANI is blind to
local GC-skew (GeneMark-ES self-training sensitivity) and repeat/TE-content differences
(AUGUSTUS intron-length/gene-density priors), both plausible failure modes in fungi even at
99% ANI. Check for these as correlates in post-hoc results; revisit only if the data shows
it matters.

## 6. Decisions (settled via user interview, 2026-07-23)

1. **Reuse ANI threshold: 99.0%.**
2. **Representative selection: BUSCO completeness, N50 tiebreak, alphabetical-`out` final
   tiebreak.** Allowed to diverge from RNASEQ_PREPARE's Trinity representative; no
   retroactive redo of legacy annotations.
3. **No minimum cluster size** — reuse applies for any species with ≥2 ANI-qualified
   strains.
4. **Scope: share AUGUSTUS + SNAP + GeneMark-ES uniformly at the single 99% gate.** Revised
   mid-session from an initial GeneMark-ES-only proposal: real stage-timing data (§1) showed
   Augustus's BUSCO-seeding path is expensive (median 26 min, 67% of sampled strains) and
   SNAP training depends on the same BUSCO-derived models when no RNA-seq evidence exists —
   sharing all three avoids running the expensive BUSCO-seed step at all for
   non-representative, non-RNA-seq strains. GlimmerHMM excluded — confirmed disabled in this
   pipeline (§1), produces nothing to share.
5. **No repeat-content/GC-skew gate for now** — check post-hoc only.
6. **Add SPECIES as a valid `--compare` rank** in `compare_ANI.nf`/`query_ANI.nf` for
   future/new ANI runs (uncovered genera like Beauveria); keep reading existing GENUS-level
   `ani.db` (filtered to same-species pairs) for already-covered genera (e.g. Aspergillus).
7. **QC/provenance tracking: a column in `abinitio_reuse_assignments.csv`** — not a separate
   marker file.
8. **Rollout: opt-in flag + escape hatch + integrated first-N BUSCO-protein gate (§4.4);
   no additional hard pre-rollout validation pilot** beyond that — validate against the two
   pilot species as the feature is used, expand as confidence builds.
9. **First-N BUSCO check is an integrated Nextflow gate, not a decoupled post-hoc habit**:
   schedule representative + first N reuse-eligible strains, run existing `BUSCO_PEP`
   (deduplicated via its existing `storeDir`), gate remaining strains in that cluster on a
   completeness threshold. Ongoing monitoring beyond the first N stays decoupled/periodic.
10. **Pilot species: Aspergillus fumigatus (ANI+BUSCO already fully covered) and Beauveria
    bassiana (BUSCO covered, ANI run pending — genus not yet in `ani.db`).**

Residual risk accepted knowingly (raised by Fable's review, addressed by decisions 8/9 above
rather than eliminated): uniform three-predictor sharing + no minimum cluster size means a
single bad representative could still degrade an entire cluster if it slips past the first-N
check but manifests on a later, unusual strain (e.g. one with atypical repeat content).
Decision 5 explicitly defers building a real-time gate for that; the first-N gate (decision
9) and periodic post-hoc monitoring (§5) are the accepted mitigation, not a guarantee.

**Deferred, explicitly out of scope for this work**: refactoring `funannotate.nf` (2,500+
lines, no `modules/`/`subworkflows/` separation) toward nf-core-style structure. Raised
during this design session; user decided to hold it until the ab-initio-reuse work ships,
to avoid touching the same processes for two unrelated reasons at once. Revisit via the
`request-refactor-plan` skill afterward if still wanted.

## 7. Non-goals

- Not touching PASA/transcript evidence handling — stays fully per-strain.
- Not touching Trinity-GG sharing — already correctly done via `RNASEQ_PREPARE`/`species_tag`.
- Not building a repeat-content/GC-skew reuse gate now (§6 decision 5).
- Not building a hard pre-rollout validation pilot beyond the first-N gate (§6 decisions 8/9).
- Not refactoring `funannotate.nf`'s structure (see §6 residual note) — separate initiative.

## 8. Next steps

1. **DONE (2026-07-23)**: `nextflow/bin/species_reuse_clusters.py` written and run against
   Aspergillus fumigatus. Representative: `Aspergillus_fumigatus_Z5` (BUSCO 99.3%, picked
   from ANI-covered candidates only — see bugfix below). 373/375 non-representative strains
   `reuse_eligible` (ANI 99.14-99.93%), 2 fail-closed (no ANI pair). Store backfilled at
   `genome_annotation/_shared_abinitio/Aspergillus_fumigatus/` (augustus species dir +
   genemark.mod + snap.hmm + parameters.json with absolute paths + provenance.json).
   Assignments at `genome_annotation/_reuse_assignments/abinitio_reuse_assignments.Aspergillus_fumigatus.csv`.
   **Bugfix discovered during testing** (logged to `.living/learnings.md`): naive
   highest-BUSCO representative selection can pick a genome with zero `ani.db` coverage
   (added after the last ANI run), which orphans the entire cluster (0/375 eligible in the
   first attempt) — fixed by restricting representative candidates to those with >=1 ANI
   pair before ranking by BUSCO/N50.
2. **IN PROGRESS**: SPECIES-scoped ANI comparison for Beauveria — infra ready (`SPECIES`
   rank added to `compare_ANI.nf`/`query_ANI.nf`, §3.1); real production run pending
   (only stub-tested so far, and `COMBINE_ANI_TABLE` itself needed a bugfix along the way
   — see `.living/learnings.md` 2026-07-23/24 and `todo/TODO_REGISTRY.md` T-008/T-009).
   Once the real run completes, repeat step 1 for it
   (`species_reuse_clusters.py --species "Beauveria bassiana"`).
3. **DONE (2026-07-24, corrected 2026-07-26)**: `funannotate.nf` wiring for the
   **already-backfilled** join path (§4.3) — `share_abinitio_params`/`abinitio_reuse_csv`/
   `force_independent_species` params, `loadAbinitioReuseMap()`/`sharedParamsJsonFor()`/
   `staleSharedParams()` helpers, `-p <parameters.json>` added to the `FUNANNOTATE_PREDICT`
   command, provenance marker (`${out}.predict.abinitio_reused`), staleness extended to
   shared-params mtime.
   **Correction**: this was previously logged as "not yet validated end-to-end", but the
   12-strain T-004/T-013 comparison batch (§5, `.living/findings/funannotate-abinitio-reuse-validation.md`,
   F-006) *did* exercise this exact code path for real — `share_abinitio_params: true` +
   the per-species CSV via `-params-file`, real `nextflow run funannotate.nf` invocations
   (not standalone `funannotate predict`), all 12 strains produced a genuine
   `.predict.abinitio_reused` marker. That satisfies step 6 below; it was just never
   connected back to this doc. This work had also gone uncommitted (stashed) across a
   branch switch to `refactor_modules` — recovered and landed on `feat/abinitio-reuse-wiring`
   (based on `main`) 2026-07-26; re-lint + `-preview` re-verified for both pilot species
   post-recovery (373/376 Aspergillus fumigatus, 32/224 Beauveria bassiana reuse_eligible
   loaded correctly, zero errors). **Still not implemented**: the in-run join path (deferred
   — not needed for either current pilot, whose representatives were already predicted in
   past runs; see T-011).
4. Wire the first-N `BUSCO_PEP` integrated gate (§4.4) into the `reuse_eligible` branch
   (T-010 — deferred as a separate change from step 3 to keep that wiring reviewable).
5. Decide a refresh cadence for `species_reuse_clusters.py` (re-run as `samples.csv`/`ani.db`
   grow) — not yet defined; currently a one-off invocation with no described staleness
   handling of its own. Also **commit `nextflow/bin/species_reuse_clusters.py` and
   `nextflow/bin/combine_ani_table.py`** — both are still untracked in git despite being
   in real production use.
6. **DONE (2026-07-25)**: validation pilot (§5) — 12 strains (6 *A. fumigatus*, 6
   *B. bassiana*), BUSCO completeness within noise (mean Δ -0.13pp) of independently-trained
   legacy baselines, no quality regression attributable to `-p` reuse. See F-006.
7. **Gap found during branch-recovery testing (2026-07-26)**: the fail-closed/independent
   fallback branch (representative + no-ANI-pair strains) has only been exercised
   indirectly — those 3 Aspergillus fumigatus rows (`Z5`, `F2`, `SPS-2`) already had legacy
   `predict_results` pre-dating this feature, so the `predict_ch` filter correctly skipped
   them as "already done" rather than actually running the independent-training branch
   through this code. Low risk (that branch is the unmodified pre-existing PREDICT command,
   gated by an empty bash array), but not directly proven. Exercise it for real the next
   time a genuinely new, non-eligible strain shows up in either pilot species.
8. Also default `abinitio_reuse_csv` (`nextflow/conf/profile_funannotate.config`) points at
   a combined `abinitio_reuse_assignments.csv` that has never been generated — only the
   per-species files exist (`abinitio_reuse_assignments.{Aspergillus_fumigatus,Beauveria_bassiana}.csv`).
   Every real invocation so far has overridden the path explicitly; either teach
   `species_reuse_clusters.py` to also write/append the combined file, or change the default
   to glob per-species files.
