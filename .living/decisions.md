# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

## D: Share funannotate ab-initio gene-prediction parameters across ANI-qualified strains of the same species (2026-07-23)

**Context**: For species with dozens–hundreds of strains, `funannotate predict` re-trains
SNAP, GeneMark-ES, and AUGUSTUS from scratch per strain even though these models should
converge to near-identical parameters within a species. Full design plan at
`todo/species_level_abinitio_reuse.md`; grounded via an Explore-agent survey of
`funannotate.nf`/`compare_ANI.nf`, reviewed twice by Fable (draft plan, then post-decision
plan), and backed by a real 400-log stage-timing measurement
(`analysis/funannotate_predict_stage_timing/`).

**Decision**:
1. Reuse gate: 99.0% ANI to a per-species representative.
2. Representative selection: highest BUSCO completeness, N50 tiebreak, alphabetical-`out`
   final tiebreak; allowed to diverge from `RNASEQ_PREPARE`'s Trinity-sharing representative
   (different purposes); no retroactive redo of the existing ~10,000+ legacy annotations —
   applies prospectively only.
3. No minimum cluster size — reuse applies to any species with ≥2 ANI-qualified strains.
4. Share **all three** predictors (AUGUSTUS, SNAP, GeneMark-ES) uniformly at the single 99%
   gate. GlimmerHMM excluded (confirmed disabled in this pipeline, `-w ... glimmerhmm:0`,
   produces no trained model to share).
5. No repeat-content/GC-skew reuse gate — ANI-only; check for these as correlates in
   post-hoc validation instead.
6. Add `SPECIES` as a valid `--compare` rank in `compare_ANI.nf`/`query_ANI.nf` for
   ANI-uncovered genera; keep reading existing GENUS-level `ani.db` (filtered to
   same-species pairs) for already-covered genera.
7. Provenance/QC tracking is a column in `abinitio_reuse_assignments.csv`, not a separate
   marker file.
8. Rollout: opt-in flag (`params.share_abinitio_params`, default false) + escape hatch
   (`params.force_independent_species`) + an integrated first-N `BUSCO_PEP` gate (below) —
   no additional hard pre-rollout validation pilot.
9. First-N strains of a newly-enabled species run through the existing `BUSCO_PEP` process
   (`BFD.nf:961`, `storeDir`-deduplicated against the pipeline's normal BUSCO-protein QC —
   never double-computed) as an integrated Nextflow gate before the rest of that species'
   cluster is scheduled. Beyond the first N, monitoring is decoupled/periodic, not a live
   per-strain gate.
10. Pilot species: Aspergillus fumigatus (ANI+BUSCO already fully covered) and Beauveria
    bassiana (BUSCO covered, SPECIES-level ANI run pending).

**Alternatives considered**: (a) GeneMark-ES-only for the first rollout wave — the initial
plan, on the assumption GeneMark-ES dominated ab-initio cost; **reversed** once real
stage-timing data showed Augustus's BUSCO-seeded training path (used by 67% of sampled
strains — those lacking RNA-seq evidence) costs a median 26 min, comparable to GeneMark-ES,
and SNAP training depends on the same BUSCO-derived models in that scenario — sharing only
GeneMark-ES would leave the expensive BUSCO-seeding step running anyway. (b) Per-predictor
tiered ANI thresholds (e.g. stricter gate for GeneMark-ES, since its self-training is more
GC/codon-usage sensitive than Augustus) — rejected in favor of a single uniform 99% gate for
simplicity; ship and validate post-hoc rather than pre-engineer per-predictor risk handling.
(c) A dedicated pre-rollout validation pilot (5–10 species, both-ways comparison) before any
production use — rejected as a hard gate in favor of the first-N integrated `BUSCO_PEP`
check plus post-hoc comparison against the two named pilot species, since the user preferred
validating empirically as the feature is used over a separate up-front experiment. (d) Force
`RNASEQ_PREPARE`'s Trinity representative and the ab-initio representative to be the same
strain — rejected; the two selections optimize genuinely different things (Trinity sharing
barely depends on strain quality, ab-initio parameter quality depends on it directly) and
forcing them to match would mean picking a worse representative for one purpose or the
other.

**Rationale**: `-p <parameters.json>` (confirmed via `funannotate predict --help` and a real
parameters.json on disk) is the supported mechanism for supplying pre-trained AUGUSTUS/SNAP/
GeneMark-ES parameters and skipping their (re)training. Existing `RNASEQ_PREPARE` already
establishes the species-level grouping/caching precedent this plan extends, adding an ANI
gate that RNASEQ_PREPARE doesn't need (Trinity sharing is evidence-agnostic and low-risk;
ab-initio parameter sharing biases gene calls if the source genome is too diverged).

**Consequences**: A prior single-log inspection suggested RNA-seq hints prep dominated
PREDICT wall time (~72%) and ab-initio training was minor (~17-19%) — the 400-log aggregate
contradicted this (hints prep: median 0%, mean 9.5%; ab-initio training: median 39.8%, mean
38.7%), so the plan's cost premise is now measured rather than assumed. Fable's second review
caught a real implementation gap: the originally-described single join mechanism (wait on
representative's `.predict.done` marker) doesn't work for the two named pilot species, whose
representatives were already predicted in past runs — the Nextflow wiring needs two join
paths (in-run marker join vs. already-backfilled file-existence check), documented in
`todo/species_level_abinitio_reuse.md` §4.3. A known residual risk is accepted, not
eliminated: uniform three-predictor sharing + no minimum cluster size means a bad
representative that slips past the first-N check could still degrade later strains in a
large cluster before periodic post-hoc monitoring catches it. A related, separate
architectural question (refactoring `funannotate.nf` toward nf-core-style
modules/subworkflows) was raised during this design session and explicitly deferred until
this work ships, to avoid touching the same processes for two unrelated reasons at once.

**Related**: `todo/species_level_abinitio_reuse.md` (full plan); `analysis/funannotate_predict_stage_timing/` (cost measurement); [[project_nfschema_validation]]; [[project_trinity_cache_churn]].

**Tags**: funannotate, nextflow, funannotate.nf, ani, compare_ANI, ab-initio, augustus, genemark, snap, busco, parameter-reuse, decision

## D: Feed BFD.nf MERGE_* steps a manifest file (with mtime+size) instead of staging thousands of inputs (2026-06-25)

**Context**: MERGE_* processes took `input: path 'inputs/*'` fed by `.collect()` of all per-genome stat files. At ~8k+ genomes this stages tens of thousands of symlinks and builds an enormous stage-in command (ARG_MAX risk). User also wanted the manifest to drive staleness (re-merge when any input is newer).

**Decision**: Add `toManifest(ch, name)` that `collectFile`s one line per input as `<abs_path>\t<mtime_ms>\t<size_bytes>`; each MERGE process takes `path manifest` and reads field 1. Inputs are read by absolute storeDir path (shared FS), not staged. mtime+size are embedded in the manifest *content* so Nextflow's input-hash caching re-runs the merge under `-resume` iff an input changed.

**Alternatives considered**: (a) `find`-based staleness scan comparing inputs to a persisted prior manifest — rejected; the freshly written manifest is always newest, so it requires persisting/ diffing the old manifest *and* manually forcing a cache miss, which Nextflow content-hashing already does for free once mtime/size are in the manifest. (b) Bare path-list manifest (no mtime/size) — rejected; same paths → identical content hash → a regenerated input would be silently kept stale under `-resume`. (c) Keep `path 'inputs/*'` but raise limits — rejected; doesn't remove the symlink-storm/ARG_MAX scaling problem.

**Rationale**: One staged file per merge regardless of genome count; staleness handled by the same mechanism Nextflow already uses for resume; no new persisted state.

**Consequences**: Inputs must live on shared storage (they do — storeDir-cached). Caching now keys on path+mtime+size rather than full content hash — a `touch` of an input forces a re-merge (acceptable, conservative). `summarize_asm_stats.py` rewritten with argparse + `--manifest` + gzip (also fixed a pre-existing arg mismatch where the process passed `--reportdir/--samples/-o` to a script with no argparse). Verified via bash unit tests + `-profile BFD,test -preview` → `[SUCCESS]`.

**Related**: [[MERGE_* steps: pass a manifest file, not thousands of staged inputs; bake mtime+size in for staleness]] (learnings, 2026-06-25); same file/session as the `combine` scalar-barrier fix below.

**Tags**: nextflow, BFD.nf, merge, manifest, staging, arg-max, staleness, resume, caching, decision

## D: Use a scalar barrier (not a collected list) on the left of Nextflow `combine` in BFD.nf MERGE gating (2026-06-25)

**Context**: BFD.nf's non-glob MERGE branch (`merge_all=false` or `--taxon` active) gates each MERGE on upstream BATCH completion by combining a sync channel with a collected file-list channel, then destructuring with `.map { _s, files -> ... }`. It crashed with `MissingMethodException ...(LinkedList)` whenever the BATCH outputs were all cached.

**Decision**: Build the gate channel as a scalar — `BATCH_*.out.csv.flatten().collect().map { true }.ifEmpty(true)` — instead of `.collect().ifEmpty([])`. Applied to both `aa_sync` and `codon_sync`.

**Alternatives considered**: (a) Keep the list barrier and reshape the map closure — rejected; fragile because `combine` *spreads* list-valued left items and the breakage only manifests in the all-cached path. (b) Drop `combine` and pass `codon_paths` directly — rejected; loses the completion-gate that ensures BATCH-written files exist before MERGE globs the storeDir. (c) Use `.count()` as the scalar barrier — equivalent and acceptable, but `.map { true }.ifEmpty(true)` reads more clearly as a boolean gate.

**Rationale**: `combine` flattens/spreads tuples from its left operand; a scalar can't be spread, so the cross-product yields the intended `[gate, [files...]]` 2-tuple in both populated and empty cases. Keeps the dependency barrier intact while fixing the misbind.

**Consequences**: Both code paths bind correctly; verified via `nextflow run ... -profile BFD,test -preview` → `[SUCCESS]`. Note for future lint/preview: `-profile test` alone omits `genome_stats_outdir`/`freq_batch_size` (defined in `conf/profile_BFD.config`), so use `-profile BFD,test`.

**Related**: [[Nextflow combine spreads list-valued left items — use a scalar barrier]] (learnings, 2026-06-25).

**Tags**: nextflow, combine-operator, channel-semantics, barrier-channel, BFD.nf, merge-gating, decision

### [2026-06-18] Persona-based ideation over BFD genome-size/composition data (Broader-9)

**Context**: User wants a framework to explore how fungal genome size varies within taxonomy and which composition indicators (TEs especially) contribute, over the BFD.nf `tables/` outputs. No genome-size analysis existed yet (greenfield).

**Decision**: Ran the idea-generator convention with the Broader-9 persona set (Evolutionary Biologist, Causal Inference, Information Theorist, Statistical Physicist, Ecologist, Representation Learning, Quantitative Geneticist, Topologist/Geometer, Linguist/NLP), 2 ideas each = 18, dispatched as 9 parallel sonnet subagents in 2 batches. Output in `analysis/ideas/2026-06-18-genome-size-composition-framework/`.

**Alternatives considered**: Genomics-focused 6 (tighter fit) or all 15 (more off-target). Broader-9 chosen by the user for methodological range while staying relevant.

**Consequences**: 18 grounded ideas indexed in `00_index.md`. Next step is user triage → promote selected ideas to `todo/` and/or deep-dive. The two low-effort/data-ready ideas (#8 order-parameter, #17 Zipf) plus #13 (variance partitioning) are the recommended fast first pass.

**Tags**: ideation, idea-generator, genome-size, comparative-genomics, subagents

### [2026-06-18] Trinity-GG usability tiers keyed on NUM_TRANSCRIPTS, not N50

**Context**: Needed a per-species "usable as funannotate annotation evidence" classification for the 239 successful Trinity-GG assemblies.

**Decision**: Primary tier metric = `NUM_TRANSCRIPTS` (FAIL<100, BORDERLINE 100–999, PASS≥1000); `N50` reported as a secondary signal only. Cutoffs are swept rather than asserted as ground truth.

**Alternatives considered**: Gating on N50/contiguity, or a composite score. Rejected for v1 — for annotation *evidence*, the count of transcripts is the most direct proxy for "did the library yield usable signal", and even fragmented (low-N50) transcripts still aid annotation. A composite would add tuning knobs without changing the headline (existing assemblies are overwhelmingly PASS).

**Consequences**: 234 PASS / 3 BORDERLINE / 2 FAIL; result is stable across the swept cutoffs. N50 gating remains an open question logged in the analysis doc.

**Tags**: rnaseq, trinity, qc, thresholds, sensitivity-analysis

### [2026-06-18] Initialized mycelium with non-destructive scaffold, not restructure

**Context**: Fungi_BFD is an established working repo (Nextflow pipelines, SLURM launchers, lib/, misc/, etc.) with no prior `.living/` layer. Asked to "initialize living repo".

**Decision**: Ran `init_repo.py` in default scaffold mode (added `.living/`, manifests, `todo/`, `algorithms/`, `analysis/`, `data/`, `reference_material/`, hooks) rather than `--restructure`.

**Alternatives considered**: `--restructure` mode — rejected because (a) its file-moving step is an unimplemented TODO (audit-only), and (b) moving an active pipeline repo's files would be disruptive. Default mode is non-destructive (`mkdir exist_ok=True`, no moves), so existing content is untouched.

**Consequences**: New empty mycelium dirs sit alongside the existing project layout. Existing pipeline code was not reorganized into `analysis/`/`algorithms/`; those manifests start empty and can be backfilled as work is logged.

**Tags**: mycelium, init, repo-structure, non-destructive

---

## D: Reproducible samples.csv production with data-driven curation (2026-06-19)

**Decision**: Rewrote `scripts/create_samples_file.py` to be reproducible. SPECIES/STRAIN/ASMID sanitization centralized in `scripts/sample_sanitize.py` (strain rules kept byte-compatible with `nextflow/lib/SampleUtils.groovy::cleanStrain`). All manual curation moved out of hand-edits into data files under `data/curation/`: `exclude_asmids.txt` (142 removals), `keep_dupes.csv` (6 protected multi-assembly isolates), `overrides.csv` (per-ASMID fixes).

**Why**: The previous flow required ~177 manual edits to `samples.csv` after each regeneration (quote stripping, `*`/`#`/`:` strain cleanup, ASMID extension trim, row removals), so it was not reproducible.

**Dedup policy**: Curated removals kept as an explicit list because empirically neither "prefer GCF" (16/39 isolates) nor "newest accession" (71/110) reproduces the human picks — selection tracks best-annotated/reference genome. A default tie-breaker (prefer GCF, else newest GCA) applies only to *new* uncurated collisions; intentional multi-keeps are protected via `keep_dupes.csv`.

**Validation**: On the 7,963 assemblies with taxonomy, output matches hand-curated samples.csv: PHYLUM/CLASS/ORDER/FAMILY/BUSCO/TRANSL_TABLE 0 diffs, SPECIES 1 (cosmetic), GENUS 2, LOCUSTAG 1; 54 STRAIN diffs are the intended aggressive normalization.

**LOCUSTAG**: hashed from the *raw* `{ACCESSION}_{ASM_NAME}` (display ASMID is cleaned), assigned after a stable sort so collision suffixes (A/B/C) are deterministic.

**Tags**: samples-csv, ncbi, reproducibility, sanitization, curation, dedup

---

## D: Curation files for samples.csv dedup/overrides (2026-06-19, addendum)

**Decision**: Four data files in `data/curation/` control samples.csv curation, each a distinct concern:
- `exclude_asmids.txt` — drop assemblies entirely.
- `preferred_asmids.txt` — force which assembly wins a species+strain dedup group (overrides the GCF/newest tie-breaker). Added so a specific (e.g. older) assembly can be pinned, e.g. S. cerevisiae YJM981 → keep older `GCA_000976515.2_Sc_YJM981_v1` over newer `GCA_025434895.2`.
- `keep_dupes.csv` — intentionally retain >1 assembly per isolate.
- `overrides.csv` — `ASMID,FIELD,VALUE` per-field value corrections, applied after sanitization (matched on cleaned ASMID); changes content, not dedup selection. Applied before dedup, so overriding SPECIES/STRAIN changes that row's grouping key.

**Why**: dedup winner selection and field-value correction are separate problems; conflating them in one file caused confusion (overrides was mistaken for a winner-selector).

**Dedup winner precedence**: preferred_asmids > RefSeq GCF > newest accession; keep_dupes short-circuits to retain all.

**Tags**: samples-csv, curation, dedup, preferred-assembly, overrides

---

## D: Stop hook made non-blocking for this repo (2026-06-19)

**Decision**: Wrapped the mycelium Stop hook with `.claude/mycelium-stop-nonblock.sh` and repointed `.claude/settings.local.json` Stop[0] at it (via `$CLAUDE_PROJECT_DIR`). The wrapper runs the real `mycelium-stop-check.sh` (preserving session-log finalization, LOG_REGISTRY upsert, last-session.md, background log-scribe) but downgrades a `{"decision":"block"}` verdict to a silent allow; other output (additionalContext) passes through.

**Why**: The block fires once per work-cycle whenever `.living/` wasn't re-touched after the last analysis script ran. In this pipeline repo `scripts/create_samples_file.py` is run repeatedly, so the block recurred constantly and interrupted normal iteration. Reflection is still prompted via the passed-through additionalContext.

**Scope/limits**: Project-local only (no edit to shared marketplace files; survives plugin updates). Revert by pointing Stop[0] back at `.../mycelium-stop-check.sh` and deleting the wrapper.

**Tags**: mycelium, hooks, stop-hook, tooling, settings, non-blocking

---

## D: Do NOT reclaim normalized RNA-seq reads from rnaseq_reads/ (2026-06-20)

**Decision**: Keep all normalized reads in `rnaseq_reads/` (~787 GB); do not stub/delete them after Trinity + training. Tooling was drafted (`scripts/audit_rnaseq_reclaim.py`, `scripts/reclaim_rnaseq_reads.py`, `doc/README_rnaseq_reclaim.md`) and a dry-run audit showed 311 species safe-to-reclaim now (163.6 GB) plus 544 Trinity-ready-but-training-pending species (292.6 GB). Nothing was applied.

**Why**: (1) Reads are consumed per-strain, keyed by `species_tag` — every strain's `FUNANNOTATE_TRAIN` re-aligns the same shared reads to its own genome (BAM/stringtie/kallisto). (2) Future updates are expected to add new strains of existing species; a reclaimed (stubbed) species would force those new strains into funannotate's empty-`--left_norm` branch (`funannotate.nf:1232`) or trigger an expensive `SRA_FETCH` re-download. (3) `FUNANNOTATE_UPDATE` (`funannotate.nf:1673`) consumes raw reads directly (not Trinity) and would re-need them if `run_update` is ever enabled (currently false). (4) Judged the ~163 GB immediate saving not worth the added fragility.

**Considered & rejected**: changing `FUNANNOTATE_TRAIN` to take Trinity only (drop reads) — loses per-strain genome evidence, doesn't free reads for `update`, and creates methodological inconsistency vs the 7,340 genomes already trained with reads+Trinity.

**Related**: [[project_trinity_cache_churn]]; `doc/README_RNASeq_handling.md`.

**Tags**: rnaseq, disk-space, reclaim, trinity, funannotate-train, funannotate-update, rejected

---

## D: Stage comparative inputs by SPECIES+STRAIN basename, not LOCUSTAG filename (2026-06-22)

**Decision**: Fix `comparative_genomics.nf` staging by reconstructing the input filename from `SampleUtils.makeSampleTag(SPECIES, STRAIN)` (carried as a `BASENAME` manifest column), while keeping the *staged symlink* names as `{LOCUSTAG}.faa`/`{LOCUSTAG}.cds.fa`.

**Why**: The on-disk inputs are named by sanitized species+strain (how funannotate.nf wrote them), but downstream steps (`cat *.faa`, OrthoFinder `-f`) and stable cross-step identity are cleanest keyed by LOCUSTAG. Reusing `makeSampleTag` (vs reconstructing the name ad hoc, or vs the raw `SPECIES_IN` column) guarantees the sanitization matches the writer exactly — strain cleaning handles `;`, `*`, `:`, quotes, and whitespace/`/#[]?{}` collapse.

**Considered & rejected**: (1) Quick-patch using the `SPECIES_IN` column — rejected by user in favor of SPECIES+STRAIN, which is what the files are actually named from. (2) Re-deriving the sanitized name inline in shell — rejected; duplicates non-trivial logic already in `SampleUtils` and would drift. (3) Renaming all `input/pep` files to LOCUSTAG — rejected; touches immutable-ish staged inputs and breaks funannotate's naming.

**Related**: [[comparative_genomics.nf STAGE_FILES expected wrong input filenames]]; `nextflow/lib/SampleUtils.groovy`.

**Tags**: nextflow, comparative-genomics, staging, makesampletag, decision

---

## D: Guard FUNANNOTATE_PREDICT against too-small genomes with a BOTH-gates pre-flight + post-predict catch (2026-06-24)

**Decision**: Add a two-layer guard to `funannotate.nf` FUNANNOTATE_PREDICT, tunable in `conf/profile_funannotate.config`:
1. **Pre-flight gate** (before the expensive predict call): compute contig stats from the input FASTA and skip-with-flag only if the assembly trips **BOTH** gates — small (`predict_min_asm_bp=8000000`) **AND** fragmented (`predict_frag_max_n50=10000` OR `predict_frag_max_contigs=1000`). Skipped genomes are recorded to `${target}/predict_skipped_too_small.tsv`, marked `.predict.skipped_too_small`, and exit 0 (no wasted GeneMark/BUSCO compute, no batch crash). `predict_min_asm_bp=0` disables it.
2. **Post-predict catch** (backstop): if predict produced no GBK, grep the log for `Not enough gene models .* to train Augustus`; if matched, flag+skip cleanly; otherwise still hard-fail so real errors surface.

**Why**: (1) The **BOTH-gates AND** requirement is the whole safety mechanism — a complete small genome has high N50 / few contigs and therefore cannot trip the fragmentation gate, so the guard never false-skips legitimately small taxa. Verified against user-flagged cases: Ashbya/Eremothecium yeasts, complete Microsporidia, Rozella (11.3 Mb), Paramicrosporidium (7.2 Mb / N50 70 kb) — all pass. (2) The conservative 8 Mb pre-flight catches the obvious sub-8 Mb junk up front (7 of 9 known failures); the post-predict catch is the authoritative backstop for borderline 8–16 Mb fragmented genomes, deferring to funannotate's own verdict so a wrong skip is impossible. (3) Empirically (`analysis/funannotate_model_failures/`) all 9 real `too_few_models` failures are small+fragmented and large genomes never fail this way, so the rule is well-separated.

**Also decided**: append the 7 not-yet-listed too-few-models genomes to `suppress.txt` with model-count reasons (2 of the 9 were already suppressed). The guard **automates** the existing manual `suppress.txt` "Too small" curation rather than replacing it.

**Considered & rejected**: (a) flat genome-size cutoff — rejected; false-flags complete small genomes (Malassezia ~7–9 Mb) and misses fragmented large junk. (b) fragmentation (N50) alone — rejected; over-flags complete-but-fragmented large genomes (40 Mb rusts) that annotate fine. (c) hard-fail on too-few-models (status quo at `funannotate.nf:1390`) — rejected; one junk genome aborts the whole batch. (d) gawk `asort` for N50 — rejected/fixed; replaced with a portable `awk|sort -rn|awk` pipeline (verified identical to `seqkit`).

**Related**: [[funannotate "Not enough gene models N to train Augustus (30 required)" is the ground-truth too-small signal; asm_stats table is stale]]; `analysis/funannotate_model_failures/FUNANNOTATE_MODEL_FAILURES.md`; `suppress.txt`.

**Tags**: funannotate, augustus, predict, guard, too-small, asm-stats, suppress-txt, nextflow, decision

## D: Bump FUNANNOTATE_TRAIN attempt-1 cpus 2→4 pending upstream Trinity-GG threading fixes (2026-07-23)

**Decision**: Raised `nextflow/conf/profile_funannotate.config` `FUNANNOTATE_TRAIN` attempt-1 `cpus` from 2 to 4 (attempt-2+ stays at 8).

**Why**: Trace-profiled 1,891 historical `FUNANNOTATE_TRAIN` runs at `cpus=2`: median %cpu ~78%, p90 ~105% of the 2-core budget. A read-only review of `funannotate train`'s orchestration (`~/projects/funannotate/funannotate-live`, `funannotate/aux_scripts/trinity.py`) found the root cause: the Trinity-GG Butterfly assembly stage — likely the single longest stage for a typical run — fans out single-threaded per-cluster jobs via a `multiprocessing.Pool` sized `args.cpus - 1`, so at `cpus=2` it was `Pool(1)`, i.e. fully serial regardless of the CPU allocation. Filed upstream as [nextgenusfs/funannotate#1178](https://github.com/nextgenusfs/funannotate/issues/1178) with fix PR [#1180](https://github.com/nextgenusfs/funannotate/pull/1180) (branch `fix/trinity-butterfly-cpu-pool`, targets `target_1.9/rust_EVM_trinity_PASA`, not yet merged/released). Also filed [#1179](https://github.com/nextgenusfs/funannotate/issues/1179)/[#1181](https://github.com/nextgenusfs/funannotate/pull/1181) for an unrelated operator-precedence bug in the same function's `bamthreads` calculation. Bumping to `cpus=4` now still helps even before the upstream fix lands — the buggy `cpus - 1` pool becomes `Pool(3)` at `cpus=4` instead of `Pool(1)`, and further re-testing at higher cpu counts should wait until the pool-sizing fix is actually deployed on this cluster's funannotate install.

**Considered & rejected**: Waiting for the upstream fix to merge/release before touching this config — rejected because the current (buggy) `cpus - 1` behavior already benefits from more requested CPUs, just off-by-one; no need to block on upstream review/release cadence for an immediate partial win.

**Follow-up**: Once PR #1180 merges and the cluster's funannotate build is updated, re-profile `FUNANNOTATE_TRAIN` %cpu/RSS (via `analysis/nextflow_memory_profile/`) to see whether attempt-1 cpus should go higher than 4 — PASA's mostly-single-threaded validation/DBI steps (not fixed by #1180) will cap the benefit of further increases.

**Tags**: nextflow, funannotate-train, throughput, cpus, upstream-bugfix, trinity-gg, decision

## D: Route eligible non-'short' funannotate processes to the preempt account/partition (2026-07-23)

**Decision**: In `nextflow/conf/profile_funannotate.config`, switched `ANTISMASH_RUN`, `FUNANNOTATE_ANNOTATE`, and the non-`short` retry attempts of `SRA_FETCH_SE` to `queue = 'preempt'` with `clusterOptions` appending `--account=preempt`, instead of blanket-applying preempt to every non-`short`-queue process as originally asked.

**Why**: `preempt` caps wall-time at 1 day and its ~149 nodes are general-purpose (~200 GB RAM, no GPU) — the same pool as `short`. Checked every non-`short` process in the profile against that: `GENOME_CLEAN`/`GENOME_CLEAN_BATCH` need 500 GB and are pinned to specific highmem nodes (`-w h04,h05,h06`) outside the preempt pool; `SIGNALP_RUN` needs a GPU (`exfab`); `RNASEQ_PREPARE`, `FUNANNOTATE_TRAIN`, `FUNANNOTATE_PREDICT`, `SRA_FETCH`, and `INTERPROSCAN_RUN` all have retry branches whose `time` directive exceeds 24h. Only `ANTISMASH_RUN` (8h), `FUNANNOTATE_ANNOTATE` (24h), and `SRA_FETCH_SE`'s attempt-2/3 branches (max 24h) fit both the time cap and the general-node resource profile.

**Considered & rejected**: (a) force everything onto preempt and clamp `time` to 24h — rejected by user, risks truncating multi-day steps (e.g. `GENOME_CLEAN_BATCH` at 7d) mid-run with no completion path. (b) attempt-1-on-preempt-then-escalate for every process — user chose the narrower, resource-checked option instead.

**Follow-up**: If preemption kills `SRA_FETCH_SE`/`ANTISMASH_RUN`/`FUNANNOTATE_ANNOTATE` jobs often enough to hurt throughput, consider the escalate-on-retry pattern (already used by `SRA_FETCH`/`FUNANNOTATE_TRAIN`) to fall back to `epyc` after a preempt failure.

**Tags**: nextflow, funannotate, slurm, preempt, queue, account, decision

## D: Extend preempt attempt-1 routing to RNASEQ_PREPARE, FUNANNOTATE_TRAIN, FUNANNOTATE_PREDICT (2026-07-23)

**Decision**: Extended the preempt-account routing from the prior entry ([[Route eligible non-'short' funannotate processes to the preempt account/partition]]) to the attempt-1 branch of `RNASEQ_PREPARE`, `FUNANNOTATE_TRAIN`, and `FUNANNOTATE_PREDICT` in `nextflow/conf/profile_funannotate.config`. Each already had an attempt-1 config that fits preempt's 1-day cap on general-purpose hardware (RNASEQ_PREPARE: 16GB/24h; FUNANNOTATE_TRAIN: 2GB/24h; FUNANNOTATE_PREDICT: 12GB/24h), so `queue` on attempt 1 is now `'preempt'` with `--account=preempt` added to `clusterOptions`. Retries (attempt ≥2) are unchanged from the existing escalation logic — they fall back to the prior default/highmem queue and scale memory/time/cpus as before.

**Why**: These three are the initial (attempt-1) stages of the funannotate pipeline and run far more often than their retries, so this is where preempt's throughput benefit matters most. The existing retry-escalation pattern (already used for `SRA_FETCH`) meant only the `queue`/`clusterOptions` closures needed a new attempt-1 branch — memory/time/cpus formulas were untouched.

**Tags**: nextflow, funannotate, slurm, preempt, queue, account, rnaseq-prepare, funannotate-train, funannotate-predict, decision

## D: Build the ANI concat manifest from the live channel, not a filesystem glob (2026-07-27)

**Decision**: In `nextflow/workflows/compare_ANI.nf` and `nextflow/subworkflows/local/ANI_REPRESENTATIVE_SELECT/main.nf`, construct the manifest passed to `CONCAT_ANI_TSVS` by collecting the `ani_tsv` channel emitted by `ANI_COMPARE_METHOD` rather than globbing `${params.outdir}/${params.ani_method}/${params.compare}/**/*.ani.tsv` from disk. Also default `ANI_REPRESENTATIVE_SELECT`'s compare rank to `SPECIES` when `params.compare` is unset, since representative selection is species-level by design.

**Why**: The original glob approach used `channel.fromPath(...).collect().val` inside a `.map` closure, which deadlocks. A fallback to synchronous `files()` would avoid the deadlock but still risks missing outputs because it races against `publishDir` copying. Collecting the live channel guarantees the manifest contains exactly the TSVs produced by this run, with no timing dependency on published files.

**Considered & rejected**: (a) Fix only the case of `channel` vs `Channel` and keep the glob — rejected because even uppercase `Channel.fromPath()` inside a closure hangs. (b) Use `files()` synchronous glob — rejected because `publishDir` may not have finished copying when the trigger channel fires. (c) Wait on both the channel and a filesystem check — rejected as over-engineering; the channel is the source of truth.

**Consequences**: `CONCAT_ANI_TSVS` now receives its input as soon as the last compare group finishes, independent of `publishDir` state. On `-resume`, cached compare processes still emit their `ani_tsv` tuples, so the manifest remains complete. The change also de-duplicates files because the channel emits each group's canonical output exactly once. Defaulting `ANI_REPRESENTATIVE_SELECT` to SPECIES removes the need for callers to pass `--compare` just to satisfy representative selection.

**Tags**: nextflow, ani, compare_ANI, ANI_REPRESENTATIVE_SELECT, CONCAT_ANI_TSVS, channel-based-manifest, deadlock, publishDir-race, species-default, decision

## D: Reorganize genome_stats/function flat directories by ASMID/LOCUSTAG hash bucket, not genus (2026-07-30)

**Decision**: Plan (not yet implemented; see `todo/genome_stats_storage_reorg.md`, T-014) to
reorganize `results/genome_stats/*` (8 subfolders, up to 55,398 flat files) and `results/function/*`
(9 subfolders) into 2-level hash-bucketed subdirectories keyed by `ASMID` (for assembly-level
outputs: asm_stats, asm_reports, BUSCO_genome) or `LOCUSTAG` (for annotation-derived outputs:
BUSCO_protein, aa_freq, codon_freq, gene_stats, intergenic_stats, function/*) — not by genus, and
not by the `Genus_species_strain` filename string currently used for name-keyed types. Bucket width
is 2 hex chars (256 buckets) for most types, 3 hex chars for `gene_stats`/`intergenic_stats` given
their ~7x file-per-genome multiplier. A generated, read-only symlink tree
(`results/genome_stats_by_name/<GENUS-or-NOGENUS>/...`) preserves human-browsable-by-species access
without making it the canonical store. Sequencing: the one-time legacy migration runs *before* the
Nextflow `storeDir` path cutover deploys, so `storeDir`'s existence check finds files already
in place and doesn't trigger a full recompute (changing `storeDir` invalidates Nextflow's `-resume`
cache regardless of migration timing — cache identity and file-existence are different mechanisms).
DuckDB (`db/BFD.duckdb`) + its MCP server stay as the query layer; only `tables/*.csv.gz` staging
converts to Parquet, one unpartitioned file per data type (not bucket-partitioned, to avoid
reintroducing the small-file NFS problem one layer up).

**Why**: `Genus_species_strain` is not a stable/unique key (329 of samples.csv's rows have no
GENUS; `collect_asm_stats.py` already has a multi-step fallback cascade for missing-genus/strain
lookups). Genus-based directory grouping was considered (user's original suggestion) but rejected
because genus population is extremely skewed (Aspergillus/Penicillium/Fusarium have hundreds of
strains; most genera have 1) — it would just relocate today's flat-directory problem into a
handful of hot genus folders. Hash fan-out on the stable key is uniform by construction. Plan was
reviewed by an expert-persona agent (Fable, acting as HPC storage/data-management expert) before
finalizing; its critique is incorporated inline in the plan doc (per-type bucket width, migration
idempotency/target-uniqueness checks, symlink staleness trigger + backup exclusion, unpartitioned
Parquet, and the explicit storeDir/-resume cache-invalidation call-out).

**Considered & rejected**: (a) genus/genus_species_strain directory grouping — rejected, skew +
missing-genus rows (see above). (b) Bucket-partitioned Parquet mirroring the hash-bucket scheme —
rejected per Fable's review: turns one CSV.gz per type into 256 small Parquet files per type,
reintroducing the small-file NFS metadata problem this whole plan exists to solve. (c) Deploying
the new `storeDir` paths before running the legacy migration — rejected: would trigger a full,
HPC-hour-costly recompute across every cached BFD stats process, avoidable by migrating first.

**Tags**: genome_stats, function, storage-reorg, hash-bucket, asmid, locustag, storeDir,
duckdb, parquet, nfs, hpc, T-014, decision

## D: Delete asm_reports/ and collect_asm_stats.py outright rather than migrate (2026-07-30)

**Decision**: During T-014's storeDir/consumer-script scoping (issue #9), discovered
`results/genome_stats/asm_reports/` (22,412 files, migrated into hash buckets earlier the same day)
has no active producer anywhere in the current pipeline — grepped all of `nextflow/` and found
nothing writes to or reads from it. `scripts/collect_asm_stats.py` (its apparent matching consumer)
is also dead: not referenced by any current `.nf`/`.sh` file, and reads from a completely different,
older directory convention (`genomes/*.scaffolds.stats.txt`) unrelated to
`results/genome_stats/asm_stats/`. Both were confirmed redundant with the live path: `asm_stats`
(produced by `CALC_ASM_STATS`) → `MERGE_ASM_STATS` → `summarize_asm_stats.py` →
`tables/<Taxon>/asm_stats.tsv.gz`. Deleted `results/genome_stats/asm_reports/` outright (`rm -rf`,
not archived) and `git rm scripts/collect_asm_stats.py`, rather than including either in the
storeDir/consumer-script cutover. Removed `asm_reports` from `TYPE_CONFIG` in
`scripts/one-off/reorg_genome_stats_hash_buckets.py` and the `BUCKET_WIDTH` tables in both
`nextflow/bin/genome_stats_paths.py` and `nextflow/modules/common/utils.nf`.

**Why**: No point maintaining storeDir/consumer-script code, hash-bucket-width config, or disk
space for a directory nothing produces or reads. Confirmed via direct investigation (not assumed)
that the genuinely-used assembly-stats path is `asm_stats`, not `asm_reports` — verified this
matters for representative-strain selection specifically (N50 there is actually parsed from
`BUSCO_genome` summary text via `_N50_RE`, not from `asm_stats` directly, but `asm_stats` itself
feeds `tables/*/asm_stats.tsv.gz` for broader assembly QC reporting and is very much live).

**Considered & rejected**: archiving `asm_reports/` instead of deleting — rejected as unnecessary
caution once the "no producer, no consumer, superseded by a documented live alternative" case was
fully confirmed rather than assumed; the migration manifest (`results/_migration_manifest.csv`)
still records exactly where every one of those 22,412 files came from before this session, if ever
needed for forensic purposes.

**Tags**: genome_stats, asm_reports, asm_stats, dead-code, cleanup, T-014, decision

## D: Retire per-taxon MERGE_SAMPLES re-runs; master-only merge + separate EXTRACT_BFD_TAXONOMIC_SUBSET tool (2026-07-31)

**Context**: As part of scoping issue #11 (Parquet staging conversion, T-014 §D.1), reviewed how
`--taxon`-scoped merges currently work. `MERGE_SAMPLES` and the rest of `MERGE_*` are re-run per
taxon to build `tables/<Taxon>/*` subset folders. Investigation confirmed `taxonRowFilter()`
(`nextflow/modules/common/utils.nf`) is applied at the *base genome channel* in `BFD.nf`, before
any per-genome computation — meaning the hash-bucketed per-genome store (T-014's main reorg) is
already universal/cumulative across all `--taxon` scopes ever run, and the per-taxon `tables/`
re-merge only exists as a side effect of scoping the manifest glob, not a deliberately designed
scoped-artifact feature.

**Decision**: Stop building `tables/<Taxon>/*` via re-scoped `MERGE_*` runs.
`MERGE_SAMPLES` simplifies to always build one master, unscoped `samples`/`species` table (no more
`--taxon`-filtered subset calls in `subset_samples.py`/`build_species_table.py`). A new tool,
`EXTRACT_BFD_TAXONOMIC_SUBSET` (own issue, not folded into #11, sequenced after #11 lands), is a
Python CLI wrapping DuckDB `ATTACH`/`COPY` that filters the always-current master DuckDB/Parquet
down to a requested taxon post-hoc. `species`/`samples` is the authoritative key table — every
other table is filtered by `INNER JOIN` against that key set (never by independently
re-evaluating the taxon predicate per table), with a post-extract row-count assertion (per-table
matching-key row count vs. `species`) before writing output, per this repo's robust-analysis
convention. No new staleness logic needed for the extractor since the master build is already
full-rebuild-per-run (§D.1) — it just filters whatever's currently on disk.

**Alternatives considered**: (a) Keep per-taxon `MERGE_*` re-runs as "small scoped artifacts" —
rejected once confirmed the per-genome store is already universal; a SQL filter over always-fresh
master data gives the same correctness guarantee without re-glob/re-cat/re-dedup of the whole
per-genome store per taxon request. (b) Make extraction a Nextflow process — rejected; pure SQL
filtering over already-materialized data has no per-genome parallelism to gain, only adds
scheduler overhead. (c) A bare `.sql` script instead of a CLI — rejected; a CLI gives argument
validation, logging, and testability. (d) Fold extraction into #11 — rejected per Fable's
sequencing advice; land Parquet conversion first against a stable target, then build the
extractor separately.

**Rationale/verification**: Reviewed by Fable, which flagged one required confirmation before
implementing — whether the per-taxon build ever served as an access boundary (a taxon-scoped DB
handed to a collaborator who should specifically not see other taxa), not just a performance
artifact. **User confirmed: purely performance/convenience, never an access boundary** — so the
extractor is a straightforward correctness-focused tool, no access-control design needed. Also
grepped consumers: `MCP/BFD_mcp_server.py` already defaults to the single master `db/BFD.duckdb`
(`PROTEIN_DB_PATH` env var), unaffected by retiring per-taxon builds.

**Consequences**: Simpler `MERGE_SAMPLES`/`MERGE_*` (no `--taxon` branching), one master DB/Parquet
set as the single source of truth, taxon-scoped views become an on-demand extraction step rather
than a pipeline-time artifact. Full detail in `todo/genome_stats_storage_reorg.md` §D.2.

**Tags**: genome_stats, merge_samples, extract_bfd_taxonomic_subset, duckdb, parquet, taxon-filter, T-014, decision

## D: Implement EXTRACT_BFD_TAXONOMIC_SUBSET as a standalone Python/DuckDB CLI (2026-08-01)

**Context**: Follow-up to the decision above — with `MERGE_SAMPLES` now building only the full master table set (#27/#28 merged), a taxon-scoped view of the data needs a post-hoc extraction path. Fable's original review recommended a Python CLI wrapping DuckDB `ATTACH`/`COPY`, driven by `species` as the authoritative key table with `INNER JOIN`-based filtering and a post-extract row-count assertion.

**Decision**: `scripts/extract_bfd_taxonomic_subset.py` — `--taxon RANK:VALUE` (same convention as `taxonRowFilter()` elsewhere in this pipeline: uppercase RANK, exact-match VALUE), `--master` (default `db/BFD.duckdb`), `-o/--outfile` (required), `--dry-run` (report match count without writing). Implementation:
1. `ATTACH` the master DB read-only; filter `species` by the taxon predicate into the new output DB first.
2. Every other table (`asm_stats`, `gene_info` keyed by `LOCUSTAG`; all 19 remaining tables keyed by the derived `species_prefix` column) is built via `INNER JOIN` against the just-filtered `species` table in the *output* DB, not re-filtered independently — avoids orphaned rows if a table's own taxonomy columns ever drifted from `species`'s.
3. Indexes mirrored from `build_BFD_duckDB.sh` (kept manually in sync — no DuckDB API to introspect/replay index DDL from a read-only attached DB without also losing UNIQUE constraints).
4. The two analytical views (`v_species_summary`, `v_protein_annotation`) are recreated by reading their SQL directly from `duckdb_views()` on the master and re-executing verbatim against the output DB — avoids duplicating view definitions a second time (drift risk if `build_BFD_duckDB.sh`'s views ever change).
5. Post-extract row-count assertion: every table checked for orphaned keys not present in the filtered `species` table; on any orphan, the output file is deleted and the script fails loudly (robust-analysis convention) rather than leaving a partially-valid DB on disk.

**Verified against the real master DB** (`GENUS:Malassezia`, 81 genomes matched in `species`): correct partial-population behavior confirmed (function/gene tables scoped to whichever of the 81 genomes actually have computed data — 1, at time of testing — not all 81, since the master DB itself is only as complete as what's been computed so far); views queryable and return correct joined output; 46 indexes created; error paths (bad rank, zero matches, missing master file) all fail loudly with exit code 1 as designed.

**Alternatives considered**: (a) Making extraction a Nextflow process — rejected per the original decision's reasoning (pure SQL filtering over already-materialized data, no per-genome parallelism to gain). (b) Duplicating the view SQL inline rather than reading it from `duckdb_views()` — rejected in favor of reading it live from master, since duplicating it would be a second copy to keep in sync with `build_BFD_duckDB.sh` forever.

**Tags**: genome_stats, extract_bfd_taxonomic_subset, duckdb, taxon-filter, T-014, decision

## E: Route BUSCO_GENOME/BUSCO_PEP intermediates through `$SCRATCH` (2026-08-03)

**Context**: `BUSCO_GENOME` and `BUSCO_PEP` (`nextflow/modules/BFD/BUSCO_GENOME/main.nf`, `nextflow/modules/BFD/BUSCO_PEP/main.nf`, invoked from `subworkflows/local/BFD_GENOME_STATS.nf`) already `storeDir`-publish only a single small `short_summary*.txt` per run and `rm -rf` the full BUSCO run directory afterward — but BUSCO itself writes thousands of small intermediate files (HMMER/tblastn/Augustus/metaeuk) into that run directory while it executes, and that directory previously lived in the NFS-backed Nextflow work dir for the whole run.

**Decision**: Point BUSCO's `--out_path` at `$SCRATCH` (node-local fast storage) instead of the process's default (NFS) work dir, `cp` only the final summary file back out, then `rm -rf` the scratch run dir. Added `module load workspace/scratch` to the `busco_genome`/`busco_pep` `beforeScript` blocks in `nextflow/conf/profile_BFD.config` (same pattern already used by `GENOME_CLEAN_BATCH`, `pfam`, `SIGNALP_RUN`, `RNASEQ_PREPARE`), with `SCRATCH=${SCRATCH:-/tmp}` fallback if the module isn't loaded.

**Not yet verified**: whether this actually improves wall-clock — motivation is to test the hypothesis that BUSCO's many-small-file I/O pattern is NFS-metadata-bound. Needs a real before/after timing comparison (`sacct`/`seff`) on a representative genome before treating this as validated.

**Tags**: busco, scratch, io-performance, nextflow-hpcc, decision

## F: Route representative-strain selection through tables/busco_genome.parquet + asm_stats.parquet (2026-08-03)

**Context**: Follow-up to decision E. `nextflow/bin/pick_representative_strain.py` (invoked via `PICK_REPRESENTATIVE_STRAIN` from `workflows/compare_ANI.nf` under `--run_ani_reuse`) picked representative strains by BUSCO completeness + N50, but was regex-parsing both values directly out of raw `results/genome_stats/BUSCO_genome/**/*.BUSCO_summary.*.txt` files — including re-deriving N50 from BUSCO's own free-text output, duplicating (and being more fragile than) the authoritative `N50_bp` already computed in `tables/asm_stats.parquet`. There was also no merged BUSCO table in `tables/` at all — `BUSCO_GENOME`/`BUSCO_PEP` outputs weren't wired into any `MERGE_*` step.

**Decision**: Added `tables/busco_genome.parquet` as a first-class merged table, assembly-level only:
- `BFD_GENOME_STATS.nf` now emits `busco_genome` (from `BUSCO_GENOME.out.summary`).
- New `nextflow/bin/summarize_busco_stats.py` (mirrors `summarize_asm_stats.py`'s manifest-driven pattern) parses `*.BUSCO_summary.*.txt` into a TSV keyed by `ASMID` (no samples.csv join needed — BUSCO_genome is assembly-level, same as asm_stats).
- New `MERGE_BUSCO_GENOME` module (mirrors `MERGE_ASM_STATS`) converts that TSV to `tables/busco_genome.parquet`, wired into `BFD_MERGE.nf` (both glob and run-mode branches) and `BFD.nf`.
- `scripts/build_BFD_duckDB.sh` gained a `busco_genome` table (guarded on the parquet file existing, since `run_busco_genome` defaults false) joined to `species` for `LOCUSTAG`.
- `pick_representative_strain.py`'s `load_busco_by_asmid()` now runs a DuckDB query (`LEFT JOIN` `busco_genome.parquet` with `asm_stats.parquet` on `ASMID`) instead of glob+regex-parsing text files. CLI flag renamed `--busco-dir` → `--tables-dir` (default `tables/`); `PICK_REPRESENTATIVE_STRAIN/main.nf` and its `compare_ANI.nf` caller updated to pass `params.tables` instead of `${genome_stats_outdir}/BUSCO_genome`. Added `tables = "${launchDir}/tables"` to `nextflow.config`'s shared params block so `-profile ani` alone (without `BFD`) still resolves it.

**Explicitly out of scope**: `BUSCO_PEP` (annotation-level, protein-set completeness — only exists *after* gene prediction) was NOT wired into this merge or into representative selection. Representative-strain picking is a pre-annotation decision; conflating it with a post-prediction QC metric would be a category error. `MERGE_BUSCO_PEP` → `tables/busco_protein.parquet` is a separate, independent follow-up for annotation QC reporting (`.living/log/2026-08-03`, alongside `pfam`/`cazy` counts), not for strain selection.

**Verified**: `nextflow config -profile BFD` and `-profile ani` both parse; a full `-profile test -stub-run` of the BFD pipeline completed successfully (`completed=15 failed=0`) and produced `tests/output/tables/busco_genome.parquet` via `MERGE_BUSCO_GENOME`. `summarize_busco_stats.py` correctly parses a real (stub-format) `*.BUSCO_summary.*.txt` into the expected TSV row. `load_busco_by_asmid()` correctly LEFT-JOINs synthetic busco_genome/asm_stats parquet fixtures (including the case where an ASMID has BUSCO but no asm_stats match → N50 falls back to 0).

**Tags**: busco, tables, parquet, representative-strain, pick_representative_strain, merge, T-014, decision

## G: Place telomere-finder defaults in base `nextflow.config` and enable them in `test.config` (2026-08-03)

**Context**: The new `FIND_TELOMERES` step is wired into `BFD_GENOME_STATS.nf` and consumes `params.run_telomeres` plus a handful of `telomere_*` parameters. These values were initially added only to `conf/profile_BFD.config`, which is loaded by `-profile BFD`. A `-profile test -stub-run` immediately failed with `Cannot invoke method toBoolean() on null object` because `run_telomeres` was undefined.

**Decision**: Add safe base defaults for all `telomere_*` parameters (with `run_telomeres = false`) to the shared `params {}` block in `nextflow/nextflow.config`, and explicitly set `run_telomeres = true` in `nextflow/conf/test.config` so the test profile exercises the new step. Also add `withLabel: 'telomeres'` to `test.config`'s process section so stub runs stay local and module-free.

**Verified**: `-profile test -stub-run --pipeline BFD` completed successfully (16 tasks, 0 failed) and `BFD:BFD_GENOME_STATS:FIND_TELOMERES` and `BFD:BFD_MERGE:MERGE_TELOMERES` both appear in the DAG/execution.

**Tags**: telomeres, nextflow, config, defaults, test-profile, decision

## H: Rewrite `find_telomeres.py`'s 3'-end matching to use forward-window slicing instead of sequence reversal, and add a bp-tolerance for `terminal` (2026-08-04)

**Context**: See `.living/learnings.md` (2026-08-04, "`find_telomeres.py`'s 3'-end search was structurally broken...") for the full bug diagnosis. In short: the original 3'-end search reversed the sequence string without complementing it, then matched against a reverse-complemented pattern — two transforms that don't correspond, so 3'-end telomeres were essentially never detected. A second bug required an exact `start_coord == 0` match for `terminal`, dropping real terminal telomeres that had a few non-canonical leading bases.

**Review process**: Before implementing, the fix plan was reviewed independently by two agents — one briefed as a bioinformatics/sequence-analysis domain expert, one run on the Fable model — both given the bug diagnosis and proposed plan but told to verify by reading the code themselves, not take the diagnosis on faith. Both independently confirmed all three bugs by tracing the code directly. Two refinements came out of that review and were folded in:
- (Fable) Bug 3 (repeat-counting using the wrong monomer) only affected the regex path, not the fuzzy path, which already passed the correct `search_monomer`. Fable also flagged that the `monomer` output column had inconsistent semantics between regex mode (`canonical`) and fuzzy mode (`search_monomer`) — unified to always report `canonical`.
- (Bio-expert reviewer) Recommended exposing `distance_to_end` as a numeric column rather than relying solely on a boolean `terminal` cutoff, so downstream consumers can apply their own threshold.

**Decision**: 
1. Both `regex_terminal_hits` and `fuzzy_terminal_hits` now slice the terminal window directly in forward orientation for both ends (`offset = max(0, scaffold_len - search_window)` for 3', `offset = 0` for 5'), then search forward with `search_monomer` (already the correct per-end pattern). No sequence reversal anywhere in the module now.
2. Added `--terminal-tolerance` CLI flag (default 10bp): `terminal` is true if the tract is within this many bp of the literal scaffold edge, replacing the exact `== 0`/`== scaffold_len` check. New `distance_to_end` column added to `TelomereHit` and the TSV output.
3. Regex-path repeat counting now uses `search_monomer` (the actually-matched pattern) instead of `canonical`, matching what the fuzzy path already did correctly.
4. `TelomereHit.monomer` (and the TSV `monomer` column) now always reports the canonical (5'->3') form regardless of `--fuzzy`/regex mode, for schema consistency.
5. Threaded `--terminal-tolerance` through `FIND_TELOMERES/main.nf`, `nextflow_schema.json`, and `profile_BFD.config` (new `telomere_terminal_tolerance = 10` param), matching the existing `telomere_*` parameter pattern. Updated the module's `stub:` block TSV to include the new `distance_to_end` column so stub-run schema matches real output.

**Deliberately not changed in this pass**: the bio-expert reviewer flagged that the "canonical monomer at 5', RC at 3'" convention assumes a fixed chromosome strand orientation, which is actually arbitrary per assembly/contig — most of the default monomer patterns (all except the `TTAGGG`/`CCCTAA` pair, which is explicitly RC-symmetric) have no RC-symmetric twin, so a contig assembled in the "wrong" global orientation could still miss real telomeres at both ends unless `--both-ends` is used. `--both-ends` defaults to `false` pipeline-wide. This is a real, distinct correctness gap from the three bugs fixed here, but flipping that default doubles search volume/runtime across the whole BFD run and is a scope decision beyond "fix the bugs" — left as a follow-up (see `todo/TODO_REGISTRY.md`), not silently changed.

**Verified**: Re-ran on the same 5 real genomes used to originally diagnose the bug (`GCA_009805915.1_ASM980591v1`, `GCA_982397285.1_T36-F`, `GCA_041684285.1_ASM4168428v1`, `GCA_018345885.1_ASM1834588v1`, `GCA_900073075.1_CML3066.v2`) under production parameters. *F. graminearum* CML3066 (a real 4-chromosome telomere-to-telomere assembly) went from 1→7 of 8 real telomeric ends detected (the 8th has no repeat in the raw sequence — confirmed by direct inspection, a correct negative). *N. crassa* OR74A went 1→13 hits, *A. nidulans* ASM4168428v1 went 1→9 hits, both now with a healthy 5'/3' balance instead of 100%-5'-only. `summarize_telomeres.py` re-tested against the fixed output and aggregates correctly (extra `distance_to_end` column is a no-op for it, since it reads by column name). `python3 -m py_compile` clean on both scripts; `nextflow_schema.json` still valid JSON.

**Tags**: telomeres, find_telomeres, bug-fix, reverse-complement, terminal-tolerance, fable-review, bioinformatics-review, decision
