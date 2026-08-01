# Learnings

Append-only log of gotchas, surprises, and insights.

**Entry template:** copy from `skills/core/templates/learning-entry.md` (includes Category, What happened, Why it matters, Resolution, Tags fields). The `**Tags**:` line is consumed by `generate_index.py --summary-heuristic` to build the cluster summary in INDEX.md — use them.

### [2026-06-25] MERGE_* steps: pass a manifest file, not thousands of staged inputs; bake mtime+size in for staleness

**Category**: insight

**What happened**: The BFD.nf MERGE_* processes (aa_freq, codon_freq, intergenic, gene_stats, asm_stats) declared `input: path 'inputs/*'` and were handed the full `.collect()` of per-genome stat files. At BFD scale (~8k+ genomes × up to 7 gene_stats categories) that means Nextflow stages tens of thousands of symlinks per merge and builds a huge stage-in command — slow and risking "Argument list too long". Refactored to a single **manifest file**: a `toManifest(ch, name)` helper flattens the file channel and `collectFile`s one TAB-delimited record per line — `<abs_path>\t<mtime_ms>\t<size_bytes>` — and each MERGE process takes `path manifest`, reading the path from field 1 (`while IFS=$'\t' read -r f _mtime _size`). Files are read directly from their stable storeDir absolute paths (shared FS), not staged.

**Why it matters**: (1) Eliminates the massive stage-in / symlink-per-file overhead and the ARG_MAX risk for any kingdom-wide merge. (2) **Staleness**: with `path 'inputs/*'` Nextflow hashed file contents, so a regenerated input forced a re-merge under `-resume`. A bare path-list manifest would *lose* that — same paths → same content hash → stale merge silently kept. Baking mtime+size into the manifest *content* restores it: any input regenerated (newer mtime / different size) changes the manifest bytes → merge task input-hash changes → Nextflow re-runs on `-resume`; unchanged → byte-identical manifest → cache hit/skip. This is the staleness check the manifest enables, done via content hashing rather than a separate `find` mtime scan (which is fragile because the freshly-written manifest is always newest, so you'd have to diff against a *persisted* prior manifest and still force a cache miss manually).

**Resolution**: `toManifest()` helper + `path manifest` inputs in all five MERGE_* processes (`nextflow/BFD.nf`); both `use_glob` and non-glob workflow branches wrap their MERGE args in `toManifest(...)`. `scripts/summarize_asm_stats.py` rewritten with argparse to accept `--manifest` (reads path from field 1) and gzip output — it previously had no argparse at all yet the process called it with `--reportdir/--samples/-o`, a latent always-fail. Verified: merge bash logic over a tab manifest (header-once + type-suffix filtering for gene_stats), the python script over a manifest, and `nextflow ... -profile BFD,test -preview` → `[SUCCESS]`.

**Tags**: nextflow, BFD.nf, merge, manifest, staging, arg-max, staleness, resume, caching, collectFile, asm-stats

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped — the manifest pattern is the structural fix. Any future many-file merge in this repo should reuse `toManifest()` rather than `path 'inputs/*'` + `.collect()`.

### [2026-06-25] Nextflow `combine` spreads list-valued left items — use a scalar barrier

**Category**: gotcha

**What happened**: BFD.nf crashed at runtime with `MissingMethodException: ...closure217.call() is applicable for argument types: (LinkedList) values: [[.../codon_freq/Aaosphaeria_arxii...csv.gz, ...]]`, followed by `CALC_ASM_STATS` / `MERGE_*` "Cannot obtain the semaphore" + `InterruptedException` (just the teardown cascade). Root cause was in the non-glob MERGE branch (`merge_all=false` or `--taxon` active): the sync/barrier channel was `BATCH_CODON_FREQ.out.csv.flatten().collect().ifEmpty([])`, which emits a **List** value. Nextflow's `combine` operator *spreads* a list-valued left item, so `codon_sync.combine(codon_paths)` did not produce the expected 2-tuple. When all BATCH outputs were already cached, `codon_sync` resolved to the empty list `[]`, spread to zero elements, and `combine` emitted just the single `codon_paths` list `[[files...]]`. The downstream 2-arg closure `.map { _s, files -> ... }` couldn't bind a single LinkedList → fatal. The same latent bug existed in the `aa_freq` site.

**Why it matters**: `combine` flattening tuples is documented but easy to forget when the left operand is meant only as a completion gate. Any "wait for upstream, then attach a collected file list" pattern that puts a `.collect()`/`.ifEmpty([])` channel on the *left* of `combine` and then destructures with a multi-arg map closure will silently misbind — and only fails in the all-cached path, so it can pass early tests and break later.

**Resolution (CORRECTED 2026-06-25 same day — first attempt was incomplete)**: `combine` spreads list-valued items on **BOTH** sides, not just the left. The first fix made only the barrier scalar (`.collect().map { true }.ifEmpty(true)`) but left the file channel as a `.collect()` list on the right; runtime then produced one FLAT emission `[true, f1, f2, …, fN]` (verified: `MissingMethodException ...closure216 ... (LinkedList) values: [[true, /…/aa_freq/Aaosphaeria…csv.gz, …]]`) which `{ _s, files -> }` still can't bind. The correct idiom is `scalarBarrier.combine(perItemChannel)` — the file channel must emit **one file per item** (do NOT `.collect()` it before combine) — yielding clean `[true, f]` 2-tuples for `.map { _s, f -> f }.filter { it.exists() }`. Confirmed with a standalone probe: `Channel.of(true-value).combine(Channel.of('f1','f2','f3'))` → `[true,f1] [true,f2] [true,f3]`. Fixed both `aa_freq` and `codon_freq` non-glob blocks in nextflow/BFD.nf (collected list downstream via the [[merge-manifest]] `toManifest()` instead). Note: `-profile test` alone is insufficient for `-preview` — `genome_stats_outdir`/`freq_batch_size` live in `conf/profile_BFD.config`, so lint/preview must use `-profile BFD,test`. Also note: `-preview` does NOT execute operators, so it does not catch this class of runtime misbind — a real (or stub) run is required to validate `combine`/`map` destructuring.

**Tags**: nextflow, combine-operator, channel-semantics, barrier-channel, BFD.nf, caching, preview-lint

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped for both sites. Repo convention: "to gate items on a completion barrier, use `scalarBarrier.combine(perItemChannel)` — never `.combine()` a `.collect()`ed list on either side; collect AFTER the combine/map/filter." Plus `-profile BFD,test` for lint/preview. A stub-run smoke test would catch operator-level misbinds that `-preview` misses.

### [2026-06-19] Any BFD genome-size analysis must control two confounds + needs the All_Taxa merge

**Category**: insight

**What happened**: A 9-persona ideation session over `tables/asm_stats.tsv.gz` (22,412 assemblies / 8,061 species) surfaced the same two confounds independently from nearly every disciplinary lens, plus a shared data prerequisite: (1) **assembly-quality** — `total_length_bp` and `masked_pct` both depend on contiguity/fragmentation (`contig_count`/`N50`/BUSCO), so raw size/repeat comparisons are biased; (2) **taxonomic non-independence** — many strains per species, phylogenetically clustered, so per-assembly regression pseudoreplicates. Also: the kingdom-wide `All_Taxa` compositional merge (aa_freq, codon_freq, intergenic, pfam/cazy/merops, idp…) is **not yet materialized**, so most composition-vs-size analyses need one BFD.nf merge run first. And `asm_stats` genome-size tails (~0 Mb and ~1711 Mb) are likely quality artifacts, not biology.

**Why it matters**: These apply to *every* future genome-size/composition analysis in this repo. Skipping the quality control inflates apparent size↔repeat coupling; ignoring phylogenetic structure overstates significance; analyzing the size tails without filtering pollutes results.

**Resolution**: Documented in `.living/findings/genome-size-architecture.md` and baked into the ideation `CONTEXT.md` so future analyses inherit the constraints. Recommended framework deliverable (idea #12) is explicitly a *quality-residualized* embedding.

**Tags**: genome-size, comparative-genomics, confounds, assembly-quality, phylogenetic-non-independence, asm_stats

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A shared loader for `asm_stats` that (a) filters/flag the implausible size tails, (b) attaches a quality score from contig_count/N50/BUSCO, and (c) attaches derived taxonomic ranks — used by all genome-size analyses — would make this structural.

### [2026-06-18] scilintr needs a python3.12 venv here; flags "old" inside "threshold"

**Category**: tip

**What happened**: `scilintr` (mycelium analyze step) requires Python ≥3.11 and is not installable into the cluster's pandas-bearing conda 3.9; `/usr/bin/python3.12` has no pip. Created a venv: `python3.12 -m venv /tmp/scilintr_venv && /tmp/scilintr_venv/bin/pip install scilintr`, then lint via `/tmp/scilintr_venv/bin/scilintr <path>`. Separately, its `implicit-file-selection` rule false-positived on the output filename `threshold_sensitivity.csv` — it matched the substring "old" inside "thresh**old**". Renamed the file to `sensitivity_sweep.csv` to avoid a vacuous waiver.

**Why it matters**: `/tmp` is ephemeral — the venv won't persist across sessions, so the install must be redone (or moved somewhere durable). And any analysis output/variable containing the substrings "old"/"latest"/"backup" will trip `implicit-file-selection` regardless of intent; prefer names that avoid them over waivers.

**Resolution**: Documented venv recipe; renamed offending output. Consider a durable venv location (e.g. under the repo or `~/.local`) if scilintr is used regularly.

**Tags**: scilintr, mycelium, python-version, linting, tooling, false-positive

**mitigation_type**: ambient-awareness

**structural_mitigation_candidate**: A committed `skillpacks/` or `tools/` bootstrap script that creates a durable scilintr venv and is referenced by analysis run scripts would make this structural. Until then it is awareness-only.

### [2026-06-18] TrinityGG_failed.tsv LOW_COUNT rows double-list successes, not failures

**Category**: gotcha

**What happened**: A strict success∩failure disjointness assertion in the QC pipeline fired: 5 species appeared in BOTH `misc/TrinityGG_summary.tsv` (successes) and `misc/TrinityGG_failed.tsv`. Investigation showed these are exactly the 5 `REASON=LOW_COUNT` rows — the assembly FASTA *exists* (hence it's a success) but is also flagged in the failed report. Their `NUM_TRANSCRIPTS` match exactly across both files. The other failure reasons (MISSING=874, EMPTY=65) are true no-assembly failures and are disjoint from successes.

**Why it matters**: Treating the failed report as "all failures" double-counts the 5 LOW_COUNT species, inflating the failure denominator (1183 vs the correct 1178 attempted) and the failure rate. Any downstream "how many assemblies failed" stat must split LOW_COUNT (a quality flag on an existing FASTA) from MISSING/EMPTY (no FASTA).

**Resolution**: `reconcile_landscape()` now partitions the failed report into no-assembly (MISSING/EMPTY) vs LOW_COUNT, asserts no-assembly is disjoint from successes, asserts LOW_COUNT ⊆ successes, and cross-checks that LOW_COUNT NUM_TRANSCRIPTS agree between the two reports. LOW_COUNT species are reconciled into the success set and handled by the FAIL/BORDERLINE tiers.

**Tags**: rnaseq, trinity, qc, data-reconciliation, double-counting, denominator

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped — the disjointness + subset + cross-report-agreement assertions in `analysis/rnaseq_trinity_qc/scripts/run.py::reconcile_landscape` now fail loudly if the success/failure sets overlap in any way other than LOW_COUNT.

### [2026-06-18] Mycelium scripts need python3.12, not the default python3

**Category**: gotcha

**What happened**: On the UCR HPCC the default `python3` is 3.9 (system miniconda). Mycelium scripts (`init_repo.py`, `generate_index.py`, `recall_lessons.py`, `install_convention.py`) fail on import under 3.9 with `ImportError: cannot import name 'UTC' from 'datetime'` and `TypeError: unsupported operand type(s) for |: 'type' and 'NoneType'` (PEP 604). `python3.12` (`/usr/bin/python3.12`) runs them fine. Also: `init_repo.py` could not auto-locate the `network/conventions/` dir, so the three core packs had to be installed manually with `install_convention.py --network-dir`.

**Why it matters**: Anyone re-running mycelium tooling here with bare `python3` will hit confusing import errors. The SessionStart `mycelium-health.sh` hook invokes `generate_index.py` via plain `python3`, so `.living/INDEX.md` will silently fail to auto-regenerate — the index must be rebuilt manually with `python3.12` when learnings/decisions change.

**Resolution**: Always call mycelium scripts with `python3.12`. Documented in `ENVIRONMENTS_INSTALLATIONS.md`. Core packs (robust-analysis, report-generator, idea-generator) installed via explicit `--network-dir $MYC/network/conventions`.

**Tags**: mycelium, python-version, hpcc, environment, hooks, tooling

**mitigation_type**: convention

**structural_mitigation_candidate**: A wrapper that shims `python3` → `python3.12` for mycelium invocations, or a `requires-python` preflight check in the hook scripts that warns when the resolved `python3` is <3.11. Until shipped, mitigated by the ENVIRONMENTS_INSTALLATIONS.md note (convention).

---

## L: Upstream ncbi_accessions_taxonomy.csv is stale/duplicated (2026-06-19)

**Observation**: `../../1KFG/2026/NCBI_fungi/ncbi_accessions_taxonomy.csv` covers only 8,061 of 22,412 accessions and every row is duplicated ~2.8× (22,412 lines / 8,061 unique). The committed `samples.csv` predates this and was built from a complete file. The old `create_samples_file.py` would also silently corrupt rows here: its merge does `dict[asm].extend(row[4:])`, so any duplicated assembly multi-appends columns.

**Impact**: A fresh run today leaves ~14k assemblies with empty lineage + `fungi` BUSCO default and species falling back to the strain-laden verbatim name.

**Fix / mitigation**: Regenerate upstream (`cd ../../1KFG/2026/NCBI_fungi && pixi run make lib/ncbi_accessions_taxonomy.csv`). New `create_samples_file.py` dedups taxonomy by first occurrence, merges by column name (never extend), and logs coverage so the staleness is visible instead of silent.

**Tags**: ncbi, taxonomy, taxonkit, stale-data, data-validation, samples-csv

---

## L: comparative_genomics.nf STAGE_FILES expected wrong input filenames (2026-06-22)

**Category**: gotcha

**What happened**: `comparative_genomics.nf` `STAGE_FILES` symlinked inputs by `{LOCUSTAG}.faa` and `{LOCUSTAG}.cds.fa`, but the actual files in `input/pep/` and `input/cds/` are named by sanitized `{species}_{strain}` with different extensions: `{tag}.proteins.fa` and `{tag}.cds-transcripts.fa` (e.g. `Aaosphaeria_arxii_CBS_175.79.proteins.fa`). The `[ -f "$src" ]` guard meant the loop silently matched nothing, so MMseqs/DIAMOND/OrthoFinder would have clustered an empty FASTA with no error. The workflow had never been run (no prior `comparative/` output), so the mismatch was undetected.

**Why it matters**: Any first run of the comparative clustering workflow (mmseqs/mcl/orthofinder) would have produced empty or failed results that look superficially successful. The silent-skip pattern (`[ -f ] && ... && ln`) hides missing-input bugs.

**Resolution**: `PREPARE_COMPARATIVE` now derives a `BASENAME` column via `SampleUtils.makeSampleTag(SPECIES, STRAIN)` — the *same* canonicaliser `funannotate.nf` uses to write those files (`nextflow/lib/SampleUtils.groovy`, auto-loaded) — and carries it as a 3rd manifest column (`LOCUSTAG,GROUP,BASENAME`). `STAGE_FILES` resolves `{BASENAME}.proteins.fa`/`{BASENAME}.cds-transcripts.fa`, symlinks to `{LOCUSTAG}.faa`/`{LOCUSTAG}.cds.fa` (preserves stable IDs + downstream `cat *.faa`), WARNs on each missing input, and exits non-zero if zero files stage. Verified with `-stub-run` on `ORDER:Capnodiales` (3/6 staged, 0 broken; the other 3 genuinely lack proteins in `input/pep`). NB: filename reconstruction must always reuse `makeSampleTag`, never reconstruct ad hoc — strain cleaning (semicolon/asterisk/colon/quote handling) is non-trivial.

**Tags**: nextflow, comparative-genomics, mmseqs, staging, samplesutils, makesampletag, silent-failure, filename-convention

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Any new process that resolves per-genome input files by name should call `SampleUtils.makeSampleTag` rather than assume LOCUSTAG-named files; consider a shared staging helper.

---

## L: funannotate "Not enough gene models N to train Augustus (30 required)" is the ground-truth too-small signal; asm_stats table is stale (2026-06-24)

**Category**: gotcha / data-validation

**What happened**: Screened all 8,082 persisted `funannotate-predict.log` files for genomes that fail prediction because the assembly is too small/fragmented to yield Augustus training models. The authoritative, greppable failure line funannotate emits is `Not enough gene models N to train Augustus (30 required), exiting`. Only **9 genomes** hit it (plus 12 other non-completions); **8,054 completed**. All 9 are tiny/fragmented (127 kb–15.7 Mb); rusts (Uromyces/Puccinia, should be 100 Mb+) sitting at 2–7 Mb are grossly partial assemblies.

Two non-obvious gotchas surfaced:
1. **`results/genome_stats/asm_stats/` (4,094 `*.stats.txt`) is STALE** — built on an earlier genome set; **none** of the 21 failing genomes had an entry. So the asm_stats pre-screen cannot validate/catch current failures until `CALC_ASM_STATS` is re-run over the current `samples.csv`. Compute metrics on the fly (`seqkit stats`/AAFTF) for genomes missing from the table.
2. **Join key**: predict-log line 1 carries `--name <LOCUSTAG>`; LOCUSTAG is **unique** across all 22,929 `samples.csv` rows (0 dups) → clean LOCUSTAG→ASMID→PHYLUM join. Do NOT parse ASMID from the `-i` path: gzipped genomes are inflated to `genome_input.fa`, erasing the ASMID.

**Why it matters**: A flat genome-size cutoff is wrong here — *Malassezia* (~7–9 Mb, 8 contigs, N50 >1 Mb) is complete and annotates fine, while a fragmented 40 Mb rust completes too. The only assembly-stat rule that cleanly separates real failures is **small AND fragmented** (assembled bp < ~16 Mb AND (N50 < ~10–20 kb OR contigs > ~1000)). This catches 9/9 too-few-models failures with zero false-skips of complete small genomes. Verified safe against the small-genome taxa the user flagged: Ashbya/Eremothecium yeasts (8.9–9.7 Mb, N50 ~1.5 Mb), complete Microsporidia, Rozella (11.3 Mb) and Paramicrosporidium (7.2 Mb, N50 70 kb) all pass (frag=0). The *opposite* failure class — large genomes (0.5–1.3 Gb rusts; good 18–34 Mb assemblies) that show `INCOMPLETE_unknown` — is a resources/walltime problem, NOT detectable by a size filter.

**Resolution**: `analysis/funannotate_model_failures/parse_predict_failures.py` does the scan/join/crosstab; in-pipeline guard added to `funannotate.nf` (see decision). 7 of the 9 added to `suppress.txt` (2 already there); 3 of the original 21 were already manually suppressed — the guard automates that curation.

**Tags**: funannotate, augustus, training-models, asm-stats, stale-data, suppress-txt, locustag, basidiomycota, rust-genomes, malassezia, data-validation

**mitigation_type**: code-fix + analysis

**structural_mitigation_candidate**: Re-run `CALC_ASM_STATS` whenever `samples.csv` gains genomes, so the asm_stats table never lags the annotated set; the merged `asm_stats.tsv.gz` is the natural input to an automated too-small pre-screen.

### [2026-07-16] SRA_FETCH_SE must dump single-end runs with --split-3, not --split-files

**Category**: gotcha

**What happened**: `SRA_FETCH_SE` (nextflow/funannotate.nf) downloaded every SE accession with `parallel-fastq-dump ... --split-files`. `--split-files` is a *paired* operation — it separates R1/R2 mates into `_1`/`_2` files — and it fails on genuine single-end runs (SINGLE layout), which have no second mate to split out. Switched the SE process to `--split-3`, which routes reads by what each spot actually contains: genuine SINGLE → all reads in `ACC.fastq.gz` (no `_1`/`_2`); SE_trinity with 2 reads/spot → `ACC_1`/`ACC_2` (take `_1`); SE_trinity with 1 read/spot → unpaired reads in `ACC.fastq.gz`. The existing SE_FILE detection already prefers `_1` then falls back to `ACC.fastq.gz`, so it handles all three cases unchanged. The two `--split-files` calls in the PE `SRA_FETCH` were left as-is: that process genuinely wants pairs and already flags a lone `_1` as an SE candidate.

**Why it matters**: Any species routed to SE fetch (blacklist action `SE_trinity`, or query CSV `layout==SINGLE`) with true single-end data was failing the dump outright, leaving a zero-byte `_norm_SE.fastq.gz`. Separately, that missing/empty SE output is *also* what breaks `storeDir` completeness and re-triggers fetch — see the incomplete-store issue for Pseudocercospora_eumusae (reads built before the `_norm_SE` output was added on 2026-06-12, so the store looks incomplete and re-runs). Routing (5-column query CSVs with no `layout`, not blacklisted `SE_trinity`) is a distinct problem: such SE species still fall into the PE path and fail there.

**Resolution**: `--split-files` → `--split-3` in `SRA_FETCH_SE` plus comment updates (nextflow/funannotate.nf). Verified `parallel-fastq-dump --help` confirms extra args pass through to fastq-dump, and `nextflow config -profile standard` parses the edited module without error. Not yet exercised on a live SE accession end-to-end.

**Tags**: nextflow, funannotate, SRA_FETCH_SE, parallel-fastq-dump, split-3, single-end, rnaseq, storeDir, sra-layout

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Backfill zero-byte `_norm_SE.fastq.gz` stubs for the 417 species that have real R1/R2 but no SE file, so `storeDir` stops re-triggering PE fetch on pre-2026-06-12 reads; and refresh the 3013 five-column `sra_query` CSVs to the 6-column (with `layout`) format so genuine SINGLE species route to SE fetch instead of failing in the PE path.

### [2026-07-16] Normalized 3101 cached sra_query CSVs to the 6-column (layout) schema in place

**Category**: insight

**What happened**: `rnaseq_reads/sra_query/*.sra_query.csv` existed in three schema generations from successive SRA_QUERY/SRA_QUERY_BATCH versions: 4-column (`...,spots`, oldest, 88 files, all header-only), 5-column (`...,platform`, 3013 files), and 6-column current (`...,platform,layout`, 2103 files). Only the 6-column form carries `layout`, and SRA_FETCH_SE routes an accession to single-end handling only when `layout==SINGLE` (a missing column defaults to PAIRED), so pre-6-column files could never route to SE fetch. Wrote `scripts/one-off/refresh_sra_query_layout.py` (dry-run default, `--apply`, atomic per-file `os.replace`) to stamp the 6-column schema in place. The stamped layout is provably PAIRED: every pre-6-column query was esearch-filtered on `PAIRED[Layout]` (verified at git f7db039~1), so every 4/5-column data row is PAIRED by construction; the script refuses to guess (skips any 4-column *data* row, of which there were none). Result: all 5204 files now 6-column, 4350 data rows, zero ragged rows.

**Why it matters**: (1) Uniform schema — COLLECT_SRA_QUERY's merged manifest and any strict 6-column parser now see consistent input. (2) It is a **pure format normalization, not a behavior fix**: because the consumer already defaults missing `layout` to PAIRED, stamping PAIRED changes no routing. In particular it does NOT fix mislabeled-SE species like Pseudocercospora_eumusae, whose accessions SRA genuinely reports as PAIRED (they now read `...,PAIRED`) yet are single-end — those still route to PE `SRA_FETCH` and fail, and need a `SE_trinity` blacklist entry instead. (3) Discovering *new* single-end data for the header-only files (1620 old 5-col + 2004 6-col + 88 4-col = empty queries) still requires deleting those cached CSVs and re-running so SRA_QUERY re-queries with the SINGLE fallback — the conversion alone adds no accessions.

**Resolution**: `scripts/one-off/refresh_sra_query_layout.py --apply` run once; 3101 files converted (88 4→6, 3013 5→6), 2103 already current. Pairs with the earlier `--split-3` fix so genuine SINGLE accessions, once present in a 6-column CSV, both route to and dump correctly in SRA_FETCH_SE.

**Tags**: nextflow, funannotate, sra_query, schema, layout, single-end, rnaseq, one-off, data-normalization

**mitigation_type**: data-fix

**structural_mitigation_candidate**: Next, (a) re-query the ~3700 header-only species (delete their cached CSVs, re-run) to pick up SINGLE-layout RNA-seq the old PAIRED-only queries missed; (b) promote detected SE candidates (rnaseq_se_candidates.csv) to `rnaseq_blacklist.csv` as `SE_trinity` so mislabeled-PAIRED species like Pseudocercospora_eumusae route to SE fetch; (c) backfill zero-byte `_norm_SE.fastq.gz` stubs so storeDir stops re-triggering PE fetch on pre-2026-06-12 reads.

### [2026-07-16] Backfilled 417 SE stubs + quarantined 3624 empty sra_query CSVs to trigger SE re-query

**Category**: insight

**What happened**: Two coupled data-fix operations, each a one-off script (dry-run default, `--apply`):
(c) `scripts/one-off/backfill_se_stubs.py` created zero-byte `_norm_SE.fastq.gz` for the 417 species that had real (non-empty) paired R1+R2 but no SE file (reads built before the `_norm_SE` output existed, 2026-06-12). This completes their `storeDir` so `SRA_FETCH` stops needlessly re-downloading. Critically it stubs ONLY genuine-PE species: the 228 species with empty R1+R2+missing-SE are left untouched so they stay eligible to pick up single-end data (a stub would make storeDir consider them done and block SRA_FETCH_SE forever). Verified: 0 real-PE species remain with missing SE; 228 empty ones intentionally left.
(a) `scripts/one-off/requery_empty_sra_query.py` quarantined (moved, not deleted, to `rnaseq_reads/sra_query_requery_quarantine_20260716/`) the 3624 header-only sra_query CSVs so both the SRA_QUERY storeDir check and the SRA_QUERY_BATCH `[ -s ]` reuse check miss → those species re-query on the next run. Also flipped `enable_single_end = true` in conf/profile_funannotate.config so the re-query actually runs the SINGLE[Layout] fallback (was false; a bare re-query would only redo the PAIRED search). Species list saved. 1580 data-bearing CSVs remain in place.

**Why it matters**: (c) stops the expensive (32cpu/96GB/2h) PE re-fetch churn for 417 species that already have good reads, while preserving SE-discovery eligibility for the 228 empty ones. (a) is what actually gives ~3624 species a first-ever single-end query — the old PAIRED-only queries never checked for SE data. The two are sequenced deliberately: the SE stubs must NOT cover species that (a) might newly populate with SINGLE data. After the pipeline re-runs, newly-found SINGLE accessions route to SRA_FETCH_SE (now `--split-3`, so single-end dumps correctly) and produce real `_norm_SE.fastq.gz` for the previously-empty species.

**Resolution**: both scripts applied once. `enable_single_end=true` committed to profile. Next action is the user's pipeline re-run, which will (1) re-query the 3624 quarantined species with PAIRED+SINGLE, (2) SE-fetch any SINGLE hits. Mislabeled-PAIRED species (e.g. Pseudocercospora_eumusae) still need `SE_trinity` blacklist entries — unaffected by this.

**Tags**: nextflow, funannotate, sra_query, storeDir, single-end, rnaseq, backfill, re-query, enable_single_end, one-off, data-fix

**mitigation_type**: data-fix + config

**structural_mitigation_candidate**: Promote detected SE candidates (rnaseq_se_candidates.csv) to rnaseq_blacklist.csv as `SE_trinity` after the re-query run so mislabeled-PAIRED species route to SE fetch. Longer term, have SRA_FETCH write the `_norm_SE` stub unconditionally (already does) and consider a preflight that flags incomplete stores rather than silently re-running.

### [2026-07-16] Validated --split-3 SE fix: 2-sample SRA_FETCH_SE regression test passed

**Category**: insight

**What happened**: Ran a deterministic 2-sample test of the SRA_FETCH_SE `--split-3` fix on branch fix/rnaseq-single-end-handling. Species: Sporisorium_graminicola + Malassezia_globosa — both had empty reads (needed regeneration), SE_trinity blacklist entries, and cached clean genomes. Setup (do_annotation/test_se/): 2-row samples.csv; pre-staged sra_query CSVs listing each species' real SE_trinity accessions (the natural PAIRED esearch filters them out on readlength/spots, so they'd otherwise route to no_data); deleted the 0-byte R1/R2/SE stubs so storeDir was incomplete (per the funannotate.nf:850 note, re-routing PE→SE requires deleting the norm_* files). Ran with -params-file params.json (skip_sra_query, stop_after_sra_fetch, enable_single_end all true — note: nf-schema rejects `--flag true` as a string, booleans MUST come via -params-file or config, not bare CLI). Result: `[SUCCESS] completed=3 failed=0`, both species routed to SRA_FETCH_SE and produced real non-empty SE reads (Sporisorium 2,370,539 reads/71MB; Malassezia 226,096 reads/14MB) with R1/R2 correct 0-byte stubs, zero retries. Under the old `--split-files` this single-end data would have failed.

**Why it matters**: End-to-end confirmation the SE path works: SINGLE download via --split-3 → header-fix → bbnorm → fastp → complete storeDir. Two operational gotchas captured: (1) boolean params via -params-file only (nf-schema strict typing); (2) `cleanup = true` in profile_funannotate.config deletes work dirs on success — so a missing task work dir is evidence of success, and failed tasks are the ones whose dirs persist (which is why the Pseudocercospora failure dirs were still around to inspect).

**Resolution**: Fix validated on real single-end SRA data. Test scaffold kept at do_annotation/test_se/ (samples, params.json, run_test_se.sh). The two species' staged CSVs + regenerated SE reads left in place (they genuinely need SE handling). Branch ready for merge.

**Tags**: nextflow, funannotate, SRA_FETCH_SE, split-3, single-end, rnaseq, test, nf-schema, params-file, cleanup, validation

**mitigation_type**: validation

**structural_mitigation_candidate**: For future SE routing of mislabeled-PAIRED species whose accessions the readlength/spots filter drops, the SE_trinity accessions must be seeded into the sra_query CSV (the natural query won't surface them) — consider a query variant that unions blacklist SE_trinity accessions into the per-species CSV.

### [2026-07-19] COMBINE_ANI_TABLE crashed: unquoted 'group' column is a SQLite reserved keyword

**Category**: bugfix

**What happened**: `run_ANI.sh` (do_ANI/) failed at the `COMBINE_ANI_TABLE` step after all `SKANI_SKETCH`/`SKANI_COMPARE` work completed (32/35 tasks succeeded, retried twice, then gave up). `.nextflow.log` showed `sqlite3.OperationalError: near "group": syntax error` from `nextflow/bin/combine_ani_table.py:135`. The script builds `CREATE TABLE ani_pairs (...)` and `CREATE INDEX ...` by joining column names unquoted, and one output column is literally named `group` — a SQLite reserved keyword, so the DDL wouldn't parse.

**Why it matters**: Column names chosen for readability/domain-fit (`group`, `order`, `table`, `select`, etc. are all SQL reserved words) can silently break DDL even though `INSERT ... VALUES (?, ?, ...)` positional inserts elsewhere in the same script work fine — the bug only surfaces at table/index creation, not at insert time, which is why it wasn't caught until this stage of a `-resume` run.

**Resolution**: Double-quoted all column identifiers in the `CREATE TABLE`/`CREATE INDEX` statements in `nextflow/bin/combine_ani_table.py` (kept as defensive practice), and renamed the `group` column itself to `taxon_group` throughout the script (CSV header, SQLite DDL, index list) to permanently avoid colliding with the reserved word rather than relying on quoting at every future query site. Verified end-to-end with a synthetic 1-pair/2-genome fixture (CSV + SQLite output both correct, index created as `idx_ani_pairs_taxon_group`). No other file in the repo referenced the old `group` column name. Since `run_ANI.sh` uses `-resume`, re-running `sbatch run_ANI.sh` reuses cached SKANI_SKETCH/SKANI_COMPARE results and only reruns COMBINE_ANI_TABLE + REPORT_ANI.

**Tags**: nextflow, sqlite, do_ANI, combine_ani_table, reserved-keyword, bugfix, column-rename

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Any future script that builds SQL DDL from a dynamic column list should double-quote identifiers by default rather than only when a collision is discovered; also prefer domain-specific column names (`taxon_group`, not `group`) over generic ones that happen to be SQL keywords.

### [2026-07-19] all_pairs.csv/ani.db had blank genus/species/strain for every row — names-file naming mismatch

**Category**: bugfix

**What happened**: After fixing the `taxon_group` reserved-keyword bug and re-running `run_ANI.sh` (job 26475151, COMBINE_ANI_TABLE served from a valid `-resume` cache hit), spot-checking `results/ANI/skani/GENUS/ani.db` showed every one of 263,560 rows had blank `query_asmid`/`query_genus`/`query_species`/`query_strain`/`ref_*` fields — only `taxon_group`/`ani`/filenames were populated. Root cause: `compare_ANI.nf` has two independent conventions for the per-group genome-names lookup file — `REPORT_ANI` copies it to `{group}_genome_names.tsv` (matches `combine_ani_table.py`'s `NAMES_SUFFIX`), but `COMBINE_ANI_TABLE`'s actual input channel (`names_map`, compare_ANI.nf:596/708) is fed the raw `names_{group_name}.tsv` file created directly via `file("${workflow.workDir}/names_${group_name}.tsv")` — a different prefix-based naming convention the Python script's `group_from_names_path()` and its names-dir glob (`*_genome_names.tsv`) never matched, so every names file was silently skipped and every genome fell back to `blank_names()`. This is exactly the class of silent-failure bug the mycelium robust-analysis convention warns about: the join failure produced no error, just empty columns.

**Why it matters**: The `.ani.tsv` glob (`*.ani.tsv`) still matched correctly, so `COMBINE_ANI_TABLE` "succeeded" (exit 0) and passed validation for months/runs — only a manual spot-check of populated columns caught it, not the pipeline itself.

**Resolution**: `nextflow/bin/combine_ani_table.py` `group_from_names_path()` now recognizes both `{group}_genome_names.tsv` and `names_{group}.tsv`; the names-dir glob widened from `*_genome_names.tsv` to `*.tsv` (filtered through `group_from_names_path`) so it actually finds the `names_*.tsv` files COMBINE_ANI_TABLE receives. Verified against the real staged `names_Aspergillus.tsv` from work dir `35/7293b5dc3f191a0de015a7f0268f25`: 0/263560 blank rows after the fix (all 263560 were blank before). Re-running `run_ANI.sh` will re-execute COMBINE_ANI_TABLE (script changed → cache invalidated) and republish corrected `all_pairs.csv`/`ani.db`.

**Tags**: nextflow, do_ANI, combine_ani_table, names_map, silent-failure, bugfix, data-integrity

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Add an assertion/sanity check to `combine_ani_table.py` (e.g. warn or fail if >X% of rows have blank names fields) so a names-file join failure surfaces immediately instead of requiring a manual spot-check. Also consider unifying the two names-file naming conventions in compare_ANI.nf into one to eliminate the divergence at its source.

### [2026-07-23] nextflow_memory_profile analysis was scaffolded but never actually run; fixed two bugs and reran

**Category**: bugfix

**What happened**: `analysis/nextflow_memory_profile/` was created 2026-07-17 with a doc, `run.sh`, and `scripts/profile_nextflow_memory.py`, but `outputs/` was empty and the doc's Key Findings said "To be updated after first run" — it was never actually executed, so the manifest's `status: active` and report doc were aspirational, not real. Running it now surfaced two bugs: (1) `parse_traces()` globbed `root_path.glob('**/*trace*.txt')` — a fully recursive walk from each run root, which also contains `work/`, `rnaseq_reads/` and other multi-hundred-GB trees, so it never finished inside a 100s timeout; (2) `run.sh` passed `--roots .` (relative to `analysis/nextflow_memory_profile/`, which has no trace files) instead of `../..` (the actual repo root where `logs/nextflow/*trace.txt` lives) — so even a successful run would have silently profiled nothing from the root pipeline.

**Why it matters**: A scaffolded-but-never-run analysis with a "status: active" manifest entry and a doc that reads like a real report is worse than no analysis — it looks authoritative. Also a reminder that unbounded recursive globs (`**`) from a pipeline run root are a common slow-path trap in this repo, since run dirs always contain huge `work/`/`rnaseq_reads/` trees.

**Resolution**: Replaced the recursive glob with explicit non-recursive patterns (`logs/nextflow/*trace*.txt`, plus `logs_archive/*/logs/nextflow/*trace*.txt` for legacy copies), fixed `run.sh`'s root to `../..`, and added a `--min-realtime-min` filter (used at 30 min for `run.sh`) to exclude fast-path/no-op task runs from the RSS percentile stats — e.g. `FUNANNOTATE_TRAIN`'s "already trained, just extract shared files" early-exit branch, which otherwise dilutes the real memory requirement (unfiltered p90 RSS 989 MB vs 1.33 GB after excluding 977/1891 short runs). Full run now completes in ~1s (was timing out at 100s+) and populated `outputs/memory_profile_summary.csv` and the doc's Key Findings for real, from 15 trace files / 10,904 deduped task hashes.

**Tags**: nextflow, trace-analysis, memory-profiling, bugfix, glob-performance, funannotate-train, silent-failure

### [2026-07-23] Highest-BUSCO representative selection can pick a genome with zero ANI coverage, orphaning the whole reuse cluster

**Category**: bugfix

**What happened**: While implementing `nextflow/bin/species_reuse_clusters.py`
(`todo/species_level_abinitio_reuse.md`), the first version of representative selection
ranked candidates purely by BUSCO completeness (then N50, then alphabetical tiebreak) with
no regard for ANI data availability. Dry-running against Aspergillus fumigatus picked
`Aspergillus_fumigatus_F2` (BUSCO 99.3%, the highest in the group) as representative, and
the resulting `abinitio_reuse_assignments.csv` showed 0/375 other strains `reuse_eligible`
— every single one fail-closed on "no ANI pair found". Root cause: `F2`'s ASMID
(`GCA_051942955.1_DHU_Fsp_f2_v1`) has **zero rows in `ani.db`** — it's a genome that was
BUSCO-scored after the last `compare_ANI.nf` run and was never ANI-sketched. Since ANI is
only looked up *to the representative*, picking an ANI-uncovered genome as representative
orphans the entire cluster even when most other strains have excellent mutual ANI coverage
among themselves (confirmed: after the fix, the real representative `Aspergillus_fumigatus_Z5`
gave 373/375 strains reuse-eligible at 99.14-99.93% ANI).

**Why it matters**: This is a silent, plausible-looking failure — the script ran without
error and produced a complete, well-formed CSV; only inspecting the eligibility numbers
(0/375) revealed something was wrong. A representative-selection rule that only accounts
for assembly/annotation quality (BUSCO, N50) and ignores whether the *reuse-gating dataset
itself* (ANI) actually covers that candidate will silently degrade to "nothing qualifies"
whenever a high-quality-but-recently-added genome outranks everything else on quality alone.

**Resolution**: `pick_representative()` now restricts candidates to those with at least one
within-species `ani.db` pair to another strain in the same group *before* ranking by
BUSCO/N50/tiebreak, falling back to the full (quality-ranked) candidate set only if none of
the candidates have any ANI coverage at all (the genuinely-species-not-yet-ANI-covered case,
e.g. Beauveria bassiana before its pending SPECIES-level ANI run — there every strain will
fail-closed regardless of which one is nominally "the representative", so the fallback
doesn't mask anything).

**Tags**: funannotate, ani, species_reuse_clusters, representative-selection, silent-failure, bugfix, ab-initio-reuse

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Any future "pick the best X" selection logic that feeds
into a downstream join against a second, independently-populated dataset (here: BUSCO quality
ranking feeding into an ANI lookup) should filter candidates to those covered by the second
dataset before ranking — otherwise the highest-ranked candidate by the first dataset alone
can silently be absent from the second, and every downstream consumer fails closed.

### [2026-07-23] do_ANI/run_ANI.sh never loads the singularity module; ANI.nf silently relies on ambient shell environment

**Category**: bugfix

**What happened**: The first `--compare SPECIES` run (Aspergillus, for
`todo/species_level_abinitio_reuse.md`'s Beauveria pilot prep) hit widespread
`SKANI_COMPARE` failures — exit 127, `.command.err`: `env: 'singularity': No such file or
directory`. `profile_ANI.config` sets `singularity.enabled = true`, but `do_ANI/run_ANI.sh`
only runs `module load nextflow` before `nextflow run ...` — it never loads `singularity`,
even though the cluster's module system has it (`module avail singularity` lists
`singularity-ce/3.9.3`/`3.9.3`/`4.3.2`) and it's not in default `PATH`
(`which singularity` finds nothing on a fresh shell). Checked every `run_*.sh` under
`nextflow/`: none of them load `singularity` either. The existing 263,560-row `ani.db`
proves a GENUS-level run succeeded before, which only works if `sbatch`'s default
environment-export happened to inherit `singularity` from whatever interactive shell
submitted that job at the time — an ambient dependency the script never guaranteed.

**Why it matters**: A launcher script that "usually works" because of an inherited shell
module, with no explicit `module load` for a hard dependency the config declares
(`singularity.enabled = true`), is a silent time-bomb — it succeeds for whoever happens to
have the right modules loaded interactively and fails opaquely (a bare `env: ... No such
file` inside a Nextflow task, not an obvious top-level error) for anyone who doesn't. This
is the third distinct silent/opaque-failure gotcha found in the ANI pipeline area this
month (see the two 2026-07-19 entries above: reserved-keyword SQL DDL, names-glob mismatch)
— worth treating as a pattern, not three unrelated one-offs.

**Resolution**: Added `module load singularity` to both `do_ANI/run_ANI.sh` and
`nextflow/run_ANI.sh` (top-level, reached via `do_annotation/nextflow` -> `../nextflow`
symlink) right after `module load nextflow`/`source /etc/profile.d/modules.sh`. Other
`run_*.sh` launchers under `nextflow/` were not audited for the same gap — tracked as
T-007 in `todo/TODO_REGISTRY.md`.

**Tags**: nextflow, singularity, do_ANI, run_ANI, environment, module-load, silent-failure, bugfix

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Any `run_*.sh` SLURM launcher for a profile with
`singularity.enabled = true` (or any other hard container/env dependency) should explicitly
`module load` it rather than relying on inherited shell environment — audit all `run_*.sh`
scripts under `nextflow/` for the same gap, not just the two ANI launchers found here.

### [2026-07-24] AUGUSTUS requires a species directory's basename to match its parameter-file prefixes — the ab-initio shared store violated this and broke the first real validation run

**Category**: bugfix

**What happened**: The first real end-to-end validation of `todo/species_level_abinitio_reuse.md`
(T-004) — `funannotate predict -p <shared parameters.json>` for `Aspergillus_fumigatus_47-10`
— failed on its first attempt: `"ERROR: augustus --proteinprofile test failed, likely a
compilation error. This is required to run BUSCO, exiting."` A Nextflow automatic retry
happened to succeed (funannotate resumes from `predict_misc/` checkpoints), but the
resulting annotation was measurably worse than the strain's existing independently-trained
legacy annotation: BUSCO completeness 95.0% vs. 96.9% (missing BUSCOs nearly doubled,
17→31 of 1122), and 9,078 vs. 9,878 total gene models (−8.1%). Root cause, found by direct
comparison against a known-working stock AUGUSTUS species (`anidulans`): **AUGUSTUS requires
a species directory's basename to exactly match the prefix of the parameter files inside
it** — `species/anidulans/anidulans_parameters.cfg`, `species/anidulans/anidulans_exon_probs.pbl`,
etc. `nextflow/bin/species_reuse_clusters.py`'s original backfill copied the representative
strain's AUGUSTUS species directory to a species-level directory name
(`_shared_abinitio/Aspergillus_fumigatus/augustus/species/Aspergillus_fumigatus/`) **without
renaming the files inside it**, which still carried the representative's own strain-specific
prefix (`aspergillus_fumigatus_z5_parameters.cfg`, etc.) — a basename/file-prefix mismatch
that silently broke AUGUSTUS's species lookup. The retry's "success" was actually a partial
failure: AUGUSTUS's own trained-parameter contribution to the gene set was effectively
unusable, and EVM fell back on the other predictors, degrading the final consensus.

**Why it matters**: This was caught only because a real predict run happened and a careful
before/after BUSCO comparison was done — the pipeline's own retry-on-failure logic actively
*hid* the first, cleaner failure signal by resuming from a partial checkpoint state and
appearing to succeed. A less careful validation (checking only "did predict exit 0 and
produce a GBK") would have shipped a silently degraded feature. This is exactly the kind of
result the plan's S5 "run a real validation pilot before trusting this broadly" step exists
to catch — vindicates not skipping it, and argues for treating a validation run's first-vs-
retry attempt history as informative, not something to discard once a later attempt succeeds.

**Resolution**: Redesigned the shared store: (1) relocated from
`<target>/_shared_abinitio/<species_tag>/` to a new top-level location,
`gene_prediction_shared_abinitio/<species_tag>/` (a fixed absolute path, not nested under
`params.target`, since it's meant to be shared across every project tree that annotates the
same species — `do_annotation/`, `asco_annotate/`, etc. — not scoped to one run); (2) the
copied AUGUSTUS parameter files are now renamed on copy from the representative's own prefix
to `<species_tag.lower()>_*`, matching AUGUSTUS's convention exactly; (3) the resulting
species directory is also symlinked into the real `AUGUSTUS_CONFIG_PATH`'s own `species/`
directory (`lib/augustus/3.5/config/species/<species_tag_lower>` -> the shared-store copy),
so it's discoverable via the ordinary `--augustus_species <name>` mechanism too, not just the
absolute path embedded in `parameters.json` — and reuses the config tree's existing
`cgp/extrinsic/model/parameters/profile/` subdirectories rather than needing to duplicate
them per species. Verified the fix directly: compared the new store's file-naming structure
against the known-working `anidulans` stock species (structurally identical now), and
confirmed `augustus --species=<name>` no longer errors on species lookup (only on the
separate, expected "no query file" condition). Both Aspergillus fumigatus and Beauveria
bassiana stores regenerated at the new location/naming; the old, incorrectly-named
`genome_annotation/_shared_abinitio/` removed.

**Tags**: augustus, funannotate, ab-initio-reuse, species-naming-convention, silent-failure, bugfix, validation

**mitigation_type**: code-fix

**structural_mitigation_candidate**: Any future code that copies/relocates a tool's own
"species"/"model"/"profile" config directory to a new name should treat the tool's own
internal file-naming convention as a hard constraint to verify, not assume renaming the
container directory alone is sufficient — check for a real, known-working example of that
tool's own convention (as this fix did with `anidulans`) before trusting a redesigned
directory layout.

### [2026-07-24] AUGUSTUS's *_parameters.cfg hardcodes sibling filenames in its own content, not just on disk — renaming files alone still left a stale reference that crash-looped a validation run for 4+ hours

**Category**: bugfix

**What happened**: The file-renaming fix above (same day) turned out to be necessary but
not sufficient. A subsequent clean validation rerun for the same strain
(`Aspergillus_fumigatus_47-10`) appeared to "hang" on the AUGUSTUS gene-prediction step for
4+ hours with no new log output — not actually hung, but crash-looping silently: AUGUSTUS
was failing on *every* genome chunk with `"Couldn't open the file with the weight matrix:
aspergillus_fumigatus_z5_weightmatrix.txt"` (visible only in `logfiles/augustus-parallel.log`,
a different, more granular log file than `logfiles/funannotate-predict.log`, which showed no
error at all — the top-level log just looked stalled). Root cause: `*_parameters.cfg` itself
contains explicit internal directives referencing sibling filenames by name —
`/BaseCount/weightMatrixFile`, `/ExonModel/infile`/`outfile`, `/IntronModel/infile`/`outfile`,
`/IGenicModel/infile`/`outfile`, `/UtrModel/infile`/`outfile` — all still pointing at the
representative strain's original prefix (`aspergillus_fumigatus_z5_*`) even after every file
on disk had been correctly renamed to the new species-level prefix. The `.cfg` file's own
comment literally says *"change this to your species if at all necessary"* — AUGUSTUS
expects the *content*, not just the filename, to be updated when repurposing a species
profile under a new name. Confirmed no other `.cfg` file in the species directory
(`_metapars*.cfg`) has this issue — only `_parameters.cfg` cross-references siblings.

**Why it matters**: This is a second, independent way the same underlying problem (reusing
a species profile under a different name) can break — fixing the filename convention alone
looked complete (verified structurally against a working stock species, `anidulans`) but
missed a content-level cross-reference that only manifests as a runtime crash-loop, not a
parse error or missing-file error at startup. It also demonstrates a debugging-process
lesson directly from the user: `logfiles/funannotate-predict.log` alone made this look like
an indefinite hang; the actual, immediately-diagnostic error was sitting in a sibling log
file (`logfiles/augustus-parallel.log`) the whole time. **When a funannotate stage appears
stalled with no new top-level log output, check every file in that stage's `logfiles/`
directory, not just the main one, before assuming a hang.**

**Resolution**: `species_reuse_clusters.py`'s AUGUSTUS file-copy loop now text-substitutes
the representative's old lowercased prefix for the new species-level prefix inside `.cfg`
file *content* (not just the filename) before writing the copy. Verified with `grep -c
<old_prefix>` against both regenerated stores (Aspergillus fumigatus, Beauveria bassiana) —
zero stale references in either. Not yet re-verified with an actual completed predict run
at time of writing (the crash-looping job was killed once the log was found, before a
corrected rerun completed) — next validation run against the corrected store is the real
test of whether this (combined with the filename fix) fully resolves it.

**Tags**: augustus, funannotate, ab-initio-reuse, config-file-content, silent-failure, crash-loop, debugging-process, bugfix

**mitigation_type**: code-fix

**structural_mitigation_candidate**: When repurposing any tool's config/profile files under
a new identity (not just AUGUSTUS), grep the file *contents* for the old identifier, not
just check filenames — config formats commonly embed self-referential filenames or names
inside directives/headers that a pure filesystem rename won't touch. And structurally: when
investigating an apparently-hung long-running step, check ALL files in its logs directory
(`ls -t logfiles/ | head`, or grep across all of them) before spending significant wall-clock
time assuming it's merely slow.

### [2026-07-23] SKANI/MASH/SOURMASH/FASTANI_COMPARE's storeDir never re-runs a group's ANI comparison after new genomes are added to it

**Category**: bugfix

**What happened**: While debugging why `Aspergillus_fumigatus_F2` (samples.csv row added
after the last GENUS-level `compare_ANI.nf` run) never appeared in `ani.db`, found that
`results/ANI/skani/GENUS/Aspergillus/batches/Aspergillus.full.ani.tsv` had an mtime of
2026-06-20 — over a month *before* F2's sketch was even computed (2026-07-19). Multiple
`-resume` re-runs since then never refreshed it. Root cause: `SKANI_COMPARE` (and identically
`MASH_COMPARE`, `SOURMASH_COMPARE`, `FASTANI_COMPARE` in `compare_ANI.nf`, plus
`SKANI_DIST_QUERY` in `query_ANI.nf`) used `storeDir` for their per-group comparison output.
`storeDir` only checks whether the declared output file already exists at that path — unlike
Nextflow's normal task-hash-based caching, it does **not** re-evaluate whether the task's
actual inputs (here: the group's current sketch list) have changed. So once a genus/species
group's `.full.ani.tsv` exists, adding new genomes to that group later (a normal, expected
occurrence as `samples.csv` grows) silently never triggers recomputation, no matter how many
times the pipeline is re-run with `-resume` — the new genome gets correctly sketched
(`SKANI_SKETCH` is a separate, correctly-behaving per-genome `storeDir`) but never compared
against anything in that group. Confirmed via `report_ani.py`'s own output on a *fresh*
SPECIES-level run (which did re-execute, since that combination of group+rank had never been
computed before): F2 showed up correctly as a processed, explicitly-flagged singleton/outlier
(`best_ANI=0.00%`) — proving the pipeline logic itself handles a genuinely-uncomparable
genome correctly; the bug was specifically storeDir masking staleness for a group that had
*already* been computed once.

**Why it matters**: This is a structurally different failure mode from a task *failure*
(which `-resume` + normal caching handles correctly, confirmed separately with the
singularity-module fix earlier this session) — this is silent staleness from **data growth**,
which `-resume` cannot fix by design once `storeDir`'s existence-check is satisfied. Every
GENUS/SPECIES-level ANI comparison for a genus that has ever been run before will silently
exclude any genome added to `samples.csv` afterward, with no error, warning, or visible
symptom other than a missing/zero-pair genome in `ani.db` — exactly the kind of pipeline
"success" that's actually wrong that the mycelium robust-analysis convention warns about.

**Resolution**: Changed all five processes (`SKANI_COMPARE`, `MASH_COMPARE`,
`SOURMASH_COMPARE`, `FASTANI_COMPARE` in `compare_ANI.nf`; `SKANI_DIST_QUERY` in
`query_ANI.nf`) from `storeDir "..."` to `publishDir { "..." }, mode: 'copy'` (closure syntax
required for publishDir paths that reference per-task input variables like `group_name` —
plain-string interpolation only works for `storeDir`, and using it for `publishDir` throws
`No such variable: group_name` since Nextflow evaluates a non-closure publishDir path once at
process-definition time, without task-input scope). This matches `REPORT_ANI`'s own existing
`publishDir { ... }` convention already in the same file. Effect: Nextflow now hashes the
actual task inputs (the sketch/genome list) on every run, so a group whose membership changed
since the last run is correctly detected as needing re-execution, while an unchanged group
still gets skipped via normal `-resume` caching (no loss of the caching benefit `storeDir` was
originally used for). Manually force-regenerated the stale `Aspergillus.full.ani.tsv` by
deleting it (storeDir output must be removed by hand once; the new publishDir-based caching
prevents this class of staleness going forward). Verified both `compare_ANI.nf` and
`query_ANI.nf` still parse and `-preview`-run cleanly after the change.

**Tags**: nextflow, storeDir, publishDir, caching, ani, compare_ANI, query_ANI, silent-failure, bugfix, staleness

**mitigation_type**: code-fix

**structural_mitigation_candidate**: `storeDir` should be treated as suspect by default for
any process whose semantic inputs are a *set* that can grow over time (a group's genome
list, a species' strain list, etc.) — it is only safe for processes whose output is a pure,
permanent function of immutable per-item inputs (e.g. `SKANI_SKETCH`'s one-sketch-per-genome
output, which is correctly safe since a genome's own sketch never changes once computed).
Audit the rest of the codebase for other `storeDir` usages gating on a *collected/grouped*
input rather than a single immutable item.

### [2026-07-23] COMBINE_ANI_TABLE's merged ani.db reflects only whichever group(s) were active in the latest run, not the full historical set — OPEN INVESTIGATION, not yet resolved

**Category**: bug (unresolved)

**What happened**: After the `storeDir`->`publishDir` fix above (same day, same session) and a
successful Beauveria SPECIES-level + GENUS-level ANI computation (confirmed correct: 24,753
real within-species pairs, per-group `Beauveria_ANI_report.txt`/`batches/Beauveria.full.ani.tsv`
present and internally consistent), the GENUS-level `results/ANI/skani/GENUS/ani.db` was
checked after the pipeline's next full run (which recomputed the previously-deleted, stale
`Aspergillus.full.ani.tsv`, now covering 1,541 genomes) completed. The resulting `ani.db`
contained **263,560 rows, but `SELECT COUNT(DISTINCT query_genus)` = 1 ("Aspergillus" only)**
— Beauveria's correctly-computed group data never made it into the merged table, despite its
own per-group files being present and correct on disk in `results/ANI/skani/GENUS/Beauveria/`.

**What's confirmed, not guessed**:
- Beauveria's own `SKANI_COMPARE`/`REPORT_ANI` outputs are correct and present on disk,
  independent of the combine step.
- Aspergillus's recomputation (1,541 genomes) completed correctly this run.
- `input_clean_genomes/` has ~46,172 files (roughly 2 per genome: `.fa.gz` +
  `.masked.fasta.gz`), consistent with the full ~22,930-genome `samples.csv`, not just
  Aspergillus/Beauveria — ruled out "missing genome input files for other genera" as the
  cause.
- `COMBINE_ANI_TABLE`'s inputs are `ani_tsv_ch.map{...}.collect()` and
  `names_map.map{...}.collect()` (compare_ANI.nf:708-711) — a live `.collect()` over the
  *current run's* channel emissions, published via `publishDir ..., mode: 'copy'`
  (unconditional overwrite, no merge-with-existing-file logic). This process was **not**
  touched by today's `storeDir`->`publishDir` fix (it already used `publishDir`).

**What's genuinely unresolved and needs its own investigation, not assumed**:
- Whether Nextflow's `-resume` reliably re-emits **cached** (unchanged-input, skipped) groups'
  prior outputs into `ani_tsv_ch`/`names_map` so `.collect()` can gather the full historical
  set alongside freshly-recomputed groups, or whether only actively-scheduled groups in that
  invocation flow into the channel. This is the crux mechanical question and was not directly
  confirmed either way this session.
- Whether the *original*, session-start `ani.db` (the one the 2026-07-19 blank-columns bugfix
  was verified against) ever actually had genuine multi-genus coverage, or was *already*
  effectively single-genus and this was never caught because nobody queried
  `COUNT(DISTINCT query_genus)` before — only `query_species`/blank-column checks were done
  at the time. Not verifiable retroactively without an independent record of that file's
  content.
- Behavior under back-to-back/concurrent invocations of `do_ANI/run_ANI.sh` (this session
  alone launched it at least 5 times within ~90 minutes, 19:44-20:46, while a previous
  invocation's tasks may still have been in flight) — whether `COMBINE_ANI_TABLE` from one
  invocation can race with or get clobbered by another, and whether Nextflow's own
  `-resume` session/cache database (`.nextflow/cache/`) is safe under that usage pattern for
  this pipeline. Not tested.

**Resolution**: Implemented the disk-glob approach, then hit a *second*, subtler bug while
validating it: `gatedGlobIn` (ported from `BFD.nf`, see the `storeDir`->`publishDir` entry
above) globs the *published* location, but `publishDir`'s copy is asynchronous relative to a
process's channel emission (a documented Nextflow behavior — publishDir is a side effect, not
a synchronization point). Gating the glob purely on "this run's channel emitted" let the glob
race ahead of `publishDir`'s own copy for a group *this run just finished* — confirmed via a
`-stub-run --n_test 3 --taxon GENUS:Beauveria` test where `COMBINE_ANI_TABLE` never fired at
all (`completed=10`, missing entirely from the summary) because the glob found zero files.

Final fix: **union** the disk glob with the live channel directly
(`gatedGlobIn(...).mix(ani_tsv_ch.map{...}.collect())`, same for names via
`REPORT_ANI.out.names`) — the live channel item is available immediately with no publish
race and guarantees *this* run's groups are always captured, while the glob still supplies
full historical coverage for groups untouched this run. The same group can now legitimately
appear in the manifest via both paths, so `combine_ani_table.py`'s `read_manifest()`/
`dedupe_by_group()` keep the entry with the latest manifest mtime per group — verified with a
standalone test (two paths for the same group, old content vs. fresh content) that the fresh
one always wins and no group is double-counted. Re-ran the `-stub-run` test after this fix:
`completed=11` (now includes `COMBINE_ANI_TABLE`). `REPORT_ANI` gained explicit
`emit: report`/`emit: names` labels to support this (previously unlabeled, positional-only
output). Both `compare_ANI.nf` and `combine_ani_table.py` changes verified against real
Beauveria data on shared storage (not `/tmp` — an earlier test attempt also hit a
node-local-`/tmp`-isn't-shared-across-SLURM-nodes failure, unrelated to this bug, worth
remembering for future ad hoc pipeline tests on this cluster).

Genuinely still open, not tested: whether two `do_ANI/run_ANI.sh` invocations running
truly concurrently (not just back-to-back) from the same working directory are safe —
Nextflow's own session-lock (`.nextflow/cache/<id>/db/LOCK`, tied to CWD not `-work-dir`)
was observed to correctly reject a second concurrent `nextflow run` from the same directory
during this session's testing, which suggests full concurrent execution is not actually
possible from one CWD (a second attempt just fails fast with a clear lock error) — but this
was an incidental observation, not a deliberate test, and wasn't the scenario this fix
targets (sequential back-to-back invocations, which *are* the scenario this fix resolves).

**Tags**: nextflow, ani, compare_ANI, combine_ani_table, resume, collect-operator, publishDir-race, data-loss-risk, bugfix, session-lock

**mitigation_type**: code-fix

**mitigation_type**: code-fix

**structural_mitigation_candidate**: When an analysis folder's `outputs/` is empty and its report doc says "to be updated after first run," treat that as a signal to actually run it before trusting its manifest status or citing its (nonexistent) findings.

## funannotate predict's glimmerhmm weight=0 does not skip its BUSCO-seeded training (2026-07-25)

While testing T-013 (`todo/species_level_abinitio_reuse.md`'s idea of skipping
`FUNANNOTATE_TRAIN`/PASA entirely and feeding `funannotate predict --transcript_evidence`
the raw Trinity-GG assembly directly, combined with T-004's `-p <shared_params.json>`
ab-initio reuse for augustus/genemark/snap), the run still paid for a full genome-mode
BUSCO step (~10-25 min) despite passing `-w codingquarry:0 glimmerhmm:0` — the same weight
flag production `FUNANNOTATE_PREDICT` always passes. Root cause (confirmed by reading
`predict.py` source directly, not guessing): `-w glimmerhmm:0` only zeroes glimmerhmm's EVM
consensus weight. `RunBusco` is set `True` whenever *any* predictor's `RunModes[...] ==
"busco"` (predict.py ~line 1051-1058), independent of that predictor's weight. Our shared
ab-initio store (`gene_prediction_shared_abinitio/<species>/parameters.json`) only has
`augustus`/`genemark`/`snap` pretrained — `glimmerhmm` is `[{}]` (empty) — so glimmerhmm
always needs *some* training mode, and without PASA gene models available it falls back to
`"busco"`, triggering the expensive genome BUSCO run regardless of its zero weight.

**Production T-004 runs are NOT affected**: with `FUNANNOTATE_TRAIN`'s PASA data present,
glimmerhmm trains via the cheap `"pasa"` mode instead (confirmed directly in
`genome_annotation/Aspergillus_fumigatus_47-10/logfiles/funannotate-predict.log`:
`{'augustus': 'pretrained', 'genemark': 'pretrained', 'snap': 'pretrained', 'glimmerhmm':
'pasa'}`). This is specific to no-TRAIN (T-013-style) runs, where there is no PASA data at
all to give glimmerhmm a cheap fallback.

**How to apply**: for any future no-TRAIN / `--transcript_evidence`-only profiling run
intended to measure the true minimum wall-clock (no from-scratch ab-initio training at
all), either (a) add a pretrained `glimmerhmm` entry to the species' shared
`parameters.json` (same backfill mechanism `species_reuse_clusters.py` already uses for
augustus/genemark/snap), or (b) accept the BUSCO-seeding cost as a known, documented part
of the no-TRAIN tradeoff rather than mistaking it for a bug. Weight=0 flags in funannotate
predict generally control EVM consensus contribution only, not whether a predictor is
trained — don't assume weight 0 is equivalent to "skip this predictor" for cost purposes.

**Tags**: funannotate, predict, glimmerhmm, ab-initio-reuse, T-013, T-004, busco-fallback, runtime-profiling

## Raw background processes on the login node can OOM-kill a concurrently-running nextflow driver (2026-07-25)

While running the 12-strain T-004/T-013 comparison batch, launched 9 concurrent `nohup
busco ... &` background processes directly on the login/interactive node (each BUSCO
spawning many hmmsearch/multiprocessing workers) at the same time condition A's sequential
`nextflow run` batch (`do_annotation/.claude/run_abinitio_validation_batch2.sh`) was
running its own lightweight orchestrator process on the same node. The resource spike
killed the nextflow driver process mid-strain (`Killed`, confirmed in the batch log) — and,
separately, all 9 of the raw background BUSCO processes were also lost (no output produced
for any of them), suggesting the whole process group was reaped, not just one process.

A second bug compounded the silent failure: the batch script's completion log line was
`echo "[$(date)] Finished ${out} (exit $?)"` — `$?` at that point reflects `date`'s exit
status (always 0), not nextflow's, since `$(date)` runs and resolves *before* the `$?`
expansion. So the log claimed `exit 0` for a run that was actually killed. Confirmed via
the trace file: the killed strain's trace was 52 bytes (header only, no TRAIN/PREDICT rows)
vs 228+ bytes for a real completed run.

**Production data was NOT corrupted** — the driver died before `FUNANNOTATE_PREDICT`'s
stale-clear step (`rm -rf predict_results predict_misc`) executed, so the existing GBK was
left untouched. Lucky, not guaranteed — a kill landing slightly later (mid-clear, or during
a real repredict with cleared old data and no new data written yet) would have destroyed
production data with nothing to restore it from.

**How to apply**:
1. Always capture a command's exit status immediately (`cmd; status=$?`) before running
   anything else (even `date` inside the same echo) — never trust `$?` after any
   intervening command substitution.
2. Don't run resource-heavy validation jobs (BUSCO, etc.) as raw background processes on a
   shared login/interactive node while a nextflow orchestrator (or any other long-running
   driver process) is active on the same node — submit them via `sbatch` instead so SLURM's
   scheduler enforces resource limits without competing unbounded for the login node's
   memory/CPU. See `do_annotation/.claude/busco_validation/busco_run.sbatch` (added this
   session) for the corrected pattern.
3. A sequential batch driver (required here by the `.nextflow` session lock tied to CWD —
   see earlier learnings) should hard-stop and surface an error on any non-zero exit
   instead of silently continuing to the next strain, so a killed strain doesn't get
   miscounted as done.

**Tags**: nextflow, slurm, oom, background-process, exit-status-bug, T-004, T-013, bugfix, resource-contention

### [2026-07-27] Never call Channel.fromPath() inside a Nextflow operator closure — it deadlocks

**Category**: gotcha

**What happened**: `COMPARE_ANI:CONCAT_ANI_TSVS` stayed ACTIVE with its `manifest` input queue permanently OPEN. Upstream `SKANI_COMPARE` and `REPORT_ANI` had completed, but the manifest channel feeding `CONCAT_ANI_TSVS` never emitted. The manifest was built inside a `.map { ... }` closure by calling `channel.fromPath(...).filter { ... }.collect().val`. In DSL2, `channel` (lowercase) is an unreliable alias inside closures; even with uppercase `Channel`, creating a new channel factory inside an operator closure and blocking on `.val` does not participate in the outer execution graph and hangs forever.

**Why it matters**: This pattern was present in two places: `nextflow/workflows/compare_ANI.nf` and `nextflow/subworkflows/local/ANI_REPRESENTATIVE_SELECT/main.nf`. Both feed `CONCAT_ANI_TSVS` / `PICK_REPRESENTATIVE_STRAIN`, so any ANI-driven representative-selection run would stall indefinitely.

**Resolution**: Replaced the closure-internal globbing with collection of the live `ani_tsv` channel (`.map { _group, tsv -> tsv }.collect().map { ... }.collectFile(...)`). This avoids the deadlock and also removes a race against `publishDir` copying. Where synchronous filesystem listing is truly needed, use Nextflow's built-in `files()` helper instead of `Channel.fromPath()`.

**Tags**: nextflow, channel-semantics, deadlock, Channel.fromPath, files(), ani, compare_ANI, CONCAT_ANI_TSVS, operator-closure

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped — replaced both occurrences. Repo convention: do not use `Channel.fromPath()` (or `channel.fromPath`) inside `.map`/`.filter`/operator closures; collect the relevant channel or use `files()` for synchronous globs.

### [2026-07-28] Unescaped `$` in a Nextflow `"""..."""` script block silently mangled the shell command — CONCAT_ANI_TSVS retried 3x and failed with exit 123

**Category**: gotcha

**What happened**: `COMPARE_ANI:CONCAT_ANI_TSVS` failed with exit status 123 on every retry (`.command.sh: line 7: unexpected EOF while looking for matching '''`). The `script:` block (a Groovy GString, `"""..."""`) contained a literal shell regex `grep -v '^$'` intended to drop blank lines when building the asmid manifest. Because the block is interpolated by Groovy, not passed through verbatim, the unescaped `$'` was consumed during GString parsing and the rendered `.command.sh` came out as `grep -v '^ |` — a truncated, syntactically broken pipeline. This wasn't caught at pipeline-compile time because `'^$'` parses as a valid (if wrong) Groovy expression.

**Why it matters**: Any process script that embeds shell/regex syntax containing a bare `$` (blank-line regexes, `awk '{print $1}'`, `"$HOME"`, positional shell params) is at risk of this exact silent-mangling failure mode. It costs real wall-clock time here — 3 SLURM submissions (~4 min) before the pipeline gave up.

**Resolution**: Rewrote `nextflow/modules/ani/report/CONCAT_ANI_TSVS/main.nf` to do the asmid dedupe/sort/blank-filtering in Groovy (`asmids.findAll{it}.collect{it.trim()}.findAll{it}.unique().sort()`) before the script string, then interpolate the finished list into a quoted heredoc (or `touch` if empty). This removes the `printf | tr | grep | sort` shell pipeline and the `$`-escaping requirement entirely, and is also more correct for asmid strings containing whitespace (the old `tr ' ' '\n'` would have exploded a multi-word asmid). Added a rule to `~/.claude/skills/nextflow-hpcc/SKILL.md` and this repo's `CLAUDE.md` to audit every bare `$` inside `"""..."""` script blocks before shipping.

**Tags**: nextflow, groovy, gstring, escaping, ani, compare_ANI, CONCAT_ANI_TSVS, exit-123, bugfix

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped — documented in `nextflow-hpcc` skill and project `CLAUDE.md`. Repo convention: any literal shell/regex `$` inside a Nextflow `"""..."""` script block must be escaped as `\$`, or (preferred) computed in Groovy before the script string so there's no `$` to escape.

### [2026-07-29] `join(by:, remainder:true)` only satisfies ONE left-side item per key — wrong operator for a one-to-many gate

**Category**: gotcha

**What happened**: `FUNANNOTATE_PREDICTION.nf`'s `--predict_scope all` gating (does an eligible non-representative sibling's species have shared ab-initio params available yet?) used `eligibleSiblingKeyed.join(availableSpeciesCh, by: 0, remainder: true)`, keyed by species. With a real test case of two eligible siblings of the same species (`Neurospora crassa`), only one (`Neurospora_crassa_73`) was correctly matched against the single "species available" signal; the other (`Neurospora_crassa_FGSC_2489`) was silently treated as unmatched and routed to the blocked-species report, even though its species genuinely had shared params available. Confirmed directly with a standalone probe (`left = [tuple('sp1','row-A'), tuple('sp1','row-B'), tuple('sp2','row-C')]`, `right = [tuple('sp1', true)]`, `.join(by:0, remainder:true)`) — output was `[sp1, row-A, true]`, `[sp1, row-B, null]`, `[sp2, row-C, null]`. `row-B` shares `sp1`'s key with `row-A` but still came back `null`, identical to the genuinely-unmatched `sp2`/`row-C`.

**Why it matters**: `join()` implicitly assumes a one-to-one (or many-to-one on the *left*, one-to-one on the *right*) relationship per key, consuming the right-side item once. Any gating problem shaped "one key can have MANY per-item entries waiting on it" (a species with several sibling strains; a group with several members) is the wrong shape for `join()`, even with `remainder:true` — it will silently misclassify every duplicate-keyed item after the first as unmatched, with no error or warning.

**Resolution**: Replaced with `perItemChannel.combine(collectedSet)` where `collectedSet` is built via `.collect().map { it as Set }.ifEmpty([] as Set)` on the availability channel. Probed directly and confirmed `combine()` broadcasts a *collected* Set/List value as one atomic tuple element to every left-side item, including duplicate keys — `left.combine(right.collect().map{it as Set})` on the same left fixture produced `[sp1, row-A, [sp1]]`, `[sp1, row-B, [sp1]]`, `[sp2, row-C, [sp1]]` — both `sp1` rows got the broadcast set, `sp2` correctly did not. This is a different shape from the earlier documented `combine()`-spreads-list-valued-items bug (that one was two *collected-list* channels combined together, one possibly empty, causing degenerate flattening); here one side is a real per-item channel (never collected) and the other is a single, fully-collected scalar value — the documented "safe" `combine()` shape.

**Tags**: nextflow, channel-semantics, join, combine, one-to-many, gotcha, funannotate, FUNANNOTATE_PREDICTION, ab-initio-reuse, operator-probe

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped in `FUNANNOTATE_PREDICTION.nf`. Repo convention: for any "gate N per-item channel entries on a key that's only fully known once an upstream stream closes," use `combine()` with a `.collect()`-ed Set/Map value on one side, never `join(..., remainder: true)`, whenever more than one per-item entry can share the same key. When in doubt about `join`/`combine`/`collectFile` behavior with duplicate keys or empty/collected inputs, write a tiny standalone probe script (`Channel.of(...).join/combine(...).view { ... }`) rather than reasoning about it from documentation alone — this is now the second time in this repo doing so caught a real, silent misbehavior.

## An O(#candidates) matching algorithm can look fine on the login node and stall for an hour on a slow compute node (2026-07-30)

**Category**: gotcha / performance

**What happened**: `scripts/one-off/reorg_genome_stats_hash_buckets.py`'s `longest_prefix_match()` tried every candidate ASMID/basename (up to ~46k/~31k) as a literal string prefix of each filename -- O(#candidates) comparisons per file. Interactive dry runs on the login node (`r11`) completed the full 171k-file scan in ~3 minutes, giving no hint of a problem. The first real production `--apply` run (job 26992449) was submitted via plain `sbatch -p short` with no node constraint, landed on `c18` (the `short` partition spans several CPU generations, including the old `c[01-30]` `amd,abu_dhabi` nodes alongside much newer `ryzen`/`cascade` nodes), and stalled for ~57 minutes with the log showing zero output past the initial samples.csv load and no manifest file ever created -- i.e. it hadn't even finished the read-only planning phase, let alone started moving files.

**Why it matters**: An algorithm whose cost scales with candidate-pool size (tens of thousands here) rather than input size can be "fast enough" on a fast/idle node and silently pathological on a slow/contended one, with no error, no crash, just a job that never progresses. `sinfo -o "%N %f"` on this cluster shows the `short` partition includes `c[01-30] amd,abu_dhabi` (an old AMD generation) alongside `r[01-51] ryzen,amd,{rome,milan,genoa}` and `x[01-06] intel,cascade` -- a plain `sbatch -p short` with no `--constraint` can land on any of them.

**Resolution**: Two independent fixes, both real and both necessary going forward:
1. Algorithmic: replaced the O(#candidates) candidate-list scan with an O(len(filename)) approach -- enumerate the *filename's own* delimiter positions (typically a handful) and do a dict/set membership check at each length, instead of iterating the candidate pool. Verified: full-scope dry run on the same data went from an unbounded stall to 6.8 seconds.
2. Operational: `scripts/one-off/submit_reorg_genome_stats_apply.sh` now requests `--constraint='ryzen|cascade'` on every submission, and splits the single monolithic `--apply` invocation into one job per type (`--only <type>`), chained via `--dependency=afterok`, smallest file count first -- so a stuck/slow type doesn't consume the whole time budget, and the smallest type acts as a live canary before the largest (`gene_stats`, 55k files) commits real production time.

**How to apply**: For any future one-off script whose cost scales with `len(some_lookup_table)` per input item (not just `len(inputs)`), profile it against realistic production-scale candidate-pool sizes, not just realistic input counts -- and don't trust a clean interactive test run on the login node as a proxy for compute-node behavior under `sbatch` without a `--constraint`, since this cluster's `short` partition spans CPU generations differing by roughly a decade. Prefer dict/set membership tests over "iterate all candidates and check startswith" whenever the candidate pool can be large.

**Tags**: performance, algorithmic-complexity, slurm, sbatch, constraint, node-architecture, genome_stats, one-off-script, T-014

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped in PR #17 (algorithmic fix) and `scripts/one-off/submit_reorg_genome_stats_apply.sh` (constraint + per-type batching). Repo convention candidate: any new `sbatch` submission for a script whose runtime is uncertain/untested at production scale should default to `--constraint='ryzen|cascade'` (or profile on the intended node class first) rather than an unconstrained `-p short`.

## Migrating a directory layout also changes the filename, not just the path — easy to fix the storeDir and still break everything (2026-07-30)

**Category**: gotcha

**What happened**: While scoping issue #9 (bucket-aware `storeDir` updates, split into #19/#20/#21), an expert-review pass (Fable) caught that the naive plan — "just add a hash-bucket subdirectory to each `storeDir`" — was incomplete in a way that would have silently broken the whole point of migrating data first. The already-completed legacy migration (T-014) didn't just move files into bucket subdirectories; for every LOCUSTAG-keyed type it also *renamed* them, from the old `${meta.id}` (Genus_species_strain) basename to a `LOCUSTAG`-based one (confirmed: `results/genome_stats/gene_stats/fac/F8099293.gene_CDS.csv.gz`, not `Species_strain.gene_CDS.csv.gz`). Every current producer module (`BUSCO_GENOME`, `BUSCO_PEP`, `CALC_GENE_STATS`, `CALC_INTERGENIC`, `IDP`, `MEROPS`, `PFAM`, `PREDGPI`, `SIGNALP`, `TARGETP`, `TMHMM`, `WOLFPSORT` — 12 of 13 checked, only `CALC_ASM_STATS` was already correct) still emits `${meta.id}.<ext>` in both its `output:` declaration and its script/stub body.

**Why it matters**: If only the `storeDir` directive is updated (bucket added) and the emitted filename is left alone, a post-cutover run's output won't match the existing migrated file at the same bucket path *at all* — `storeDir`'s existence check looks for `bucket/Species_strain.ext`, finds nothing (only `bucket/LOCUSTAG.ext` exists), and every genome gets needlessly recomputed. Worse, the recompute doesn't even fix the mismatch — you end up with two different filenames for the same genome sitting in the same bucket forever. This is silent: no error, no warning, just an unnecessary full-pipeline recompute (HPC-hour-costly for BUSCO-heavy types) followed by permanent orphaned duplicates. Exactly the kind of thing a code read alone can miss, because "update storeDir to add a bucket" reads as obviously sufficient until you actually diff the migrated filename convention against the current module's emitted filename.

**Resolution**: Issues #19/#20 were rewritten to require BOTH changes per module: the `storeDir` bucket subdirectory, AND the emitted filename itself (output declaration + script/stub body) switching from `${meta.id}` to `${meta.locustag}` (or `${meta.asmid}` for the two ASMID-keyed types, `CALC_ASM_STATS`/`BUSCO_GENOME`). Their acceptance criteria now explicitly require confirming both zero-recompute *and* zero-duplicate-filename-per-genome, not just that the pipeline runs without error.

**How to apply**: Whenever a one-time data migration changes a file's *identity* (not just its location) — renaming, re-keying, reformatting — any code that will later regenerate that file must be checked for whether it reproduces the *exact* new identity, not just the new location. "Where does it go" and "what is it called" are separate questions; a migration plan/review that only asks the first one will miss this class of bug. Worth a standing check for any future storage migration in this repo.

**Tags**: genome_stats, storeDir, migration, filename-vs-path, locustag, meta.id, T-014, gotcha, fable-review

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped in issues #19/#20 (not yet implemented as code — this was caught during pre-implementation scoping, before any of the actual storeDir changes were written). Repo convention candidate: any issue describing a "make code match a completed data migration" task should explicitly diff the migration's output naming convention against the current producer's naming convention as a required scoping step, not just the directory structure.

## `beforeScript` module loads do not reliably survive into a process's `script:`/`stub:` block on this HPC — put `module load` inline instead (2026-07-31)

**Category**: gotcha

**What happened**: While implementing issue #27 (T-014 step 5, Parquet conversion), added `module load duckdb` to the `merge` label's `beforeScript` in `nextflow/conf/profile_BFD.config` (mirroring the existing `module load miniconda3; module load biopython` pattern already there) plus a matching override in `test.config`. A `-stub-run` against `-profile BFD,test` failed on every `MERGE_*` process with `duckdb: command not found`, even though `beforeScript`'s content was correctly embedded in `.command.run` and `module load duckdb` succeeds fine when run standalone in an interactive shell.

Root cause, confirmed by direct experimentation (not guessed): Nextflow's generated `.command.run` executes `beforeScript`'s content inside `nxf_main()`, a **non-login** `bash` process. `nxf_main()` then calls `nxf_launch()`, which spawns `/bin/bash -l .command.run nxf_trace` — a **brand-new login shell** — and the actual `script:`/`stub:` payload (`.command.sh`) is *also* shebanged `#!/bin/bash -l`, i.e. its own fresh login shell. On this HPC, a login shell's startup chain (`/etc/profile` → the user's own profile scripts) unconditionally reloads a **fixed baseline module list** regardless of what was exported/loaded by the parent process — confirmed via `bash -l -c 'echo $LOADEDMODULES'` showing a fixed set (`slurm`, `openmpi`, `java`, `texlive`, `pandoc`, `R`, `hpcc_user_utils`, `miniconda3`) that drops anything loaded interactively beyond that baseline (tested: `tmux`, `helix`, `neovim`, `duckdb` all vanished across the login-shell boundary; `miniconda3` "survived" only because it was already part of the baseline, not because `beforeScript` actually propagated it). So `beforeScript`'s `module load duckdb` was silently discarded before `.command.sh` ran — and the pre-existing `module load miniconda3`/`module load biopython` in the same `beforeScript` block have almost certainly had the *same* propagation failure all along, just invisibly masked because miniconda3 happens to reload anyway by default and the merge scripts that used biopython apparently didn't hard-fail without it.

**Why it matters**: Any future `module load X` (or bare `export PATH=...`, or `conda activate ...`) added to a process's `beforeScript` will silently fail to reach the actual command on this cluster unless `X` happens to already be part of the login shell's default baseline — with no error at the `beforeScript` stage itself, only a downstream "command not found" (or worse, a silently-wrong tool version if a *different* copy happens to be on the default `PATH`). This is exactly the kind of two-shell-boundary crossing that's invisible from reading the Nextflow config alone.

**Resolution**: Reverted the `beforeScript` module-load approach entirely. Put `module load duckdb 2>/dev/null || true` **inline, as the first line of each `script:`/`stub:` block**, immediately before the `duckdb -c ...` invocation — i.e. in the *same* shell as the command that needs it, with no shell boundary in between. This mirrors the already-existing, already-working `export PATH="${projectDir}/bin:\$PATH"` line several `MERGE_*`/`RUN_*` modules already have inline in their own `script:` blocks (e.g. `MERGE_CAZY`, `MERGE_MEROPS`) — that pattern was already the correct one; `beforeScript` module loads were the anomaly. Verified: after moving to inline `module load duckdb`, the same `-stub-run` that failed 27/30 processes on `duckdb: command not found` succeeded past that point (remaining failure was an unrelated `params.scripts` path issue caused by launching from the wrong directory, not this bug).

**How to apply — standing rule for this repo's Nextflow modules**: **Never put `module load`/`conda activate`/raw `export PATH=` in a process's `beforeScript` and assume it reaches `script:`/`stub:`.** `beforeScript` and the process payload run in different login-shell instances on this cluster, and login-shell startup here resets to a fixed module baseline, discarding anything loaded in between. If a process needs a tool not already in that baseline (check via `bash -l -c 'echo $LOADEDMODULES'`), load it **inline at the top of that same `script:`/`stub:` string**, right before first use — not via `beforeScript`, not via `process.beforeScript` in a shared label config. `beforeScript` is still fine for things that only need to be true *before* staging/setup (rare in this repo), just not for anything the command itself needs on `PATH`.

**Tags**: nextflow, beforeScript, module-load, login-shell, bash-l, hpc, slurm, environment-modules, duckdb, T-014, gotcha

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped for `duckdb` across all 15 `MERGE_*` modules touched by #27. Repo convention candidate: promote the "inline module load, never `beforeScript`, for anything a `script:`/`stub:` block needs on `PATH`" rule into CLAUDE.md's Script Conventions alongside the existing Nextflow `$`-interpolation rule — same class of "looks obviously correct, fails silently at runtime" Nextflow/HPC gotcha.

## The `nextflow run ...` *driver* process needs its own tools on `PATH` — a shell-launched pipeline silently fails if you forget to `module load` before invoking it (2026-08-01)

**Category**: gotcha

**What happened**: While doing a real (non-stub) verification run of `--pipeline compare_ani` (T-014 follow-up: testing BUSCO + ANI against one real assembly), launched `nextflow run ...` from an interactive shell that hadn't loaded the `singularity` module. Every `SKANI_COMPARE` task failed with `env: 'singularity': No such file or directory`, `exit 127`. This is a *different* failure class from the `beforeScript` gotcha logged above (2026-07-31): that one was about a process's own `script:`/`stub:` block losing an environment change across a login-shell boundary *inside* Nextflow's task-launch machinery. This one is simpler and more upstream — `singularity.enabled = true` (`conf/profile_ANI.config`) means the Nextflow *driver* JVM itself needs `singularity` on `PATH` to construct the wrapped container-exec command in the first place, and the driver only has whatever `PATH` the shell it was launched from had at that moment.

**Why it matters**: There's no config-file fix for this — it's not a `beforeScript`/inline-script placement question, it's "did the shell you typed `nextflow run` into have the right modules loaded before you hit enter." Easy to forget when switching between pipelines that need different toolchains (`duckdb` for merges, `singularity` for containerized ANI/functional steps, etc.), and the failure mode (every task in the run dies immediately with `exit 127`/`command not found`) looks superficially like the `beforeScript` bug even though the fix is completely different (load the module in your own shell, not in any Nextflow config).

**Resolution**: `module load singularity` in the launching shell before `nextflow run -profile ani ...`. Re-ran cleanly once done.

**How to apply**: Before launching any real (non-`-stub-run`) pipeline invocation by hand, check what the target profile's `singularity.enabled`/`module load X` assumptions are (grep the relevant `conf/profile_*.config`) and make sure *your own launching shell* has those tools loaded — this is separate from, and in addition to, whatever `beforeScript`/inline-script module loads the processes themselves need. A `command not found`/`exit 127` failing every single task in a run (not just one) is the signature of this class of problem, not a per-process bug.

**Tags**: nextflow, singularity, module-load, driver-process, launch-shell, ani, compare_ani, T-014, gotcha

**mitigation_type**: operational

**structural_mitigation_candidate**: Not shipped as a structural fix (this is inherently a launching-shell habit, not something a config file can enforce) — could be mitigated by wrapping common real-run invocations in a small launcher script (`bin/run_ani.sh` etc.) that does the `module load` itself, matching the pattern `scripts/build_BFD_duckDB.sh` already uses (`module load duckdb` at the top of the script, not left to the caller's shell).

## `-profile ani`'s singularity bind-mount paths must pre-exist on the host — no setup step creates them (2026-08-01)

**Category**: gotcha

**What happened**: Same real-run verification session as above. After fixing the `singularity` module-load gap, `SKANI_COMPARE` tasks still failed: `WARNING: skipping mount of .../work/ANI/pip_cache: stat ...: no such file or directory`, `FATAL: container creation failed: mount ...`. `conf/profile_ANI.config` sets `singularity.runOptions = "--bind ...,${launchDir}/work/ANI/pip_cache:/root/.cache/pip,${launchDir}/work/ANI/python_packages:/tmp/python_packages"` — singularity's `--bind` requires the *host-side* source path to already exist; neither directory is created by any setup process, Nextflow directive, or `mkdir -p` anywhere in the ANI pipeline. It only "worked" before because some prior run had already created them once as a side effect (probably manually, or because a previous run's work dir happened to get these paths made by something else first) and nobody hit a genuinely fresh `work/ANI/` until this session.

**Why it matters**: Any first-time ANI run against a clean `work/ANI/` directory (new checkout, new launchDir, or after `rm -rf work/ANI`) hits this immediately with a somewhat opaque singularity mount error, not an obvious "missing directory" message — easy to mistake for a real container/environment problem rather than a two-`mkdir -p`-away fix.

**Resolution**: `mkdir -p work/ANI/pip_cache work/ANI/python_packages` before the run. Re-ran cleanly (`[SUCCESS] completed=25 failed=0`) once created.

**How to apply**: Before a from-scratch ANI run (fresh `launchDir`/`work/ANI`), pre-create `work/ANI/pip_cache` and `work/ANI/python_packages` (or wherever `conf/profile_ANI.config`'s `singularity.runOptions` currently points — check it hasn't drifted). A structural fix (an explicit `mkdir -p` step, e.g. folded into a lightweight setup process or the launcher script) would remove this trap entirely; flagged as a candidate below rather than applied here since it's outside T-014's scope.

**Tags**: nextflow, singularity, bind-mount, ani, compare_ani, pip_cache, work-dir, first-run, gotcha

**mitigation_type**: operational

**structural_mitigation_candidate**: Not yet shipped. Repo convention candidate: add an explicit `mkdir -p ${launchDir}/work/ANI/pip_cache ${launchDir}/work/ANI/python_packages` guard (e.g. as a cheap `beforeScript` on the affected label — safe here since it's a plain `mkdir`, not a `PATH`-mutating module load, so the earlier `beforeScript` gotcha doesn't apply) so a from-scratch ANI run doesn't hit this opaque singularity mount failure.

## `stub:` blocks drift from the real schema silently — only a downstream consumer with real column-name assertions catches it (2026-08-01)

**Category**: gotcha

**What happened**: While implementing #28 (`build_BFD_duckDB.sh` → Parquet), built a real-data test fixture set to validate the rewritten SQL end-to-end. `MERGE_IDP`'s `stub:` block (`printf 'protein_id,idp_status,disordered_residues,total_residues\n'` / `printf 'protein_id,idp_status\n'`) turned out to never have matched the real `idp.csv.gz`/`idp_summary.csv.gz` schema at all — real production output (checked against `results/function/aiupred/*/*.idp.csv.gz`) is `species_prefix,protein_id,IDP_start,IDP_end,IDP_length,mean_score` for `idp` and `species_prefix,protein_id,IDP_residues,IDP_fraction,length` for `idp_summary`. Every sibling module's stub block (`MERGE_CAZY`, `MERGE_MEROPS`, `MERGE_PREDGPI`, `MERGE_SIGNALP`, `MERGE_TARGETP`, `MERGE_TMHMM`, `MERGE_WOLFPSORT`) correctly includes `species_prefix`; only `MERGE_IDP`'s was wrong, and predates the #27 Parquet-conversion work entirely (the stub `printf` content itself was never touched by #27, only wrapped in a `duckdb COPY` step).

**Why it matters**: `-stub-run` succeeded every time regardless of this gap, because nothing in the stub-run wiring check actually asserts column names — it just confirms every process produces *a* file at the expected path. The bug was invisible until something downstream (here, `build_BFD_duckDB.sh`'s `CREATE INDEX ... ON idp_summary(species_prefix)`) tried to use a column the stub never had. This is the same root cause class as the pfam `hmm_id`/`hmm_acc` vs real `pfam_id`/`pfam_acc` mismatch already fixed in `nextflow/tests/validate_outputs.py` during #27 — stub content and real-schema content are two independently-maintained sources of truth with nothing enforcing they agree, and `-stub-run` alone cannot detect drift between them.

**Resolution**: Fixed `MERGE_IDP`'s stub block to match the real schema (see #28's commit). Found by deliberately building a fully real-data test set (ran the real per-genome merge, not stub placeholders) rather than trusting the stub fixtures — the stub-only fixtures from #27's own testing had this exact gap sitting invisible in `nextflow/tests/output/tables/idp_summary.parquet` (2 columns, no `species_prefix`) the whole time.

**How to apply**: `-stub-run` verifies wiring (does the DAG connect, does every process produce *a* file), not schema correctness. Before trusting a stub fixture's columns for anything downstream (a `CREATE INDEX`, a `validate_outputs.py`-style assertion, a joined view), diff it against one real per-genome output file for that type. When adding or touching a `stub:` block, treat its column list as needing the same scrutiny as the real `script:` block's — it's not "just a placeholder," it's the thing every future `-stub-run` will silently trust.

**Tags**: nextflow, stub-run, schema-drift, merge-idp, idp-summary, T-014, validate-outputs, gotcha

**mitigation_type**: structural

**structural_mitigation_candidate**: Fixed for `MERGE_IDP` in #28. Repo convention candidate: when writing or reviewing a `stub:` block that emits a CSV/Parquet header, require it be copy-pasted or diffed from one real per-genome output file of that type, not hand-typed from memory of what the columns "should" be — this is now the second instance of exactly this drift (pfam in #27, idp in #28) caught only by accident (a downstream consumer with real assertions), not by `-stub-run` itself.

## PFAM/HMMER real-timing A/B results (T-015): scratch-copy doesn't help, MPI needs `--cpu-bind=none` and peaks around 4 tasks (2026-08-01)

**Category**: insight

**What happened**: Ran real (non-stub) `hmmsearch` A/B timing tests on this cluster (`highclock` partition, real Pfam-A DB, real `Malassezia brasiliensis` protein set, 3786 proteins) to answer T-015's open questions (`todo/pfam_hmmer_performance.md`) with data instead of guessing. Full setup/scripts in `analysis/pfam_hmmsearch_perf/`.

**Scratch-copy the Pfam-A DB — no benefit.** Copying the ~3.7GB DB from `/srv/projects/db/pfam/...` to node-local NVMe scratch took only 2.7s (fast NFS→NVMe transfer), and `hmmsearch` wall time was statistically identical whether reading from shared storage or the scratch copy: shared=5:00.35, scratch=5:02.54, shared-repeat=5:02.84 (within noise, ~1% spread). This task is compute-bound (HMM scan against a real protein set), not I/O-bound — the DB read itself is a negligible fraction of total wall time. **Decision: not worth implementing scratch-copy staging for PFAM.**

**MPI needs `--cpu-bind=none` on this cluster — otherwise every multi-task `srun` fails outright.** First attempt at MPI timing (`srun -N 1 -n 4 --mpi hmmsearch ...`, matching `RUN_PFAM`'s exact current invocation pattern) failed immediately: `srun: error: CPU binding outside of job step allocation ... Unable to satisfy cpu bind request`. This happened identically regardless of `--mpi=pmix`/`--mpi=pmi2`/omitted — ruled out via a minimal compiled MPI hello-world (`mpicc`, `MPI_Comm_rank`/`MPI_Comm_size`) that reproduced the exact same failure with plain `srun -n 4`, then succeeded with real distinct ranks reporting in once `--cpu-bind=none` (or `--cpu-bind=threads`) was added. Root cause: SLURM's default `--cpu-bind` mode conflicts with this node's cpuset/allocation mask for multi-task job steps — not an MPI-plugin selection problem. **This means `nextflow/modules/BFD/PFAM/main.nf`'s existing MPI code path is currently broken as written** (`srun -N ${pfam_nodes} -n ${pfam_tasks} --mpi hmmsearch`, no `--cpu-bind` override) — dormant/never-hit today only because `pfam_tasks`/`pfam_nodes` both default to 1 in `conf/profile_BFD.config`, so the MPI branch is never actually exercised in production.

**MPI task scaling is real but non-monotonic — peaks around `-n 4`, gets worse beyond that.** With the `--cpu-bind=none` fix applied, real `hmmsearch --mpi` timing vs. the single-task 4-thread baseline (5:00.35):

| `pfam_tasks` (`-n`) | real workers (rank 0 is a dispatcher, not a worker) | wall time | vs. baseline |
|---|---|---|---|
| 1 (current default, non-MPI, 4 threads) | — | 5:00.35 | baseline |
| 2 | 1 | 9:15.56 | **-85% (worse)** — coordination overhead, nothing to actually parallelize against |
| 4 | 3 | 4:16.58 | **+15% faster** — best result |
| 8 | 7 | 4:49.62 | +4% faster |
| 16 | 15 | 5:11.42 | -3% (worse than baseline) |

Domain-count sanity check (`grep -vc '^#' *.tblout`) confirmed identical output (7819 domains) across every configuration — the scaling differences are pure wall-time, not partial/incorrect results.

**Why it matters**: The relationship is not "more MPI tasks = faster" — for this genome's protein-set size, useful parallel work runs out well before 16 workers, and beyond ~4 tasks the master/worker coordination and scheduling overhead outweighs the added parallelism, eventually erasing the entire gain. `-n 2` specifically is a trap: HMMER's `--mpi` mode dedicates rank 0 to dispatching work rather than searching, so `-n 2` gives exactly 1 real worker plus coordination overhead — worse than no MPI at all. Any future decision to enable `pfam_tasks > 1` in production needs the `--cpu-bind=none` fix and should default to a small task count (this data suggests ~4) rather than "more is better," and should be re-profiled against a larger, more typical protein set (3786 proteins is on the small side) before committing to a specific default.

**Resolution**: T-015 (`todo/pfam_hmmer_performance.md`) updated with these results and a concrete recommendation: don't implement scratch-copy; if enabling `pfam_tasks > 1` is ever pursued, it needs the `--cpu-bind=none` fix in `RUN_PFAM`'s `mpi_launch` construction first, and should start from `pfam_tasks≈4` rather than assuming higher is better. Findings also written into the personal `nextflow-hpcc` skill (`~/.claude/skills/nextflow-hpcc/SKILL.md`) as two new general-purpose subsections (`beforeScript` module-load gotcha review + the new `--cpu-bind=none` requirement) so any future Nextflow-on-this-HPCC session gets this automatically, not just this repo.

**How to apply**: Before assuming an I/O-bound bottleneck (shared storage, network mount) explains a slow HPC task, profile the actual wall-time split — this task turned out to be purely compute-bound despite reading a multi-GB file from network storage per run. Before assuming "more MPI tasks always helps," profile a real scaling curve on real data — the non-monotonic peak-then-decline shape here would have been missed by testing only one or two task counts. Before trusting any MPI timing number, verify the *launch mechanism itself* isn't silently broken (the `--cpu-bind` failure here would have made every "does MPI help" test either crash outright or need investigation) using a cheap, fast-to-iterate compiled sanity check rather than debugging through the slower real tool.

**Tags**: performance, pfam, hmmer, hmmsearch, mpi, cpu-bind, slurm, srun, scratch-partition, task-scaling, T-015, insight

## `hmmsearch --cpu` and `--mpi` are mutually exclusive — a second real bug found only by testing the actual production command, not a reconstructed one (2026-08-01)

**Category**: gotcha

**What happened**: After applying the `--cpu-bind=none` fix to `RUN_PFAM`'s MPI code path (previous entry) and validating it with raw `sbatch` A/B scripts, went back and applied the fix directly to `nextflow/modules/BFD/PFAM/main.nf`, then validated through a real Nextflow-launched run (`--pipeline BFD --asmid ... pfam_tasks=4`) against a real, more typical genome (*Aspergillus fumigatus* B-1-43-1, 10,079 proteins — larger than the 3786-protein genome used for the original A/B timing, per that analysis's own caveat that it should be re-checked on a bigger genome). The real Nextflow-generated command failed immediately: `hmmsearch --mpi --cut_ga --noali --cpu 4 ...` → `Failed to parse command line: Option --cpu is incompatible with option(s) --mpi`. This combination was never actually exercised by the earlier raw `sbatch` A/B tests, which manually constructed `hmmsearch --mpi ...` commands *without* `--cpu` — only the real module's `script:` block unconditionally appends `--cpu ${task.cpus}` regardless of MPI mode, so the earlier validated fix was incomplete for the actual code, not just missing this one flag interaction.

**Why it matters**: The A/B timing analysis (previous entry) was real and correct data, but it validated a hand-built reproduction of the MPI invocation, not the literal command the module would actually generate. The gap between "the mechanism I tested" and "the mechanism the code actually runs" hid a second, code-specific bug that only surfaced when testing through the real pipeline with real config wiring. Confirmed directly (not assumed): `hmmsearch --mpi --cpu 4 --cut_ga --noali /dev/null /dev/null` reproduces the exact same "incompatible" error standalone.

**Resolution**: `RUN_PFAM`'s `script:` block now only appends `--cpu ${task.cpus}` in the non-MPI (single-task, threaded) branch; the MPI branch passes no `--cpu` flag at all (`def cpu_flag = params.pfam_tasks > 1 ? "" : "--cpu ${task.cpus}"`). Re-validated end-to-end through the real Nextflow pipeline against the larger genome: baseline (`pfam_tasks=1`) 9m 9s, fixed MPI (`pfam_tasks=4`, both bugs fixed) 7m 40s — a real 16.2% improvement, consistent with (slightly better than) the ~15% seen on the smaller test genome in the original A/B analysis. Domain output validated non-empty and structurally correct (16,874 domain hits).

**How to apply**: When a real-tool timing/validation result depends on a hand-constructed reproduction of a command (raw `sbatch` script, manual shell invocation) rather than the actual generated command from the real pipeline code, don't treat that result as validating the code path itself — it only validates the *mechanism* (here: does `--cpu-bind=none` fix MPI launches at all). Before merging or trusting a fix, run it through the real code with real config wiring at least once, ideally against data more representative of production scale than whatever was convenient for the first quick test.

**Tags**: pfam, hmmer, hmmsearch, mpi, cpu-flag, incompatible-options, RUN_PFAM, T-015, gotcha

**mitigation_type**: structural

**structural_mitigation_candidate**: Shipped in `nextflow/modules/BFD/PFAM/main.nf` (branch `pfam-mpi-cpu-bind-fix`). Both fixes (`--cpu-bind=none` + conditional `--cpu` flag) now validated through the real Nextflow-launched pipeline against a real, production-representative genome, not just a hand-built reproduction.
