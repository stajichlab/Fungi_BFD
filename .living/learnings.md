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

## `GENOME_CLEAN_BATCH` passed the samples.csv PHYLUM *name* string straight to `AAFTF fcs_gx_purge -t`, which requires a numeric NCBI taxid (2026-08-03)

**Category**: bug

**What happened**: `nextflow/modules/funannotate/genome/GENOME_CLEAN_BATCH/main.nf` had a "prefer the curated PHYLUM column from samples.csv over a live taxonkit lookup" fast path (added to dodge stale/incomplete taxdump lookups for brand-new taxids). It read `samples.csv` column 7 (`PHYLUM`) via `awk` and passed the result directly to `fcs_gx_purge -t "$phylum"`. That column holds a taxonomic **name** (e.g. `Ascomycota`), confirmed by inspecting `samples.csv` directly — but `AAFTF fcs_gx_purge -t/--taxid` requires a **numeric** NCBI taxid (its own `--help` text: "NCBI Taxonomy ID for contamination matches, i.e. 4890 for Ascomycota"). The sibling non-batch module, `GENOME_CLEAN/main.nf`, never had this bug — it only ever resolves phylum via `taxonkit lineage | reformat -f "{p}" | name2taxid`, which already lands on a numeric taxid, and it validates the result against `^[0-9]+$` before use. The batch variant's fast path bypassed that conversion+validation entirely.

**Why it matters**: `taxonid` (NCBI_TAXONID) and the phylum's own taxid are different numbers — one identifies the organism, the other identifies its phylum. AAFTF's `-t` flag wants the latter, resolved to a number. Passing a name string like `Ascomycota` would either error out inside `fcs_gx_purge` or silently be misinterpreted, and there was no format check on this fast path to catch it (unlike the numeric-regex check already present in `GENOME_CLEAN/main.nf`).

**Resolution**: `GENOME_CLEAN_BATCH/main.nf` now takes the samples.csv PHYLUM value as `phylum_name`, resolves it through `taxonkit name2taxid` to get the numeric taxid, and falls back to the existing taxonid→lineage resolution if that comes back empty. Added the same `^[0-9]+$` validation-and-warn as `GENOME_CLEAN/main.nf` has, so a failed resolution is visible in logs rather than silently mis-feeding `fcs_gx_purge`. Verified the fix against real data: `samples.csv` row for `GCF_010015735.1_Aaoar1` has `PHYLUM=Ascomycota`; `echo Ascomycota | taxonkit name2taxid` → `4890`, matching AAFTF's own documented example.

**How to apply**: Any time a pipeline stage reads a taxonomy column out of `samples.csv` (PHYLUM, CLASS, ORDER, etc.) and feeds it to a tool expecting a **numeric NCBI taxid** rather than a name, run it through `taxonkit name2taxid` first — don't assume a "curated" spreadsheet column is already in the format the downstream tool wants. Check both/either of the two `GENOME_CLEAN*` modules if adding a third variant, since they can drift out of sync (this bug existed only in the batch variant, not the original).

**Tags**: funannotate, GENOME_CLEAN, GENOME_CLEAN_BATCH, taxonkit, taxid, phylum, fcs_gx_purge, AAFTF, samples.csv, bug

**mitigation_type**: structural

**structural_mitigation_candidate**: Fixed in `nextflow/modules/funannotate/genome/GENOME_CLEAN_BATCH/main.nf` — PHYLUM name is now resolved to a numeric taxid via `taxonkit name2taxid` with fallback to lineage-based resolution and numeric-format validation, matching `GENOME_CLEAN/main.nf`'s existing behavior.

## `file()` with a glob pattern triggers a Nextflow deprecation warning; current code uses `files()` for globs (2026-08-03)

**Category**: gotcha

**What happened**: User reported seeing a Nextflow warning about using `file*()` (i.e., `file()` with a glob/asterisk) instead of `files()`. A full search of the current `nextflow/` tree found no `file("...*...")` calls with glob patterns; the only globbing helper, `gatedGlobIn()` in `nextflow/modules/common/utils.nf`, already uses `files("${baseDir}/${glob}")`. A clean `-profile test -stub-run --pipeline BFD` completed without emitting the warning, and `.nextflow.log` contained no `file()`/\`files()` deprecation message.

**Why it matters**: The warning is real in Nextflow 26.x — `file()` returns a single path and warns when the pattern matches multiple files; `files()` is the correct factory for glob multi-file collections. If it appears in future runs, it means some code path is passing a glob to `file()` (possibly via a variable that contains `*`), not that the current telomere modules are at fault.

**Resolution**: No code change was required for the current tree. Confirmed the telomere merge (`MERGE_TELOMERES`) and finder (`FIND_TELOMERES`) modules do not use `file()` with globs.

**How to apply**: When you see this warning, grep for `file(` calls whose argument contains `*`, `?`, or `**` (including dynamically built strings). Replace those with `files()` or, if a single path is required, ensure the pattern is unambiguous. If the warning is observed on code that currently looks clean, capture the exact `.nextflow.log` snippet — Nextflow sometimes reports the originating script line.

**Tags**: nextflow, file-vs-files, glob, deprecation-warning, telomeres, gotcha

## `find_telomeres.py`'s 3'-end search was structurally broken from initial implementation — reversed the sequence string without complementing it, then searched with a reverse-complemented pattern (2026-08-04)

**Category**: bug

**What happened**: The initial telomere finder (`nextflow/bin/find_telomeres.py`, `regex_terminal_hits`/`fuzzy_terminal_hits`) searched the 3' scaffold end by doing `search_seq = seq[::-1][:search_window]` (plain string reversal, no complementing) and then matching that against `search_monomer = reverse_complement_pattern(canonical)` (reverse **and** complement). Those two transforms don't correspond — reversing a sequence without complementing it, then searching for the reverse-complement of a motif, essentially never matches a real terminal repeat. Verified directly: `reverse("CCCTAACCCTAA...")` = `"AATCCCAATCCC..."`, which does not contain `"CCCTAA"` (the pattern the code was actually searching for at that end).

A second, independent bug compounded this: `terminal = (start_coord == 0)` required the repeat to start at literally index 0 of the scaffold, with zero tolerance. Real assemblies commonly carry a few non-canonical bases before the repeat locks into clean register (e.g. a chromosome starting `AACCCTAAACCCTAA...` — the clean register starts at index 2). With `--allow-internal` defaulting False, these genuinely-terminal telomeres were silently dropped.

**Why it wasn't caught at the time**: `.living/findings/fungal-telomeres.md`'s original validation pass only checked 2 tracts on one genome (*N. crassa* OR74A), both of which happened to be 5'-end hits — the systematic 100%-5'-only skew wasn't noticed as a red flag, and the tool was declared "ready for kingdom-wide application."

**Impact measured on real data**: on `GCA_900073075.1_CML3066.v2` (*F. graminearum*, a genuine 4-chromosome telomere-to-telomere reference where all 8 chromosome ends visibly carry telomere repeats in the raw FASTA), the buggy tool reported only 1 of 8 real telomeric ends (~12% recall) under production parameters. Across 5 test genomes total, 100% of reported hits were `5prime` — zero `3prime` hits ever, which in hindsight is itself diagnostic (real chromosome-level assemblies should show roughly balanced 5'/3' hit counts).

**Resolution**: Fixed in `nextflow/bin/find_telomeres.py` — both ends are now searched in the scaffold's native forward orientation (`search_monomer` already encodes the correct pattern per end; no sequence reversal needed at all). `terminal` now uses a `--terminal-tolerance` (default 10bp) instead of an exact 0/scaffold_len match, and a numeric `distance_to_end` column was added to the TSV output so downstream consumers can apply their own cutoff rather than trusting a single boolean. A related third bug (regex-path repeat counting used the canonical monomer instead of the actually-matched `search_monomer` for 3'-end tracts, undercounting to 1) was fixed in the same pass, and the `monomer` output column was unified to always report the canonical (5'->3') form regardless of `--fuzzy`/regex mode (previously the two modes disagreed on what that column meant). Re-tested on the same 5 genomes: *F. graminearum* CML3066 now reports 7/8 real telomeric ends (the 8th genuinely has no repeat in the raw sequence, confirmed by direct inspection — correct negative, not a miss); *N. crassa* OR74A went from 1→13 hits; *A. nidulans* ASM4168428v1 went from 1→9 hits; both now show a healthy 5'/3' balance.

**How to apply**: When writing/reviewing any "search the reverse-complement strand" logic, do not conflate "reverse the string" with "reverse-complement the string" — they require correspondingly different patterns to search with, and a mismatch between the two transforms produces code that runs without error but silently never matches. Also: a **100%-skewed 5'-vs-3' (or +/- strand) hit ratio in a symmetric biological search is a red flag worth investigating before declaring a tool validated**, even if the few spot-checked hits look individually correct. Before trusting a boundary/edge-matching condition (`== 0`, `== len(seq)`), check it against real assembly data — real sequences are rarely in perfectly clean register right at an edge; a `==` requirement where a tolerance/`<=` would do is a common source of false negatives at exactly the biologically interesting spot.

**Tags**: telomeres, find_telomeres, reverse-complement, coordinate-bug, terminal-tolerance, bioinformatics-review, bug

## The c01-c30 "Abu Dhabi" no-AVX2 SLURM exclusion was applied to `FUNANNOTATE_TRAIN`/`FUNANNOTATE_PREDICT` but missed `RNASEQ_PREPARE`, which runs the same crash-prone code path (2026-08-07)

**Category**: gotcha

**What happened**: A prior fix (see `FUNANNOTATE_TRAIN`/`FUNANNOTATE_PREDICT` entries in `nextflow/conf/profile_funannotate.config`) added `--exclude=c[01-30]` to SLURM `clusterOptions` because c01-c30 are old AMD "Abu Dhabi" (Opteron 6300 / Piledriver) CPUs with no AVX2, and funannotate's compiled binaries SIGILL-crash on them instantly. `RNASEQ_PREPARE` (`nextflow/modules/funannotate/rnaseq/RNASEQ_PREPARE/main.nf`) runs `funannotate train -i ... --stop_after_trinity`, which internally calls `hisat2-build` for the genome-guided Trinity assembly — the exact same crash-prone code path — but its `clusterOptions` in `conf/profile_funannotate.config` never got the `--exclude=c[01-30]` added, and it also runs on the `preempt` account/queue (the same queue pool that includes the abu_dhabi nodes). Result: `hisat2-build` failures inside `funannotate train` that looked like a recurrence of the "abu dhabi processor problem."

**Why it matters**: When a SLURM node-exclusion fix is scoped to specific process names because a symptom shows up there first, any other process that shells out to the same underlying crash-prone binary (hisat2, augustus, bam2hints, etc.) on the same queue pool needs the identical exclusion — the fix does not automatically propagate across process blocks in `withName:` config, even within the same pipeline.

**Resolution**: Added `--exclude=c[01-30]` to both attempt-1 and later-attempt `clusterOptions` for `RNASEQ_PREPARE` in `nextflow/conf/profile_funannotate.config`, matching the existing `FUNANNOTATE_TRAIN`/`FUNANNOTATE_PREDICT` pattern.

**How to apply**: When you add a SLURM node exclusion (or any node-targeting fix) to fix a SIGILL/illegal-instruction/random-crash pattern for one process, grep `nextflow/modules/**/main.nf` for other processes that shell out to the same binary or wrap the same external tool (`funannotate train`, `hisat2`, `augustus`, etc.) and check whether they run on the same queue pool — apply the exclusion everywhere the crash path can be reached, not just where the symptom was first observed. Considered alternatives: a SLURM `--constraint` feature filter (no `avx2` feature is registered on this cluster, only `abu_dhabi`/`amd` — a constraint would have to be phrased as `--constraint='!abu_dhabi'`, which is equivalent to but less explicit than `--exclude=c[01-30]`) and building funannotate's dependencies from source without AVX2 optimization (rejected — throws away the perf benefit on the many non-abu_dhabi nodes just to support ~30 old nodes cluster-wide).

**Tags**: nextflow, slurm, abu-dhabi, avx2, hisat2, funannotate-train, rnaseq-prepare, sigill, gotcha

### [2026-08-09] `clean_genome_fa.py` did not drop funannotate's "alphabet < 4" contigs; predict aborted on 3-base contigs

**Category**: bug

**What happened**: funannotate predict failed on `Rugonectria_rugulosa_M7C8` (`GCA_058970095.1_ASM5897009v1`) with "Found 1 bad contigs, where alphabet is less than 4 [this should not happen]" then `ERROR: funannotate predict did not produce expected GBK`. The offending contig `JBUFDE010000307.1` (2,156 bp) contained only A(1310)/G(268)/T(578) — 3 distinct nucleotides, zero C. `funannotate.library.analyzeAssembly()` counts the **number of distinct characters** in each contig (uppercased, including N/ambiguity) and flags any contig with `< 4` distinct as suspect, aborting predict unless `--force`. `clean_genome_fa.py` only filtered by `--len`, so bad-alphabet contigs sailed through into the masked genome and predict blew up downstream.

**Why it matters**: A 3-base or 2-base contig is a strong contamination/QA signal (a real chromosomal contig always has all 4 bases), and funannotate will hard-stall the whole predict on it. The clean step is the right single choke point to catch these; catching them post-hoc in the run `work/` dir wastes a full predict attempt.

**Resolution**: `nextflow/bin/clean_genome_fa.py` now mirrors funannotate's exact check — counts distinct uppercase chars per contig, and by default (`--min-alphabet 4`, matches funannotate) **drops** any contig with fewer than 4 distinct characters, emitting one stderr WARNING per dropped contig with its composition (e.g. `A:1,310, G:268, T:578`) plus a summary count. Opt out with `--min-alphabet 1`. No Nextflow changes needed (GENOME_CLEAN/GENOME_CLEAN_BATCH call it with only `--len`). The affected genome's clean/masked inputs must be regenerated (the `.masked.fasta.gz` storeDir artifact embeds the bad contig) before predict can proceed.

**How to apply**: When "preprocessing" a genome before funannotate, run the *same* alphabet census funannotate runs (`len(set(seq.upper())) < 4`) rather than only checking IUPAC membership/length — "looks like DNA" is not the same as "has a complete alphabet". Validate a clean script against the exact library.py predicate it's supposed to defend against, ideally on the real failing assembly.

**Tags**: funannotate, clean_genome_fa, alphabet, contig-filter, predict-failure, GBK, bug, bioinformatics-review

### [2026-08-10] AAFTF container: `beforeScript`-free, in-script `singularity exec` wins here

**Category**: insight

**What happened**: Switched CALC_ASM_STATS, GENOME_CLEAN, GENOME_CLEAN_BATCH from `module load AAFTF` / `module load taxonkit` to the new `AAFTF-v0.7.0-beta.2.sif` (built by the AAFTF repo's release CI from its AAFTF.def; installed in the Nextflow singularity cacheDir; image hardcodes `AAFTF_DB=/opt/aaaftf_db` and ships taxonkit v0.20.0). These processes were NOT moved to Nextflow `process.container`; instead each script block defines `SING_BINDS` + `SING="singularity exec ${SING_BINDS} ${AAFTF_SIF}"` and calls `$SING AAFTF ...` / `$SING taxonkit ...`. The AAFTF repo's CI only builds the SIF/pushes releases ($ nothing to do on the BFD side), and verified empirically: `$SING AAFTF --version` → `AAFTF 0.7.0b2`, `$SING taxonkit version` → `v0.20.0`, `bash -lc` login shell sources `/etc/profile.d/aaaftf.sh` so `AAFTF` resolves + `AAFTF_DB` points at `/opt/aaaftf_db`, and the bind exposes host `/srv/projects/db/AAFTF_DB` contents.

**Why it matters**: These genome-clean processes also need host-side tools (`pigz`, R/conda, `clean_genome_fa.py`) and the FCS-GX DB staged into the HOST `/dev/shm` by `setup_fcs_shm.sh` (`/dev/shm/gxdb/all`). `process.container` would containerize the whole task (losing host /dev/shm + host modules without heavy binds) — manual `singularity exec` scopes the container to just the AAFTF/taxonkit calls while the rest of the script stays native. The container is invoked with `--bind ${taxondb}:${taxondb},/srv/projects/db/AAFTF_DB:/opt/aaaftf_db,/dev/shm:/dev/shm`.

**Resolution / How to apply**: All three modules drop their `module load AAFTF|taxonkit`, add `module load singularity` + the SING/SING_BINDS block, and prefix every in-container call (`AAFTT fcs_gx_purge`, `AAFTF assess`, `taxonkit lineage/reformat/name2taxid`) with `$SING`. Root param `params.aaftf_sif` added in `nextflow/config`. If the image is rebuilt, update only `params.aaaftf_sif` (and `_BATCH`'s copy) — no module-level code changes.

**Tags**: aaftf, singularity, container, taxonkit, fcs-gx, dev-shm, geneclean, nextflow, insight

### [2026-08-10] AAFTF-SIF genome-clean real-run validation: workflow passes end-to-end, cluster contention + cleanup quirks bite

**Category**: validation / gotcha

**What happened**: Ran the first real `GENOME_CLEAN_BATCH` under the AAFTF-SIF wiring (3 small real genomes — S. cerevisiae, S. pombe, C. neoformans — via `conf/test_clean_real.config` with `only_clean=true`, `run_ani_reuse=false`, `clean_batch_size=3`). Driver launched interactively through a detached `screen` (not sbatch; see below). Pipeline completed `Succeeded: 2` (SETUP_TAXONDB reuse + the batch): the 465 GB FCS-GX db staged into host `/dev/shm` at 520 MB/s (498.6 GB in 914 s, logged to `logs/nextflow/fcs_gx_shm_timing.tsv`), `clean_batch_1` ran on h04 (job 27332111, 18:05 runtime), and `/dev/shm` was freed on teardown. Outputs: `input_clean_genomes/<asmid>.fa.gz` + `clean/*.purge.fasta.gz` + `*.purge.fcs_gx-taxonomy.tsv.gz` with correct divisions (S288C `fung:ascomycetes`/`budding yeasts` agg-cov 0.998; C. neoformans `fung:basidiomycetes` agg-cov 0.979). `TO_ADD_TO_SUPRESS.csv` empty (nothing below min length).

**Why it matters / gotchas**:
1. **highmem "Resources" wait**: the batch asks for 16 cpus + 500 GB pinned to `-w h04,h05,h06`. All three nodes had >500 GB *already allocated* by other jobs (~655/643/800 GB), leaving <280 GB free each — SLURM therefore waited ~6 h before h04 freed enough. On a busy cluster the 500 GB request against fixed node pins is the scheduling bottleneck, not anything in the pipeline.
2. **`cleanup = true` hides per-batch outputs**: `profile_funannotate` sets `nextflow.cleanup=true`, so the task work dir (including `clean_batch_*.manifest.tsv`) is scrubbed on success — the manifest is not on disk after a clean run. Harmless for `only_clean` (no downstream consumer), but don't look for it in `input_clean_genomes/` after a successful run.
3. **sbatch launcher still broken for real Nextflow runs** (unchanged 2026-08-10): launching `nextflow run` from an sbatch job fails deterministically on `.nextflow/history.lock` (Java `NoSuchFileException`) and top-level `mkdir` (Permission denied) even though interactive/`screen` shells work identically from the same dir — the Java/Nextflow process under the job environment is the difference. Detached `screen` + `> logs/nextflow/real_clean_driver.log 2>&1` is the proven pattern for real runs (production driver runs this way too).
4. **`run_ani_reuse` must be false for `only_clean` test configs**: `workflows/funannotate.nf` hard-errors at build time if `run_ani_reuse=true` while `params.abinitio_reuse_csv` is unset, killing the run before any task starts.

**How to apply**: For real clean-step validation/acceptance, launch via `screen -dmS` from `nextflow/` (not sbatch); expect an unbounded wait if the pinned highmem nodes are contended; read timing from `logs/nextflow/fcs_gx_shm_timing.tsv` and validate taxonomy by `pigz -dc <asmid>.purge.fcs_gx-taxonomy.tsv.gz` (header-only lines mean no flagged contamination rows). Keep `conf/test_clean_real.config` + `tests/real_clean/` together as the reproducible unit (see `nextflow/tests/real_clean/REAL_CLEAN.md`).

**Tags**: aaftf, singularity, fcs-gx, dev-shm, geneclean, batch-clean, highmem, slurm, nextflow, validation, gotcha

### [2026-08-10] HMMER SIF has no MPI runtime — containerize only the non-MPI hmmsearch path

**Category**: insight

**What happened**: Containerized `RUN_PFAM`'s hmmsearch (non-MPI default) with `hmmer_3.4--hdbdd923_1.sif` via `singularity exec --bind $PFAM_DB:$PFAM_DB`, mirroring the AAFTF SIF pattern (`params.pfam_sif` in `nextflow.config`). Smoke test of the exact invocation against `Testus_fungus_STRAIN1.proteins.fa` + the real `$PFAM_DB` (2026-01-27-Pfam38.2) exited 0 in 12.8 s with real domtbl/tblout hits. The `pfam_tasks>1` MPI branch was left alone because **the SIF ships `hmmsearch`/`phmmer`/`hmmscan` but no `mpirun`** — it cannot drive HMMER's `--mpi` mode, so that branch still uses `module load hmmer/3.4-mpi` + `srun ... --mpi`. In the MPI branch `$SING` is unset (empty), so the invocation stays byte-identical to the pre-container module-load form.

**Why it matters**: Biocontainers' hmmer image is built without the MPI stack, so "containerize all of hmmer" would silently break `--mpi` (or require a different MPI-enabled image and SIMD/PMI bindings that don't carry over to SLURM srun cleanly). The correct scope is: container for the threaded non-MPI path, host module for the MPI path. Also re-confirmed the repo's Groovy-GString rule: shell vars (`SING`, `SING_BINDS`, `PFAM_SIF`, `$PFAM_DB`) must be written `\${SING}`/`\$PFAM_DB` inside the `"""..."""` script block, while Groovy defs (`${mpi_launch}`, `${params.pfam_sif}`) are plain `$`/`${}`.

**How to apply**: When containerizing an HPC tool that has both threaded and MPI modes, check the image for `mpirun`/MPI libs first (`singularity exec <sif> which mpirun`) and only containerize the mode the production config actually uses (here `pfam_tasks=1`). Gate any mode that needs MPI behind the host module, and keep the branch divergence explicit in the same `if/else` so the container hit is one isolated hunk.

**Tags**: hmmer, hmmsearch, pfam, mpi, mpirun, singularity, container, nextflow, insight
### [2026-08-10] SwissProt search + MEROPS containerization: engine choice, DB compat, and MEROPS_DB pin

**Category**: validation / gotcha

**What happened**: Built the BFD SwissProt homology module (`RUN_SWISSPROT`) with two switchable engines, both from Singularity containers (`params.swissprot_search` = `diamond` default | `blastp`):
- **diamond 2.2.5** (`quay.io/biocontainers/diamond:2.2.5--he361c42_0`, staged as `.../shared/lib/singularity_cache/diamond_2.2.5--he361c42_0.sif`). Chosen over the host `diamond` module (2.1.24) and over `getwilds/diamond:latest` because the getwilds image is stale (latest 2.1.16, 2025-12-30) — user preference: "prefer to use diamond from container over local if recent", and quay biocontainers is the canonical source.
- **NCBI blastp 2.16.0** (pre-existing `.../shared/lib/singularity_cache/depot.galaxyproject.org-singularity-blast-2.16.0--h66d330f_5.img`, verified `blastp 2.16.0+`). The same image now drives the **containerized RUN_MEROPS** (replaces `module load ncbi-blast`).

**Why it matters / gotchas**:
1. **diamond 2.2.5 reads the 2.1.24-built dmnd with no rebuild** — verified against `nextflow/lib/swissprot/uniprot_sprot.dmnd` (build version 178, 575,503 seqs). Confirmed again in the full real-proteome run.
2. **`db-merops/124` exports `MEROPS_DB=/srv/projects/db/MEROPS/124`, not /120** — I initially assumed /120; the modulefile at `/opt/linux/rocky/8.x/x86_64/modules/db-merops/124` sets /124 (the /120 symlink/default is older). `merops_scan.lib` is a symlink to `meropsscan.lib` with full makeblastdb 6-file DB present in /124.
3. **The `swissprot` label needs its own `withLabel` block** in `conf/profile_BFD.config` (8 cpus/16 GB/8 h/`epyc`) and in `conf/test.config` (local + `beforeScript=':'`, or the stub-run would try to queue SLURM and load unrelated modules).
4. **Eighteen-column custom blasttab is engine-agnostic**: both diamond `-f 6` and blastp `-outfmt "6 ..."` emit `qseqid sseqid pident positive nident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp qlen slen stitle`, so `MERGE_SWISSPROT` has one code path. Diamond 2.2.5 emits all 18 columns (verified: 1,167 HSPs from a 100-protein real subset, all width 18).
5. **Raw 80-80 is a derived flag, not a pre-filter**: cutoffs stay permissive (`-e 1e-5 -k/max_targets 20`); `func_transfer_80_80` = `pident>=80 AND query_cov>=0.8 AND hit_cov>=0.8` is computed in `merge_swissprot.py`. `query_cov=length/qlen`, `hit_cov=length/slen`.
6. **Real-proteome benchmark**: `Moesziomyces_antarcticus_T-34` proteome (6,288 seqs) → diamond `--sensitive -k 10 -e 1e-5` @ 8 threads = **2m19s** (4,609 aligned), backing the 8-cpu/8-h swissprot label.

**How to apply**: To refresh the SwissProt release, replace `nextflow/lib/swissprot/uniprot_sprot.{fasta,dat}.gz` from `https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/`, rebuild `uniprot_sprot.{dmnd,blastDB}` (makeblastdb 2.16.0 container), and bump `reldate.txt`; `BUILD_SWISSPROT_ANNOT` re-parses the flatfile each launch (no storeDir) so the annot table auto-propagates. Switch engines per-run with `--swissprot_search blastp`.

**Tags**: swissprot, diamond, blastp, merops, singularity, container, nextflow, gotcha, benchmark
### [2026-08-10] Real-SLURM smoke of containerized MEROPS + SwissProt: /scratch is node-local, and BUILD_DUCKDB wants a full table set

**Category**: validation / gotcha

**What happened**: Ran a real-SLURM smoke (single real genome, Moesziomyces_antarcticus_T-34 / FDB2C4F4) covering containerized `RUN_MEROPS` + `RUN_SWISSPROT` (diamond) + `MERGE_*` + `BUILD_SWISSPROT_ANNOT`. Two real gotchas surfaced:
1. **`/scratch/jstajich` is node-local NVMe, not shared across SLURM compute nodes** (stat shows `/dev/nvme0n1p6` XFS). Nextflow stages `path(...)` inputs as symlinks; a `samples.csv` under `/scratch/...` was invisible to the compute node the job landed on → `FileNotFoundError: 'samples.csv'` inside `MERGE_SAMPLES` (looked like a code bug, was environmental). Fix: keep ALL smoke inputs/outputs/workdir on the shared filesystem (`/bigdata`). A workdir inside `/scratch` is equally fatal.
2. **`BUILD_DUCKDB` (`run_build_duckdb`) fails on an incomplete table dir**: `build_BFD_duckDB.sh` hard-globs the full parquet set (`asm_stats.parquet`, ...), so a SwissProt/MEROPS-only smoke errors "IO Error: No files found that match the pattern .../asm_stats.parquet". Expected — that step is for complete runs; smoke sets `run_build_duckdb=false`.

**Other notes from the smoke**: (a) keep a real run's workdir/logbook out of the repo work dir by launching from a dedicated run dir and passing `-resume`, `-w`, and let `launchDir` defaults produce the production layout (`results/`, `tables/`, `db/`, `work/`) — no outdir/tables/db overrides needed; (b) with `run_setup=false` the workflow still hard-validates the input dirs at startup (workflows/BFD.nf `dirIndex(...)`), so point `pep_dir/cds_dir/gff_dir/genome_dir/trna_dir` at the real `_runs/input/{pep,cds,gff3,dna,trna}` (10,943 genomes each; naming `{id}.proteins.fa`, `{id}.cds-transcripts.fa`, `{id}.gff3`, `{id}.scaffolds.fa`, `{id}.trna.gff3`); (c) `setsid bash driver.sh </dev/null >/dev/null 2>&1 &` (not bare `nohup ... &`) so a foreground tool timeout cannot SIGTERM the process group and kill the run.

**Outcome**: smoke fully green — `merops.parquet` (13 cols), `swissprot.parquet` (68,238 rows), `swissprot_annot.parquet` (575,503 rows = full release), `samples/species.parquet`; all 68,238 swissprot rows join to an annotation row on `swissprot_acc = accession`.

**How to apply**: For pipeline smokes on UCR HPCC, keep every path on `/bigdata`, launch from a dedicated run dir with `-resume -w`, and only enable `run_build_duckdb` when the full table set is present.

**Tags**: nextflow, slurm, hpcc, scratch, node-local, smoke, swissprot, merops, duckdb, gotcha
### [2026-08-10] Blastp-engine smoke benchmark + storeDir cache gotcha: switching swissprot_search needs a fresh storeDir

**Category**: validation / gotcha / benchmark

**What happened**: Ran the full smoke with `swissprot_search='blastp'` (NCBI blast 2.16.0 container) end-to-end under real SLURM for the same T-34 genome. It completed green (Succeeded 6, Duration 59m 9s). Key numbers vs the diamond path (same proteome, both 8 threads, `-max_target_seqs`/`-k 20`, `-e 1e-5`):
- **diamond 2.2.5 `--sensitive`**: 2m19s, `swissprot.parquet` = 68,238 rows.
- **blastp 2.16.0 `-seg yes -soft_masking true`**: ~56 min (RUN_BLAST only), `swissprot.parquet` = 80,502 rows (55,295 distinct accessions).
- Both emit an identical 18-column blasttab (verified in storeDir output), both merge/annot paths green, both 100% `swissprot_acc`→`accession` join resolution. Blastp returns MORE hits than diamond here (soft-masking rescues alignments diamond's seeding misses). Both accessions resolve in the annot table. Differences are engine-inherent, not a merge bug.

**The operational gotcha (why this took two launches)**: trying to switch engines via `-resume` silently does NOT re-run `RUN_SWISSPROT` — the process has a `storeDir` (`${params.outdir}/swissprot/...`), and when the expected output (`${locustag}.blasttab.gz`) already exists there, Nextflow marks the task `skipped/stored` and never looks at the changed `swissprot_search` param; the cached downstream `MERGE_*` replay diamond output into `tables/` too. To genuinely re-execute, launch a FRESH run dir (empty storeDir) or clear/redirect outdir.

**How to apply**: `params.swissprot_search` only takes effect if the `outdir`/`storeDir` for that species is empty or the species is new. For engine A/B comparisons, use separate run/outdir trees (as here: `smoke_swissprot` diamond vs `smoke_blastp` blastp). In production, the two engines are interchangeable per-run via the same strip; diamond stays the default (56 min vs 2m19s is a 24× runtime difference — irrelevant for a single genome but huge across 10,943).

**Tags**: swissprot, blastp, diamond, benchmark, storedir, nextflow, cache, gotcha, smoke

**Tags**: swissprot, blastp, diamond, benchmark, storedir, nextflow, cache, gotcha, smoke
### [2026-08-12] Smoke-tested ghcr.io/nextgenusfs/funannotate rust container (beta.2/beta.3) for Trinity/PASA/EVM: 3 packaging bugs found, GeneMark confirmed absent, ~6.9x rust-Trinity speedup

**Category**: validation / gotcha / benchmark

**What happened**: Real `funannotate train --stop_after_trinity` smoke test (real genome+RNA-seq: `Alternaria_alstroemeriae`, GCA_053542345.1, 34 Mb / 70 contigs, real normalized PE reads) run under `singularity exec` against `ghcr.io/nextgenusfs/funannotate:v1.9.0-beta.2` and `:1.9.0-beta.3`. Confirmed via `funannotate check --show-versions`: both images have EVidenceModeler (Rust), PASA (Rust), and Trinity Rust utilities all `ENABLED`, and GeneMark is genuinely absent (`GENEMARK_PATH not set`) — matches the maintainer's design intent (GeneMark license forbids redistribution; predict must stay on a host module or a future standalone task).

Found 3 real packaging bugs while getting Trinity to actually run under `singularity exec` (as opposed to `docker run`, which sources more of the image's setup):
1. **`/venv/bin/fasta` was a dangling symlink** to nonexistent `/venv/bin/fasta3` (bioconda package `fasta3` actually installs a binary named `fasta36`, not `fasta3`) — `funannotate train`'s upfront `CheckDependencies` hard-blocks on this before Trinity even starts, even under `--stop_after_trinity`. Fixed in maintainer's Dockerfile (`ln -sf /venv/bin/fasta36 /venv/bin/fasta`); present in beta.3's binary (`fasta: 36.3.8g` in `check --show-versions`) but not yet in a rebuilt image at time of testing.
2. **`LD_LIBRARY_PATH=/venv/lib` was never set** as a Docker `ENV` (confirmed: `/.singularity.d/env/10-docker2singularity.sh` only auto-imports actual Dockerfile `ENV` lines, not conda `activate.d` hooks) — broke Trinity's `bamsifter` (`_sift_bam_max_cov: error while loading shared libraries: libdeflate.so.0`) and would affect any `/venv/lib`-linked compiled tool under `singularity exec`/`run`, not just `docker run` (which happens to source more of the env). Fixed in maintainer's Dockerfile; not yet in a built/pushed image as of beta.3.
3. **`which` is broken for every invocation inside the container**, even zero-argument-adjacent calls — Debian's `bash-completion`-provided `which()` shell function (present even under non-interactive `bash -c`, no login/BASH_ENV gating found) always forwards GNU-long-option flags (`--tty-only --read-alias ...`) to the real binary behind `/etc/alternatives/which`, but that binary is the minimal `debianutils` variant which only understands `-a`/`-s` and errors `Illegal option --` on any invocation. funannotate's own Python dependency checks are unaffected (uses `shutil.which()`, confirmed still resolves correctly), but real Perl code in the image shells out to `which` (`grep` found it in PASA's `Launch_PASA_pipeline.pl` and `seqclean`), so this is a live risk for the PASA stage specifically, not just cosmetic. Recommended fix: install the `which` apt package (GNU-compatible) in the runtime stage so the alternative resolves to a binary that understands the flags the wrapper function passes — cleaner than trying to neutralize a bash-completion-provided function.

`/opt/databases` (FUNANNOTATE_DB) doesn't need to pre-exist in the image — Singularity 3.9.3 auto-creates the bind mountpoint (`--bind host_db_dir:/opt/databases` tested working against this project's real `funannotate_db`). `AUGUSTUS_CONFIG_PATH=/venv/config` is baked read-only into the SIF though (`touch` inside it fails: `Read-only file system`) — needs either a writable bind mounted over `/venv/config`, or (preferred, matches this project's existing host-module pattern of `export AUGUSTUS_CONFIG_PATH=${params.augustus_config}`) an env-var override to a separately-seeded writable dir (`cp -r /venv/config/. <dir>/` once, then `--env AUGUSTUS_CONFIG_PATH=<dir>`).

With workarounds 1+2 applied (`fasta`-symlink shim + `SINGULARITYENV_LD_LIBRARY_PATH=/venv/lib`), both beta.2 and beta.3 completed a full real genome-guided Trinity assembly end-to-end (19,018 / 19,019 transcripts respectively, essentially identical — expected, same input). Timing: beta.2 3h23m00s, beta.3 3h19m00s (8 CPUs, 10,532 clusters). Compared against a real module-based (`funannotate/dev-1.9`, confirmed non-rust: `NOTE: Rust EVidenceModeler not found`, `NOTE: Rust PASA not found`) Trinity-GG run preserved in `rnaseq_data/Lepraria_finkii.funannotate-trinity.log` (14,600 clusters, 16 CPUs, 16h03m47s wall): normalizing to CPU-seconds/cluster gives ~63.4 (module/perl) vs ~9.0-9.2 (container/rust) — **~6.9x faster per cluster-CPU-second with rust Trinity**. Caveat: different species/genome (lichen vs. *Alternaria*), not a controlled same-input A/B, so treat as indicative not a precise multiplier.

Separately discovered: `analysis/funannotate_train_stage_timing/`'s existing 401-row module-baseline dataset has **zero** rows with any Trinity-GG time captured (`trinity_related_seconds` is 0 in every row, `trinity_gg` never the dominant stage) — consistent with F-004 (Trinity-GG never runs in 336/336 PASA-dominant sampled `FUNANNOTATE_TRAIN` logs, because per-strain runs reuse a shared pre-built Trinity). The real module-based Trinity-GG timing lives in `rnaseq_data/*.funannotate-trinity.log` (171 preserved logs from `RNASEQ_PREPARE`, the representative-strain-only process), a different log family the existing stage-timing script never profiled.

**How to apply**: For any future container-vs-module Trinity timing comparison, pull baselines from `rnaseq_data/*.funannotate-trinity.log` (`RNASEQ_PREPARE`), not `genome_annotation{,_training}/*/logfiles/funannotate-train.log` (`FUNANNOTATE_TRAIN` — PASA-only reuse, no fresh Trinity). For any container smoke test invoked via `singularity exec` (not `docker run`), explicitly set `LD_LIBRARY_PATH=/venv/lib` and verify `which`-dependent tools aren't silently failing — `docker run`'s extra environment sourcing can mask both bugs that `singularity exec` exposes. The container's own internal `logfiles/funannotate-train.log` was far sparser (2 timestamped lines total for the whole Trinity run) than the console stream; for stage-level timing capture, redirect and parse stdout (`[Mon DD HH:MM AM/PM]:` format, minute resolution) or the per-tool `training/Trinity-gg.log` (second resolution, but only covers the initial clustering/setup phase — the ~10K per-partition assembly jobs are dispatched without their own timestamps).

**Tags**: funannotate, singularity, container, rust, trinity, pasa, evm, genemark, packaging-bug, benchmark, gotcha, dockerfile

**Tags**: funannotate, augustus, augustus-config-path, gotcha, predict, evm
### [2026-08-12] funannotate predict requires AUGUSTUS_CONFIG_PATH's basename to be literally "config", or crashes with UnboundLocalError deep into the run

**Category**: gotcha

**What happened**: While building `analysis/genemark_es_contribution/` (GeneMark-ES A/B rerun of `funannotate predict`), pointed `--AUGUSTUS_CONFIG_PATH` at a private per-job writable copy named `<predictdir>.augustus_config` (to avoid writing into the live shared `Fungi_BFD_runs/lib/augustus/3.5/config` while real production jobs were concurrently using it). All 6 jobs ran for 9-25 real minutes — genome inflate, Hisat2/PASA reuse, GeneMark-ES self-training all completed successfully (`11,117 predictions from GeneMark`) — then crashed at the Augustus-training step:
```
File ".../funannotate/predict.py", line 2265, in main
    AUGUSTUS_BASE,
UnboundLocalError: local variable 'AUGUSTUS_BASE' referenced before assignment
```
Root cause (`predict.py` around line 507): `AUGUSTUS_BASE` is only assigned inside `if os.path.basename(os.path.normcase(os.path.abspath(AUGUSTUS_CONFIG_PATH))) == "config":` — an exact string-equality check on the basename, not a check for expected config contents. Any `--AUGUSTUS_CONFIG_PATH` whose leaf directory isn't literally named `config` silently skips the assignment, and the variable is unconditionally referenced ~1800 lines later once Augustus training actually starts (so the failure surfaces late, after GeneMark/PASA-filtering work is wasted, not at argument-parsing time).

**How to apply**: Any custom/private `AUGUSTUS_CONFIG_PATH` (container smoke tests, A/B analyses, anything that seeds its own writable copy instead of using the shared one) MUST have `config` as the literal final path component — e.g. `<private_dir>/config`, not `<private_dir>.augustus_config`. Verify with `[ "$(basename "$AUGUSTUS_CONFIG_PATH")" = config ]` before a long run, since the crash only surfaces after other expensive stages (GeneMark-ES, PASA filtering) have already completed.

**Tags**: funannotate, augustus, augustus-config-path, gotcha, predict, evm

**Tags**: funannotate, augustus, augustus-config-path, gotcha, predict, evm
### [2026-08-12] Real end-to-end validation of GENEMARK_RUN: all 3 tests pass, --genemark_gtf confirmed to bypass predict's internal GeneMark call

**Category**: validation

**What happened**: After Fable-reviewing and implementing the standalone `GENEMARK_RUN` Nextflow process (nextflow/docs/GENEMARK_RUN_DESIGN.md), ran real (non-stub) end-to-end tests against `Penicillium_citrinum_NRRL_1841` (same genome used in the F-008 GeneMark-ES contribution A/B test), via a throwaway standalone test workflow that includes the real module file directly (`analysis/genemark_run_validation/`):
1. **Fresh `gmes_petap.pl --ES`**: 11,116 gene models — matches almost exactly the A/B test baseline's `GeneMark: 11112` from predict's own internal call. ~15.5 min at 8 cores.
2. **Fast `gmes_petap.pl --predict_with <mod>`** (reusing a real `.mod` already on disk from the A/B test): 100,966-line GTF (near-identical total to fresh-ES's 101,000 lines, as expected — same genome/model, different code path). Completed in **~7 min at only 2 cores**, vs fresh-ES's ~15.5 min at 8 cores — real confirmation `--predict_with` genuinely skips training, not just seeds a faster convergence.
3. **Real `funannotate predict --genemark_gtf <gtf>`**: resolves the design doc's open question directly — predict's log shows `GeneMark path: .../genemarkESET/4.72_lic`, `GeneMark appears to be functional? True` (GENEMARK_PATH was available), but `RunGeneMarkES()` was never invoked; `--genemark_gtf` took priority and short-circuited it, confirmed by the "Summary of gene models" line showing `GeneMark 1 11116` — the exact count from step 1's GTF, proving the externally-supplied file was what actually got used. Final gene count: 11,202 vs the A/B baseline's 11,198 — 4-gene difference, normal EVM/tbl2asn tie-breaking noise, not a degradation.

**Test-harness-only bug found (not a production bug)**: Nextflow's CLI parser treats `--flag ''` (explicit empty string) as a bare boolean flag — `--shared_mod ''` silently became Groovy `Boolean true`, not an empty string, corrupting the rendered `.command.sh`. Confirmed this can't affect the real pipeline (FUNANNOTATE_PREDICTION.nf passes these as plain Groovy string literals inside `tuple()` calls, never through CLI arg parsing). Fixed in the test harness by omitting the flag entirely rather than ever passing `--param ''`.

**How to apply**: Any future standalone Nextflow test script exposing module inputs via `--param` CLI flags: never pass an intentionally-empty value as `--flag ''` — omit the flag and let the script's own `params.x = ''` default apply, or use a non-empty sentinel and translate it inside the script. This is specific to Nextflow's own CLI parser, not something `bash`-level quoting can work around.

**Tags**: funannotate, genemark, evm, predict, container-migration, validation, gmes_petap, predict_with, nextflow-cli-gotcha

**Tags**: funannotate, genemark, genemark-et, evm, predict, ab-initio
### [2026-08-12] Evaluated GeneMark ET mode (T-022): corrected the hints-source assumption, then hit a real branch-point-region failure on short fungal introns

**Category**: gotcha / validation (negative result)

**What happened**: Evaluating T-022 (GeneMark ET mode for `GENEMARK_RUN`), the original design doc scoping wrongly assumed `RunGeneMarkET()`'s `b2h`-tagged intron hints come from `funannotate_train.coordSorted.bam` (raw RNA-seq read alignment) via the external AUGUSTUS `bam2hints` binary. Tracing the actual filter (`predict.py`/`library.py`: `"\tintron\t" in line and "\tb2h\t" in line`) showed this is wrong — `bam2hints`'s *default* `--source` is `E`, which never matches. The real `b2h` tag comes from **transcript-alignment** evidence: funannotate's own `bam2ExonsHints()` (`library.py:1952`, a Python BAM→hints converter distinct from the AUGUSTUS binary) explicitly sets `btag = "b2h"` when run on minimap2-aligned transcript evidence, and AUGUSTUS's `blat2hints.pl` (BLAT path) uses the same convention. The correct input, `training/transcript.alignments.bam`, is already produced and retained by `FUNANNOTATE_TRAIN` — no new alignment work needed, confirmed on disk.

Verified the corrected understanding empirically: `bam2hints --intronsonly --source=b2h --in=training/transcript.alignments.bam --out=hints.gff` on a real genome (`Penicillium_citrinum_NRRL_1841`) produced 21,400 hints, every line independently confirmed to match `RunGeneMarkET()`'s filter exactly.

But the actual `gmes_petap.pl --ET <hints> --sequence genome.fa --fungus` run on those hints **failed**: `error, hash is empty: bp_seq_select.pl`, preceded by `warning, no data in specified range: histogram.pl`. Root-caused by measuring the hints' own intron-length distribution: `min=32 median=65 mean=81.6 max=2934` bp. GeneMark's branch-point-region extraction searches a 40-50bp window *within* each intron (`--bp_region_length 50`/`--min_bp_region_length 40` defaults) — for a median-65bp fungal intron, that window consumes nearly the entire intron, leaving `bp_seq_select.pl` with too little sequence to build a branch-point model. Confirmed this branch-point step is shared code between `--ET` and `--fungus`'s internal `ES_C` sub-mode (`gmes_petap.pl` lines 483/539/878 all call `bp_seq_select.pl`) — and the earlier **successful** `--ES --fungus` run (11,116 genes, same genome) went through an analogous step fine, because ES's candidate introns come from its own broad genome-wide self-discovered search, not hints externally restricted to aligned-transcript regions only.

**How to apply**: Before trusting any funannotate/GeneMark source-tracing from grepping tool names (e.g. "bam2hints" appears in the code but isn't the actual data source used) — trace the actual filter condition/tag values, not just which binaries are invoked nearby. For GeneMark-ET specifically on fungal genomes: don't assume default parameters transfer from GeneMark's typical plant/animal tuning — short fungal introns (median ~65bp here) may need `--bp_region_length`/`--min_bp_region_length` tuned down before `--ET` can run at all, not just for quality. Untested whether a smaller window still captures real branch-point signal.

**Tags**: funannotate, genemark, genemark-et, gmes_petap, bam2hints, fungal-introns, branch-point, negative-result

**Tags**: funannotate, genemark, genemark-et, braker, evm, predict, container-migration
### [2026-08-12] Fixed GeneMark-ET (T-022): missing strand assignment, root-caused via BRAKER's real pipeline, verified with a real gmes_petap.pl --ET run

**Category**: gotcha / validation (positive result, follow-up to the earlier negative-result entry same day)

**What happened**: Following up on the earlier "hash is empty: bp_seq_select.pl" failure (previous learnings entry, same date) — user asked to confirm whether the ET hints-derivation logic was real funannotate code or de-novo, and suggested checking BRAKER (`Gaius-Augustus/BRAKER`, the more battle-tested Augustus+GeneMark-ET+RNA-seq pipeline, also GeneMark-licensed and also unable to containerize) as a reference implementation. Fetched `braker.pl`, `get_gc_content.py`, and `filterIntronsFindStrand.pl` directly from BRAKER's GitHub `master` branch (real upstream code, not paraphrased).

Found the actual fix: BRAKER never feeds raw `bam2hints` output straight to GeneMark. Every intron hint first goes through `filterIntronsFindStrand.pl` (Artistic License, Simone Lange & Katharina Hoff): checks the genome sequence at each intron boundary against canonical splice-site dinucleotides (GT-AG/GC-AG/AT-AC), assigns the correct strand, and **silently drops any intron without a canonical splice site**. The earlier failed attempt's hints had `.` (unstranded) in the strand column — branch-point signal is inherently strand-specific, so GeneMark's `bp_seq_select.pl` had nothing orientable to work with. The `--bp_region_length`/short-intron-length hypothesis from the earlier entry was a plausible-sounding but wrong lead — BRAKER uses the same GeneMark defaults and never touches those flags either; `--gc_donor` (also investigated) turned out to be a fixed default (0.001), not genome-adaptive, ruling that out too.

Corrected recipe, verified with a real `gmes_petap.pl --ET` run (`analysis/genemark_run_validation/et_eval2/`, matching BRAKER's `get_genemark_hints()` exactly, `braker.pl` lines 4888-4995): `bam2hints --intronsonly` → `filterIntronsFindStrand.pl --score` → sort (4-key, matching BRAKER's exact multi-`sort` pipe) → `join_mult_hints.pl` (AUGUSTUS script, already on PATH via the funannotate module) → `gmes_petap.pl --ET`. 20,603 of 21,400 hints (96%) passed stranding/canonicalization for `Penicillium_citrinum_NRRL_1841`. **The `--ET` run completed successfully: 10,776 gene models** (vs. the ES run's 11,116 — same order of magnitude). Vendored `filterIntronsFindStrand.pl` into `nextflow/bin/vendor/` (not modified from upstream, license header preserved) since it isn't bundled with this project's funannotate/Augustus install.

Also confirmed while investigating: funannotate's protein-alignment hints (`exonerate2hints()`, used for `--protein_evidence`/SwissProt) tag with `src=XNT`, not the `src=P` BRAKER's `--EP`/`--ETP` (protein-informed GeneMark) modes expect — same category of tag-convention mismatch as the intron-hints case, so combining protein evidence into GeneMark isn't a free drop-in; noted as a future consideration, not attempted.

**Separately found and fixed a latent bug while re-verifying the wired `FUNANNOTATE_PREDICT` module**: `predict.py:567` unconditionally zeroes GeneMark's EVM weight (`StartWeights["genemark"] = 0`) whenever `gmes_petap.pl` isn't found on the *host running predict* — with no check for whether `--genemark_gtf` was supplied as an alternative. Not triggered by the earlier Test 3 validation (that host had `gmes_petap.pl` present), but would silently zero out a correctly-supplied `--genemark_gtf`'s contribution the moment `FUNANNOTATE_PREDICT` actually moves onto the container (no `gmes_petap.pl` there at all) — the precomputed evidence would be consumed by predict but then discarded by EVM with zero weight, no error, no warning. Fixed: `FUNANNOTATE_PREDICT/main.nf` now passes `-w genemark:1` explicitly whenever `genemark_gtf` is non-empty. Caught a second near-bug while fixing the first: funannotate's `-w`/`--weights` argparse option has no `action='append'` (confirmed directly against argparse), so a second `-w` flag on the command line silently *replaces* the first rather than merging — an initial fix that appended a second `-w genemark:1` would have silently dropped the existing `-w codingquarry:0 glimmerhmm:0`. Both weight overrides now go into a single merged `-w` argument group.

**How to apply**: When porting a hints-generation pipeline between two related-but-independently-evolved tools (funannotate vs. BRAKER, both wrapping the same underlying `gmes_petap.pl`), don't assume tag conventions or pipeline steps transfer — `b2h` vs `E` vs `P` source tags differ between the two, and BRAKER's `filterIntronsFindStrand.pl` canonicalization step has no funannotate equivalent at all. When any CLI tool's argument parser doesn't document `action='append'` behavior, verify empirically (a 3-line Python `argparse` script) before assuming repeated flags merge — assuming append when the real behavior is replace is an easy, silent way to lose an earlier override.

**Tags**: funannotate, genemark, genemark-et, braker, filterintronsfindstrand, evm-weights, argparse-gotcha, container-migration

**Tags**: nextflow, workflow-projectdir, genemark, genemark-et, standalone-test-harness
### [2026-08-12] GeneMark-ET wired into GENEMARK_RUN + FUNANNOTATE_PREDICTION, real end-to-end validated (T-022 closed); workflow.projectDir resolves per-launched-script, not to a fixed repo root

**Category**: implementation / gotcha (positive result, follow-up to the two earlier same-day GeneMark-ET entries)

**What happened**: Following the validated ET recipe (previous learnings entry), wired it for real: `GENEMARK_RUN/main.nf` gained a `mode`/`training_bam`-gated branch (`bam2hints --intronsonly` → vendored `filterIntronsFindStrand.pl` → `join_mult_hints.pl` → `gmes_petap.pl --ET`, falling back to `--ES` self-training when `training_bam` is empty -- e.g. a genome with no RNA-seq, confirmed via `FUNANNOTATE_RNASEQ.nf:249`'s `predict_no_rnaseq` branch that such genomes still reach `FUNANNOTATE_PREDICTION.nf`). Added `trainingTranscriptBamFor(out)` to `utils.nf` (mirrors `sharedGenemarkModFor()`'s existence-check shape) and threaded it + the (previously declared-but-unreferenced, confirmed via grep) `genemark_mode` param through all 3 `GENEMARK_RUN`/`GENEMARK_RUN_SIB` call sites in `FUNANNOTATE_PREDICTION.nf`.

Verified with a real `-stub-run` through the actual entry point (`nextflow run nextflow/main.nf ...`) first (`completed=4 failed=0` with the new 10-element tuple), then real-smoke-tested the actual wired module (not the by-hand recipe) via a standalone throwaway script.

**Real gotcha hit while building that smoke test, not the module itself**: `workflow.projectDir` (which `GENEMARK_RUN` uses to locate the vendored `filterIntronsFindStrand.pl` via `${workflow.projectDir}/bin/vendor/...`) resolves to *the originally-launched script's own containing directory* for the whole run -- not a fixed repo root, and not per-included-module-file. A first attempt placing the smoke-test script at `nextflow/tests/manual/genemark_run_smoke.nf` failed with `Can't open perl script ".../nextflow/tests/manual/bin/vendor/filterIntronsFindStrand.pl"` -- `workflow.projectDir` was `nextflow/tests/manual/`, not `nextflow/`. Fixed by moving the test script to `nextflow/genemark_run_smoke.nf`, sibling to `main.nf` -- confirmed correct once the pipeline banner started reading "BigFungiData 1.0" (it auto-picked-up the co-located `nextflow/nextflow.config`), and the real run succeeded: **10,780 gene models** (matches the standalone recipe's 10,776), plus a fast-reuse regression check (100,966-line GTF, exact match to the earlier reuse test) confirming the new 10-element `GENEMARK_RUN` input tuple didn't break the existing ES/`--predict_with` path.

**How to apply**: Any process/module that resolves a path via `${workflow.projectDir}` only works correctly when invoked (directly or via `include`) from a top-level script that lives at the intended project root -- for this repo, that's `nextflow/main.nf`, so `${workflow.projectDir}/bin/...` references are safe in production but will silently resolve wrong in any standalone test script placed elsewhere. Standalone module smoke tests that need this to work correctly must live at the same directory level as the real entry point, not in a nested test subdirectory.

**Tags**: nextflow, workflow-projectdir, genemark, genemark-et, standalone-test-harness, container-migration, T-022

### [2026-08-13] funannotate container beta.5: fasta36 symlink + LD_LIBRARY_PATH bugs fixed; `which` "bug" is a host env leak, not an image bug, and doesn't actually block PASA

**Category**: gotcha / validation (revised same day — initial root-cause was wrong)

**What happened**: Pulled `ghcr.io/nextgenusfs/funannotate:1.9.0-beta.5` (built to fix an unrelated bowtie2/AVX2 issue) and re-checked the 3 packaging bugs found in beta.2/beta.3. Two are now fixed: `funannotate check --show-versions` reports `fasta: 36.3.8g` cleanly (`/venv/bin/fasta36` resolves directly, no more fasta3/fasta36 symlink mismatch), and `LD_LIBRARY_PATH=/venv/lib:/.singularity.d/libs` is now set in-container (confirmed via `ldd $(command -v diamond)` showing no missing libs).

The third looked broken at first — `singularity exec $SIF bash -c 'which perl'` fails with `Illegal option --` / exit 2 — but tracing it further showed **this is not an image defect at all**. `env` inside the container shows `BASH_FUNC_which%%=...` alongside clearly host-only artifacts (`BASH_FUNC_module%%`, `BASH_FUNC_scl%%`, `BASH_FUNC__module_raw%%` — this HPC's environment-modules Tcl wrappers). The broken `which()` function is the **host's own bash function** (Rocky Linux's `which2`-style wrapper, written for GNU `which` v2.21's `--tty-only --read-alias --read-functions --show-tilde --show-dot` flags), leaking into the container via Singularity's default full environment passthrough. The container's actual `/usr/bin/which` is a minimal `debianutils` build that doesn't understand those flags — hence the crash, but only when something explicitly execs `bash`.

Confirmed this has **no practical runtime impact**: `BASH_FUNC_*%%`-style exported functions are a bash-only mechanism — `dash` (Debian's default `/bin/sh`, what Perl's `system()`/backticks spawn) never imports them. `singularity exec $SIF sh -c 'which perl'` → `/venv/bin/perl`, exit 0, clean. PASA's `.dbi` scripts shell out via `system()` → `sh`, not `bash`, so the codepaths that actually run at pipeline execution time are unaffected. Also reconfirmed `--cleanenv` alone makes `bash -c 'which perl'` succeed too, further isolating the host as the source.

Separately confirmed via the same `funannotate check` run: GeneMark is correctly absent with a graceful error (`ERROR: gmes_petap.pl not installed`, no crash), and all three rust-optimized engines report enabled (`EVidenceModeler (Rust): ENABLED`, `PASA (Rust): ENABLED`, `Trinity Rust utilities: ENABLED (4/4 on PATH)`).

**Why it matters**: No Dockerfile change is warranted for this — it would be fixing a problem in the wrong place (the image is fine; the host's bash function is what's incompatible with the image's `which` binary, and only bites interactive/explicit-`bash` use, not the `sh`-based paths PASA/Perl actually use). This also means the earlier concern that beta.4 was supposed to fix this and didn't is moot — there was nothing in the image to fix.

**Resolution**: No image or repo change needed. If it's ever annoying interactively, `singularity exec --cleanenv $SIF ...` or unsetting the function in the calling shell (`unset -f which`) before `singularity exec` clears it — invocation-side, not build-side.

**How to apply**: Before attributing a container behavior to "the Dockerfile," check whether it reproduces with `--cleanenv` and whether it depends on `bash` specifically vs `sh` — `BASH_FUNC_*%%` env-var-exported shell functions are a classic host-environment leak through `singularity exec`'s default passthrough (this project already documented that passthrough happens; this is the first time it produced a false-positive "bug"). Re-run the fasta36/LD_LIBRARY_PATH checks on future betas regardless — those two were genuine image fixes.

**Tags**: funannotate, container-migration, singularity, env-passthrough, packaging-bug, pasa, beta5, false-positive

### [2026-08-13] beta.5 real end-to-end confirm run: `run_sra_fetch=false` silently trains NOBODY, not just "skips fetching new data"

**Category**: gotcha (caught via a real full run, not stub-run)

**What happened**: First real (non-stub) confirmation run of the beta.5-containerized `FUNANNOTATE_TRAIN`/`GENEMARK_RUN` ES/ET wiring, on 4 real strains (2 species) in an isolated `do_container_confirm_test/` launch dir. Set `run_sra_fetch: false` intending only to skip hitting SRA/NCBI network, since both species' normalized reads already existed locally (`rnaseq_reads/Penicillium_citrinum_norm_*.fastq.gz` real data; `rnaseq_reads/Pichia_senei_norm_*.fastq.gz` genuine 0-byte "no data" markers). The run completed successfully (9/9, 0 failed, real `predict_results` with sane gene counts for all 4 strains) — but `genome_annotation_training/` ended up completely empty, and every strain's `predict_results/*.parameters.json` showed `"genemark": "selftraining ES"`, including the two Penicillium citrinum strains that genuinely have RNA-seq. GeneMark never ran ET for anyone.

Root cause: `FUNANNOTATE_RNASEQ.nf`'s own header comment says it plainly (line 9) — "With it \[run_sra_fetch\] off, samples pass straight through untrained." The entire `RNASEQ_PREPARE`/`FUNANNOTATE_TRAIN` pathway construction is inside `if (params.run_sra_fetch.toBoolean())`; when false, `predict_input_ch` is built from the no-training branch unconditionally, regardless of what's already sitting in `rnaseq_reads/`. `skip_sra_query=true` is a sub-flag *inside* that same `if` block (skips only the NCBI query step, reading cached CSVs instead) — it does nothing when `run_sra_fetch` itself is false, since the whole branch it lives in is never entered.

**Why it matters**: This is an easy trap for any future isolated/test run that reuses already-fetched species-level RNA-seq (the common case, since `SRA_FETCH` uses `storeDir` and is naturally idempotent/cheap to re-invoke once cached) — the intuitive "false = don't fetch, but still use what's there" reading is wrong; `run_sra_fetch=false` means "don't train," full stop. A silent, non-erroring wrong answer: the run doesn't fail, it just produces GeneMark-ES-only, Augustus-BUSCO-generic (not PASA-trained) predictions for everyone, which look superficially fine (real gene counts, no errors) unless you specifically check `genome_annotation_training/` or `parameters.json`'s recorded evidence sources.

**Resolution**: Set `run_sra_fetch: true` (with `skip_sra_query: true` to still avoid live NCBI queries) and re-ran. `SRA_FETCH`'s `storeDir "${launchDir}/rnaseq_reads"` recognizes the already-present `*_norm_*.fastq.gz` outputs and should reuse them rather than re-fetching/re-normalizing.

**How to apply**: To reuse already-fetched RNA-seq without hitting SRA, use `run_sra_fetch=true` + `skip_sra_query=true` together — never `run_sra_fetch=false` as a "just use what's cached" shortcut. When validating any run that's supposed to exercise the RNA-seq/training path, don't trust "completed successfully with real gene counts" alone — check `genome_annotation_training/*/training/` is actually populated and spot-check `predict_results/*.parameters.json`'s `genemark`/`augustus` `source` fields, since a fully-untrained fallback path produces plausible-looking output with no error.

**Tags**: funannotate, rnaseq, run_sra_fetch, container-migration, genemark, T-022, silent-misconfiguration, beta5

### [2026-08-14] beta.5 real end-to-end confirm run, take 2: TRAIN + GENEMARK_RUN (ES/ET/reuse) all validated; containerized `funannotate update` hits a real PASA hook-resolution bug

**Category**: validation (positive) + gotcha (new bug, unresolved)

**What happened**: After fixing the `run_sra_fetch` misconfiguration (previous entry) and clearing the stale `predict_results` it produced, a clean rerun (SLURM job 27462389, `do_container_confirm_test/`) genuinely exercised the full wired path: containerized `FUNANNOTATE_TRAIN` (confirmed real `transcript.alignments.bam` + `trinity-GG.fasta` for all 3 *Penicillium citrinum* strains), `GENEMARK_RUN`/`GENEMARK_RUN_SIB`, and `FUNANNOTATE_PREDICT`/`FUNANNOTATE_PREDICT_SIB`. Results: all 4 strains produced sane gene counts (B8014 10,380; NRRL_1841 9,899; NRRL_756 10,351; *Pichia senei* 4,650), and — the key confirmation — both sibling strains' `predict_results/*.parameters.json` show `"genemark": {"source": "ab-initio-reuse", ...}`, proving `GENEMARK_RUN_SIB`'s `--predict_with` shared-`.mod` reuse path genuinely fired against real ANI-driven representative/sibling assignments (B8014 representative, NRRL_1841/NRRL_756 siblings at 99.5% ANI) computed by a real `compare_ani` prerequisite pass, not a stub. (Direct ES-vs-ET proof for the representative/standalone strains isn't recoverable after the fact — `GENEMARK_RUN` has no `publishDir`/`storeDir` and `cleanup=true` wipes successful work dirs — but the branch is a deterministic function of `training_bam` existence, which is independently confirmed present for B8014 and absent for *Pichia senei*, so this is as strong a confirmation as the design allows without adding new instrumentation.)

Two real, unrelated setup/packaging issues surfaced in the same run:
1. **My own test-dir gap**: `FUNANNOTATE_ANNOTATE` failed — `lib/template.sbt` was never created when the isolated `do_container_confirm_test/lib/` was set up (only `augustus/`, `taxdump`, `swissprot_fungi.faa` were carried over). Fixed with a symlink to the shared production file. Not a pipeline bug.
2. **Real beta.5/PASA packaging bug, root-caused precisely**: `FUNANNOTATE_UPDATE` failed for good — PASA's `Load_Current_Gene_Annotations.dbi` (invoked during `funannotate update`'s annotation-comparison step; auto-selects the `GFF3::GFF3_annot_retriever::get_annot_retriever` hook whenever its `-P` arg is a `.gff3` file — unconditional, not something the calling conf can opt out of) dies with `Error, couldn't resolve path for GFF3::GFF3_annot_retriever` (`Pasa_conf.pm` line 149/`_get_hook`). Initial read (wrong): assumed the hook `.pm` needed to live under `PerlLib/GFF3/` and was simply missing from the rust-optimized PASA fork's `PerlLib/` tree. **Corrected after checking `~/projects/funannotate/PASApipeline` (the actual fork source, `git@github.com:hyphaltip/PASApipeline.git`, branch `rust_optimize`)**: PASA's hook modules are deliberately NOT in `PerlLib/` — they live in a separate top-level `SAMPLE_HOOKS/` dir (`SAMPLE_HOOKS/GFF3/GFF3_annot_retriever.pm`, `SAMPLE_HOOKS/GTF/Gtf_annot_retriever.pm`), loaded dynamically via a dedicated `HOOK_PERL_LIBS` config key that `Pasa_conf.pm::_load_hook_perl_libs()` reads and pushes onto `@INC` before hook resolution. This project's own `~/.pasa/pasa_conf/conf.txt` already has `HOOK_PERL_LIBS=__PASAHOME__/SAMPLE_HOOKS` set correctly (matches upstream's own `Docker/conf.txt` reference) — the config side is NOT the bug. **Actual root cause**: `PASApipeline/scripts/install.sh` copies `PerlLib`, `Launch_PASA_pipeline.pl`, `pasa_conf`, `schema`, `scripts`, `pasa-plugins`, and `misc_utilities` into the install prefix (lines ~160-195) but has **no line copying `SAMPLE_HOOKS/`** — confirmed `/venv/opt/pasa/src/SAMPLE_HOOKS` is entirely absent in the built beta.5 image, so `HOOK_PERL_LIBS` points at a directory that was never installed. This was the first-ever real (non-stub) test of containerized `funannotate update`; the identical bug would hit `.gtf`-based updates too (`GTF::Gtf_annot_retriever`, same missing directory).

**Why it matters**: The core migration goal for this session (containerized TRAIN + standalone GeneMark ES/ET + ANI-driven species-level reuse) is now genuinely validated end-to-end on real data, not just stub-run/unit-tested. `funannotate update`'s container path is NOT ready, but the fix is precise and tiny (one `install.sh` block), not a deep PASA/container incompatibility.

**Resolution**: TRAIN/GENEMARK_RUN/PREDICT confirmed working, no further action needed there. `FUNANNOTATE_UPDATE` containerization is blocked on a one-line fix to `hyphaltip/PASApipeline`'s `scripts/install.sh` (add a `cp -r "${PASA_ROOT}/SAMPLE_HOOKS" "${SRC_DIR}/"` block alongside the existing `PerlLib`/`pasa_conf`/etc. copies, ~line 162) — `pixi_install_pasa.sh` delegates entirely to `install.sh`, so no Dockerfile change is needed once that lands; just rebuild. Tracked as T-028. `run_update`/`FUNANNOTATE_ANNOTATE` were only just wired to the container in this same session (previously module-load only, untested even on host) — recommend leaving `run_update` off in production configs until the fork is patched and a new beta is built.

**How to apply (meta)**: When a "the config must be wrong" hypothesis is tempting, check the actual upstream source tree before proposing a fix — the first-pass diagnosis here (move `.pm` into `PerlLib/`) would have been a wrong, wasted fix; the real bug was one directory away, in the install script, not the Perl module layout or the Nextflow-side config.

**How to apply**: When containerizing a tool with multiple sub-modes (train/predict/update/annotate), validate each mode with a REAL run, not just a syntax/stub-run — this PASA hook bug and the earlier `run_sra_fetch` staleness gotcha were both invisible to `-stub-run` and both only surfaced by actually executing the pipeline against real strains with a genuinely dirty prior-run state.

**Tags**: funannotate, container-migration, genemark, T-022, pasa, funannotate-update, packaging-bug, beta5, ani-reuse, species-level-reuse-validated

### [2026-08-14] Correction: `predict_results/*.parameters.json`'s `genemark` field does NOT reveal whether GENEMARK_RUN used ES or ET — it's predict.py's own unrelated internal bookkeeping

**Category**: gotcha (self-correction — the previous learnings entry, same lineage, recommended exactly this check as a validation method)

**What happened**: After fixing the `run_sra_fetch` misconfiguration and re-running the beta.5 confirmation (real training now genuinely happened, `genome_annotation_training/` populated, sibling ab-initio reuse worked correctly end to end), the representative strain (*Penicillium citrinum* B8014, which has real RNA-seq and `genemark_mode=ET` set) still showed `"source": "selftraining ES"` in `parameters.json`. Read `predict.py` directly inside the beta.5 container to find the actual cause before concluding ET was broken:

```python
if genemarkcheck:
    if "path" in trainingData["genemark"][0]:
        RunModes["genemark"] = "pretrained"
        ...
    else:
        RunModes["genemark"] = "selftraining"
        trainingData["genemark"] = [{
            ...
            "source": RunModes["genemark"] + " " + args.genemark_mode,
            ...
        }]
```

`args.genemark_mode` is `funannotate predict`'s **own internal `--genemark_mode` CLI flag** (default `ES`), which `FUNANNOTATE_PREDICT/main.nf` never sets — completely unrelated to the external `GENEMARK_RUN` process's ES/ET choice or the `--genemark_gtf` file it hands to predict. This whole block only runs when `genemarkcheck` is true (gmes_petap.pl present on the host running predict, which it currently is, since predict itself isn't containerized yet) and reflects predict.py's *own* would-be internal genemark training/reuse path — it has no knowledge of what mode produced the externally-supplied GTF. The label is effectively a constant ("selftraining ES") for every run today, regardless of whether `GENEMARK_RUN` actually ran `--ES`, `--ET`, or `--predict_with`.

**Why it matters**: This invalidates a check I had just recommended in the immediately-preceding learnings entry ("spot-check `predict_results/*.parameters.json`'s `genemark`...`source` fields"). That check is fine for confirming training happened *at all* (vs. the fully-untrained fallback with generic BUSCO-lineage Augustus params), but it cannot distinguish ES from ET, and gives false confidence either way. There is currently **no persistent, trustworthy provenance** for which GeneMark mode actually ran for a given strain in a real (non-stub, non-smoke-test) pipeline run — `GENEMARK_RUN`'s own informative `[INFO] ... fresh ET self-training ...` log line only exists in the task's `.command.log`, which `cleanup = true` deletes after a successful run.

**Resolution**: Not yet fixed. The standalone smoke test (`nextflow/genemark_run_smoke.nf`, real end-to-end ET run, 10,780 genes matching the hand-validated recipe) remains the only actual evidence GENEMARK_RUN's ET branch works — this production-shaped confirmation run neither confirms nor refutes whether the *wiring* (`trainingTranscriptBamFor(out)` receiving a non-empty bam at the right time inside `FUNANNOTATE_PREDICTION.nf`) delivered ET correctly for a representative strain with real RNA-seq.

**How to apply**: Don't trust `parameters.json`'s `genemark.source` field for ES-vs-ET provenance — it's not measuring that. Future work should have `GENEMARK_RUN` write its own small persistent provenance marker (e.g. a one-line file alongside its `.genemark.gtf`/`.genemark.mod` output, published rather than left in the ephemeral work dir) so which branch actually ran is directly checkable without disabling `cleanup` or re-deriving it from `.nextflow.log`.

**Tags**: funannotate, genemark, genemark-et, provenance, container-migration, T-022, false-signal, beta5

### [2026-08-14] Verified `-w genemark:1` weight override works correctly when `gmes_petap.pl` is genuinely absent (the real containerized-predict scenario)

**Category**: validation (positive result)

**What happened**: Before wiring `FUNANNOTATE_PREDICT` into the beta.5 container, needed to verify the previously-uncertain safety net: `predict.py:567` unconditionally zeroes `StartWeights["genemark"]` when `gmes_petap.pl` isn't found on the host running predict, and `FUNANNOTATE_PREDICT/main.nf` passes an explicit `-w ... genemark:1` to force it back on — but this had never been tested against an environment where `gmes_petap.pl` is genuinely absent (predict has always run via the host module, which has it).

Ran `funannotate predict` for real, directly via `singularity exec` against the beta.5 image (no Nextflow/SLURM needed for this — a standalone test), on a small real genome (*Pichia senei*, 10 Mb) with a real precomputed `--genemark_gtf` (generated via the existing `GENEMARK_RUN` module in ES mode) and the exact `-w codingquarry:0 glimmerhmm:0 genemark:1` flags `FUNANNOTATE_PREDICT` uses. Confirmed `gmes_petap.pl`/`$GENEMARK_PATH` genuinely absent in the container first.

Traced the code path in `predict.py` before running: `StartWeights["genemark"]` initializes to 0 when `not genemarkcheck` (true here), then the `-w` argparse block runs *after* that and unconditionally sets `StartWeights[predictor.lower()] = weight` for anything the user passed explicitly — so `-w genemark:1` should survive. A second, later re-zero check (tied to `--auto-skip-genemark` on fragmented assemblies) is guarded by `if not GeneMark and ...` — and `GeneMark` (the variable) gets set truthy as soon as `--genemark_gtf` is processed, earlier in the script, so that check is skipped entirely when `--genemark_gtf` is supplied.

Real run confirmed both the debug-log `StartWeights` dict (`{'genemark': 1, ...}`, right after argparse) and the final `predict_misc/weights.evm.txt` (`ABINITIO_PREDICTION\tGeneMark\t1`) show weight 1, not 0. More importantly, the "Summary of gene models" log line showed **real GeneMark evidence flowed through**: `{'total': 12296, 'Augustus': 3404, 'HiQ': 47, 'GeneMark': 4957, 'snap': 3888}` — a substantial, plausible model count, not a non-zero weight applied to zero actual predictions. EVM produced 4,561 final consensus models.

**Why it matters**: This was the single named risk blocking `FUNANNOTATE_PREDICT`'s containerization — confirmed with real evidence, not just code-reading, that GeneMark's precomputed evidence is correctly weighted and consumed once `gmes_petap.pl` is genuinely absent from predict's own environment (the exact state predict will be in once containerized).

**Side findings while running this test** (both closed same day, see `todo/TODO_REGISTRY.md` T-027): `/opt/databases` (the container's baked-in default `$FUNANNOTATE_DB`) doesn't actually exist inside the image — confirms our explicit `export FUNANNOTATE_DB=<bind-mounted host path>` correctly takes precedence, no hidden bundled DB to worry about. Separately, all `*_odb12` BUSCO lineages in the shared `FUNANNOTATE_DB` were missing `lengths_cutoff` (funannotate's vendored old BUSCO2 script needs the legacy dataset layout; `*_odb10` and older are fine) — turned out to have zero production impact since real `samples.csv` rows all use bare/unsuffixed lineage names resolving to separate, complete legacy directories; the `*_odb12` dirs were moved out of `FUNANNOTATE_DB`'s top level afterward so nobody can select one by accident.

**Resolution**: `FUNANNOTATE_PREDICT/main.nf` wired into the container (`singularity exec` via `SING`/`SING_BINDS`, mirroring `FUNANNOTATE_TRAIN`), same day. Tagged `pre-funannotatepredimage` beforehand as a rollback point.

**How to apply**: When a safety-net code path (like a weight override) has never actually been exercised in the target environment, don't assume code-reading is sufficient before flipping a production switch that depends on it — run the real command in the real target environment first, and check for a persistent, literal, environment-derived artifact (`weights.evm.txt`, not a log line that could be a leftover/unrelated default) plus a plausible non-zero model count as positive evidence, not just "didn't error."

**Tags**: funannotate, genemark, evm-weights, container-migration, predict, validation, beta5

### [2026-08-14] `FUNANNOTATE_PREDICT` and `FUNANNOTATE_ANNOTATE` wired into beta.5, `PREDICT` validated end-to-end via Nextflow

**Category**: implementation / validation

**What happened**: With the `-w genemark:1` weight override confirmed working with real evidence (previous entry), wired `FUNANNOTATE_PREDICT/main.nf` and `FUNANNOTATE_ANNOTATE/main.nf` into the beta.5 container, mirroring `FUNANNOTATE_TRAIN`'s `SING`/`SING_BINDS` pattern exactly (binds: `target`, `training_target` for the `training/` symlink, `augustus_config`, `funannotate_db`, plus `proteins` for predict and `sbt_template` for annotate; defensive `unset -f which`). GeneMark logic left untouched -- it was already fully externalized via `GENEMARK_RUN` + `--genemark_gtf`.

Tagged `pre-funannotatepredimage` (annotated git tag) as a rollback point immediately before this change, per explicit request, since PREDICT had never run via container before and the weight-override risk was the one thing standing in the way.

Validated `-stub-run` clean for both, then forced a real re-predict for `Pichia_senei_UFMG-CM-Y531` (deleted its existing host-module `predict_results`/`predict_misc`, reran with `-resume` so only that one strain re-executed) through the real Nextflow-driven wiring -- not just the hand-rolled `singularity exec` test from the previous entry. Confirmed via the persistent predict log: `/venv/bin/funannotate predict ...` (the container binary path, not the host module's `/opt/linux/.../pixi/envs/default/bin/funannotate`), 4,647 genes (vs. 4,646 from the earlier host-module run of the same strain -- consistent, real).

**Runtime comparison** (real data, same strain/inputs, `Penicillium_citrinum_B8014` TRAIN): container **3,141s** (13:26:06-14:18:27, from the persistent `funannotate-train.log` timestamps) vs. host-module historical **3,788s** (`analysis/funannotate_train_stage_timing/outputs/per_run_summary.csv`) -- **container ~17% faster**, consistent with the rust-optimized Trinity/PASA build this image was made for.

**Resolution**: `FUNANNOTATE_UPDATE` (wired earlier, never exercised) and `FUNANNOTATE_ANNOTATE` (just wired, never exercised at all, containerized or not) both need a real test before calling container migration complete -- queued as a follow-up run (`run_update: true`, `run_annotate: true` added to `params_confirm.yaml`) against the same 4-strain isolated test dir.

**Tags**: funannotate, container-migration, predict, annotate, beta5, git-tag, runtime-comparison, validation

### [2026-08-18] LZ-ANI 1.2.3 is unusable as a genome-level ANI method for fungal genomes — silent failures + OOM quirks

**Category**: tool-evaluation / negative result

**What happened**: Evaluating `lz-ani all2all` (quay.io biocontainers `lz-ani:1.2.3--h9ee0642_0`, needs `module load singularity/3.9.3` or `4.3.2` — `singularity` is NOT on PATH by default) against real fungal scaffold files (32–52 MB) revealed hard failure modes with no clean way to detect them:

- `--multisample-fasta true` (the tool's default): produces a clean 3-column query/reference/ANI TSV, but at **scaffold level** — every scaffold becomes a sample (`NW_022983500.1`-style contig IDs), so it is semantically wrong for genome ANI and NOT comparable to skani/fastani/mash genome ANI.
- `--multisample-fasta false` (correct genome-level mode): **silently writes NO output on real fungal genomes, even single-genome runs** — `exit 0`, normal logs, ~4.9 GB RSS, nothing emitted. This is the worst failure class: silent data absence with no distinguishable log line. (On tiny synthetic contigs it does produce output, so tests that only use small inputs look green.)
- Thread-scaled genome-level runs OOM-kill (`exit 137`) at `-t >= 4` even for 1–2 genomes on the 4-core / 503 GB test host; ~5 GB peak RSS per 1–2 genomes.
- `--out-format` is ignored on small inputs — the tool emits a **sectioned `[lz_similarities]` file** indexed by sample number instead of a plain table.

**Resolution**: LZ-ANI de-scoped from the compare_ANI pipeline (module, config, schema, params, README all cleaned); evidence + reproduction recipe kept in `analysis/ani_method_evaluation/LZANI_PERFORMANCE.md`; the method-agnostic `scripts/compare_ani_methods.py` was kept. Container SIF left in the shared cache.

**Tags**: lz-ani, ani, container, singularity, tool-evaluation, negative-result, oom, silent-failure, de-scope

### [2026-08-18] `pathlib.Path.glob()` raises `NotImplementedError` on absolute patterns — use `glob.glob()` for absolute globs

**Category**: python gotcha

**What happened**: In `scripts/compare_ani_methods.py`, `Path(trace_dir).glob('/var/lib/.../ANI_trace.*.txt')` (absolute pattern) crashed with `NotImplementedError: Non-relative patterns are unsupported` — an easy-to-miss break when plumbing an absolute `--trace-glob` through a configurable parameter (stdlib `glob.glob()` handles both relative and absolute patterns).

**Resolution**: switched to `glob.glob(os.path.expanduser(trace_glob))` so the same code path serves relative CWD-relative globs, `~/`-expansion, and absolute globs. (Note: `nextflow_schema.json` uses the top-level `$defs` key, *not* `definitions` — the schema's param groups live under `$defs/ani_options/properties`.)

**Tags**: python, pathlib, glob, gotcha, compare-ani-methods

### [2026-08-28] Containerized `wgd syn` (i-ADHoRe) needs OpenMPI env hygiene inside the SIF — SLURM env leaks kill it

**Category**: hpc / container gotcha

**What happened**: Validating the WGD_SYN module (`wgd syn` → i-ADHoRe for collinear-block anchor detection) in the paralogoscope pipeline. First real failure before the fix: `FileNotFoundError: [Errno 2] No such file or directory: 'wgd_syn/iadhore-out/anchorpoints.txt'` with `i-adhore: error while loading shared libraries: libmpi.so.40` in the log. The lib IS in the container at `/opt/conda/lib/libmpi.so.40` (conda openmpi 4.1.6) but `LD_LIBRARY_PATH` isn't set (base prefix, not an activated env). With `LD_LIBRARY_PATH=/opt/conda/lib` set, i-ADHoRe then failed with `OPAL ERROR: Unreachable in file pmix3x_client.c at line 111` — the classic "launched under srun/tores, but OMPI was not built with SLURM's support" error: SLURM env vars from the job shell leak into the container and OMPI tries (and fails) to talk to SLURM's PMI. Also needed `OMPI_ALLOW_RUN_AS_ROOT=1` + `OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1` (container runs as root).

**Resolution**: In the module's `script:` block, run the command as `apptainer exec ... --env LD_LIBRARY_PATH=/opt/conda/lib --env OMPI_ALLOW_RUN_AS_ROOT=1,OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 ${wgd_sif} bash -c 'env -u <SLURM_*> -u <PMI_*> bash -c "wgd syn ..."'`. Explicitly unset all `SLURM*` and `PMI*` env vars inside the container before the wgd call (a static list of ~29 vars; the dynamic `for v in $(env | grep -E '^SLURM|^PMI' | cut -d= -f1); do unset $v; done` was verified by hand but the explicit list is baked in). Verified live: i-ADHoRe prints its banner and completes (`All Done!  Bye...`); the full dmd→ksd→syn Nextflow run succeeds 6/6 and publishes dmd/ksd outputs.

**Gotchas learned along the way**:
- `env -u FOO bar -c "cmd"` fails (`env: invalid option -- 'c'`) — the `-c` belongs to the *program* being invoked, so it must be `env -u FOO ... bash -c "cmd"`.
- **"No anchors found" is a LEGITIMATE rc=0 outcome** of `wgd syn` (i-ADHoRe ran, found no collinear anchors, `cli.py:716` warns and terminates) — in that case `anchors.csv` is NEVER written. So a module guard like `[ -f anchors.csv ] || exit 1` would falsely fail valid runs. Guard must be tolerant: require anchors.csv **or** a non-empty `iadhore-out/anchorpoints.txt` (i.e. i-ADHoRe actually ran).
- `write_genelists` **skips scaffolds with ≤2 genes** (`if len(sdf.index) <= 2: continue` in `/opt/conda/lib/python3.9/site-packages/wgd/syn.py` ~line 118) → i-ADHoRe "ERROR: Genelist files not found in settings file" on tiny test genomes.
- Overlapping (non-partition) MCL-like families cause "ERROR: There are duplicated gene IDs for given feature and attribute" — must feed a valid partition.
- **`wgd syn` needs `-f mRNA -a ID` for funannotate-style GFF3s**: CDS-transcript headers (`FF5840CF_00001-T1`) match mRNA features, NOT gene features; the default `-f gene -a ID` errors "No genes from families file ... found in the GFF file".
- wgd syn visualizations ARE named `*.dot.{svg,pdf,png}` — confirmed in cli.py:740-742 (the `all_dotplots` loop). (The `*_multiplicons_level.*`/`All_species_marcosynteny.*` names from viz.py belong to a different code path, `macrosynteny`, not `syn`.) So WGD_SYN's `wgd_syn/*.dot.pdf` and `*.dot.png` publishDir/output patterns are correct. Also confirmed: `wgd_syn/anchors.csv` IS the real output name (cli.py:720) written only when anchors are found; `*.anchors.ks.tsv` + `*.ksd.{svg,pdf}` come from a `-k/--ks-distribution` variant of syn (not used by default).
- Empirically, a "No anchors found" run exits **0** even though cli.py:716 calls `exit(1)` — the click wrapper normalizes it. Do not rely on rc; rely on output files.

**Tags**: wgd, i-ADHoRe, openmpi, apptainer, container, SLURM, PMI, nextflow, paralogoscope, gotcha

### [2026-08-28] Nextflow `publishDir` glob patterns don't cross directories — scope patterns to the publishDir root

**Category**: nextflow gotcha

**What happened**: In the paralogoscope wgd modules, publishDir was `{ "${params.outdir}/${meta.id}/wgd_dmd" }` with `pattern: '*.tsv'` — but `*` does not cross directory separators, so the published `wgd_dmd/` dirs existed while no files were ever copied (silent success, nothing published). Varying it to `{ outdir/{meta.id}/wgd_dmd }` + pattern `'wgd_dmd/*.tsv'` double-nested the dir (`outdir/{id}/wgd_dmd/wgd_dmd/file.tsv`).

**Resolution**: Set each module's publishDir root to `{ "${params.outdir}/${meta.id}" }` and express the pattern relative to THAT root: `pattern: 'wgd_dmd/*.tsv'`, `'wgd_ksd/*.ks.tsv'`, `'wgd_ksd/*.pdf'` (ksd's actual pdf is `*.tsv.ksd.pdf`). Verified real-run tree: `outdir/{id}/wgd_dmd/{id}.cds-transcripts.fa.tsv`, `wgd_ksd/{id}.cds-transcripts.fa.tsv.ks.tsv`, `wgd_ksd/{id}.cds-transcripts.fa.tsv.ksd.pdf`.

**Tags**: nextflow, publishDir, glob, pattern, paralogoscope, gotcha

### [2026-08-28] Login-node smoke tests are capped at 4 CPUs / 24 GB; local config must match

**Category**: HPCC gotcha

**What happened**: A real-data wgd run on the login node failed twice: `Process requirement exceeds available CPUs -- req: 8; avail: 4` then, after lowering cpus to 4, `Process requirement exceeds available memory -- req: 32 GB; avail: 24 GB`. The login node is 4 vCPU / 24 GB and is shared with other sessions.

**Resolution**: The local test config caps `comparative_wgd` at cpus=4, memory='16 GB', time='2 h'. SLURM runs (the real deployment) use 8 cpus/32 GB via the paralogoscope profile — no cap there.

**Tags**: HPCC, login-node, slurm, resources, nextflow, paralogoscope, gotcha

### [2026-08-28] `pkill -f <pattern>` can kill its own invoking shell

**Category**: shell gotcha

**What happened**: `pkill -f "plg_stub"` matched the current bash command line (which contained the pattern as an argument) and killed the tool's own shell, hanging the call until timeout.

**Resolution**: Never `pkill -f` a pattern that appears in the command you're typing. Scope it (`pgrep -af java | grep -v <self-excluding-token>`) or match a more specific string absent from the command line (e.g. the full java `-jar ...nextflow...jar` arg).

**Tags**: shell, pkill, pgrep, gotcha

### [2026-08-28] Nextflow `include` relative paths resolve from the module file, not the workflow

**Category**: nextflow gotcha

**What happened**: `include { hashBucketForType } from '../../common/utils.nf'` in `modules/comparative/wgd/WGD_DMD/main.nf` resolved to `modules/comparative/common/utils.nf` (Error: Invalid include source), because the wgd process sits THREE dirs under `modules/` while BFD's `modules/BFD/CALC_INTERGENIC/` sits TWO.

**Resolution**: Modules nested deeper under `modules/` need one more `../`: `from '../../../common/utils.nf'`. Compilation fails loudly (unlike the publishDir-glob silent failure), so it never ships broken.

**Tags**: nextflow, include, relative-path, modules, gotcha

### [2026-08-28] `storeDir` preserves output relative subpaths — move outputs to workdir root for clean buckets

**Category**: nextflow gotcha

**What happened**: With `storeDir { outdir/wgd_dmd/{bucket} }` and outputs left in the process's `wgd_dmd/` subdir, storeDir would place files at `outdir/wgd_dmd/{bucket}/wgd_dmd/file.tsv` (double nesting), exactly like the earlier publishDir double-nest learning.

**Resolution**: `mv` the declared outputs to the workdir root in the script (`mv wgd_dmd/*.tsv .`) before the processes complete; then storeDir drops them directly in the bucket dir. Note: with `-stub-run`, storeDir is NOT populated — stub validation checks parsing/flow, not the bucket layout.

**Tags**: nextflow, storeDir, publishDir, glob, paralogoscope, gotcha

### [2026-08-28] Backgrounded `setsid nohup nextflow &` dies when the invoking session/shell exits

**Category**: HPCC gotcha

**What happened**: Three detached nextflow launches (setsid + nohup + disown, `</dev/null` redirected) all died within ~30 s of the invoking tool call returning, with no error in the log. A foreground run under `timeout` survives fine.

**Resolution**: For anything long, submit a SLURM job (`run_paralogoscope.sh` uses `sbatch` on preempt) instead of backgrounding on the login node. The r3 real-data smoke test is the exception that was let to finish in the foreground of an interactive session.

**Tags**: HPCC, nextflow, background, setsid, nohup, slurm, gotcha
