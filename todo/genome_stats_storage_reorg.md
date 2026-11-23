# Plan: reorganize genome_stats / function flat-file stores + query DB

**Status**: legacy migration complete (2026-07-30), storeDir/consumer cutover (#9) in progress.
**Owner**: TBD. **Related**: T-014 in `TODO_REGISTRY.md`.
**Reviewed by**: Fable (acting as HPC storage / data-management expert), 2026-07-30 — critique
incorporated below (see "Review notes" callouts).

**Update (2026-07-30, post-migration)**: `asm_reports/` (22,412 files) was investigated during
implementation and confirmed dead — no current Nextflow module or script produces or consumes it
(superseded by `asm_stats`/`CALC_ASM_STATS` + `summarize_asm_stats.py`). It was **deleted outright**
rather than migrated, along with its matching dead consumer `scripts/collect_asm_stats.py`. Every
`asm_reports` reference below is left as historical record of the original investigation, not
current state — wherever you see `asm_stats/asm_reports/BUSCO_genome` grouped together, read it as
`asm_stats/BUSCO_genome` only. See `.living/decisions.md` for the pruning decision.

## Current state (verified in repo, 2026-07-30)

- `results/genome_stats/` has 8 subfolders, all flat, keyed two different ways:
  - **Accession-keyed** (`ASMID`, e.g. `GCA_000149445.2_ASM14944v2`): `asm_stats/` (15,599 files),
    `asm_reports/` (22,412 files).
  - **Name-keyed** (`Genus_species_strain`, built by `cleanStrain()`/`meta.id`): `BUSCO_genome/`
    (14,544), `BUSCO_protein/` (5,954), `aa_freq/` (8,025), `codon_freq/` (8,044), `gene_stats/`
    (55,398 — 7 files/genome, ~7,914 genomes), `intergenic_stats/` (8,022).
  - Note: the original ask named `BUSCO_genome` as accession-keyed; verified it is actually
    name-keyed (`{meta.id}.BUSCO_summary.{lineage}.txt`). Plan below handles both keying schemes
    explicitly rather than assuming.
- `results/function/` has 9 subfolders (aiupred, cazy, merops, pfam_hmmscan, predgpi, signalp,
  targetP, tmhmm, wolfpsort), also name-keyed, flat.
- `samples.csv` is the join table: `ASMID, LOCUSTAG, GENUS, SPECIES, STRAIN, ...`. 329 rows have
  no `GENUS` (need a NOGENUS bucket).
- Producers are Nextflow `storeDir`-cached processes in `nextflow/modules/BFD/*` (CALC_ASM_STATS,
  BUSCO_GENOME, BUSCO_PEP, CALC_GENE_STATS, CALC_INTERGENIC, BATCH_AA_FREQ, BATCH_CODON_FREQ),
  wired via `nextflow/subworkflows/local/BFD_GENOME_STATS.nf`, each hardcoding a flat
  `"${params.genome_stats_outdir}/<type>"` storeDir. `results/function/*` producers (BFD_FUNCTIONAL)
  presumed same flat pattern — confirm during implementation.
- Consumers reading these directories directly: `nextflow/bin/pick_representative_strain.py`,
  `nextflow/bin/species_reuse_clusters.py` (both take `--busco-dir results/genome_stats/BUSCO_genome`),
  `scripts/collect_busco_stats.py`, `scripts/collect_asm_stats.py`.
- **A database layer already exists** — this is not a greenfield build:
  - `tables/All_Taxa/*.csv.gz` — per-type gzipped CSV snapshots, rebuildable per taxon subset
    (`tables/<Taxon>/`).
  - `scripts/build_BFD_duckDB.sh` — rebuilds `db/BFD.<SUBSET>.duckdb` from `tables/<SUBSET>/*.csv.gz`
    **from scratch every time** (`rm -f "$DB"` then full reload).
  - `db/BFD.duckdb` — 3.3 GB, in production use, queried by `MCP/BFD_mcp_server.py`.
  - So the "build a database" ask is largely already done; the gap is that the collect_*.py
    scripts walk directories by filename pattern, which breaks once files move into subfolders
    unless updated to be bucket-aware.

## A. Naming/framework decision

**Canonical on-disk key: `ASMID` for assembly-level outputs, `LOCUSTAG` for annotation-derived
outputs — never `Genus_species_strain`.**

- `ASMID` (assembly accession) for: `asm_stats`, `asm_reports`, `BUSCO_genome`.
- `LOCUSTAG` for: `BUSCO_protein`, `aa_freq`, `codon_freq`, `gene_stats`, `intergenic_stats`, and
  all of `results/function/*` — these are keyed off the annotated gene set, not the raw assembly.
- Rationale: both keys are unique and stable in `samples.csv`; `Genus_species_strain` is not
  (329 rows have no genus; `collect_asm_stats.py` already has a multi-step fallback cascade for
  exactly this reason). Every consumer already has `samples.csv` for the join, so switching the
  *storage* key loses no information — only the at-a-glance filename, which is solved separately.

**Directory fan-out: hash bucket, not genus grouping — width chosen per type, not uniformly 2 hex
chars.**

- Bucket = first N hex chars of `sha1(KEY)`.
- **Review note (Fable): a flat 2-hex-char (256-bucket) scheme undersizes `gene_stats`/
  `intergenic_stats` at 10x growth** (gene_stats: ~7,914 genomes × 7 files today → ~55k; at 10x
  genomes, ~550k files → ~2,150/bucket at 2 hex chars, back to slow-directory territory). Widening
  later is itself a second migration (every path helper, symlink, and manifest encodes the bucket
  width). **Decision: pick bucket width per type up front:**
  - 2 hex chars (256 buckets) for `asm_stats`, `asm_reports`, `BUSCO_genome`, `BUSCO_protein`,
    `aa_freq`, `codon_freq`, and all of `function/*`.
  - 3 hex chars (4,096 buckets) for `gene_stats` and `intergenic_stats`, given the known ~7x
    file-per-genome multiplier for gene_stats.
- Path: `results/genome_stats/<type>/<hashN>/<KEY>.<ext>`, e.g.
  `results/genome_stats/asm_stats/4e/GCA_000149445.2_ASM14944v2.stats.txt`.
- **Rejected: genus/-based grouping**, even though first suggested, because genus population is
  extremely skewed (Aspergillus/Penicillium/Fusarium have hundreds of strains; most genera have 1)
  — genus/ just relocates the flat-directory problem into a handful of hot genus folders, and 329
  genomes have no genus at all. Hash fan-out is uniform by construction.

**Human-browsable secondary view: the genus/genus_species_strain idea survives, as a generated,
read-only symlink tree — never the canonical store.**

- `results/genome_stats_by_name/<GENUS-or-NOGENUS>/<Genus_species_strain>.<type>.<ext>` →
  symlink into the real hashed path. Same idea for `results/function_by_name/`.
- **Review note (Fable) — staleness trigger must be explicit**: regenerate whenever a hash of
  (`samples.csv` mtime+checksum, source-tree top-level mtimes) changes — not "any time it's stale"
  with no defined check. Script compares this hash to a stamp file and rebuilds only on mismatch.
- **Review note (Fable) — exclude from backup/rsync**: `*_by_name/` trees are regenerable and
  should be explicitly excluded from any backup/rsync/tar pass over `results/` (rsync with
  symlink-following defaults would otherwise needlessly re-touch every target file during backup).
  Document this exclusion in `DATA_MANIFEST.md` when this ships.

**Shared path-resolution helper, one implementation per language, kept in sync** (repo already has
this pattern for `cleanStrain()` — Groovy in `common/utils.nf`, Python in
`collect_busco_stats.py::build_basename_map()`):
- Groovy: `hashBucket(key, width)` added to `nextflow/modules/common/utils.nf`.
- Python: `hash_bucket(key, width)` added to a new shared `nextflow/bin/genome_stats_paths.py`,
  imported by every collect_*.py / pick_representative_strain.py / species_reuse_clusters.py.
- Cross-check test: generate N random keys, assert Groovy and Python buckets match bit-for-bit —
  a mismatch here is a silent path miss, not a crash, same class of risk as the `join()`/`combine()`
  drift bugs already logged in `.living/learnings.md`.

## B. Production plan (Nextflow BFD + funannotate + related)

1. Add `hashBucket()`/`genome_stats_paths.py` (per-type width table baked in, per §A).
2. Update every BFD module's `storeDir` (`CALC_ASM_STATS`, `BUSCO_GENOME`, `BUSCO_PEP`,
   `CALC_GENE_STATS`, `CALC_INTERGENIC`, `BATCH_AA_FREQ`, `BATCH_CODON_FREQ`, and the
   BFD_FUNCTIONAL producers) to
   `storeDir "${params.genome_stats_outdir}/<type>/${hashBucket(key, width)}"`.
   `clearIfStale()` staleness logic is path-agnostic, no change needed there.
3. Update `BFD_GENOME_STATS.nf`'s `planFreq()` expected-output path construction to include the
   bucket.
4. Update every consumer to take a pre-built manifest (key → path) rather than globbing — the repo
   already hit this exact glob-vs-race class of bug (`.living/learnings.md`,
   `join()`/`combine()`/`Channel.fromPath` entries). `pick_representative_strain.py --busco-dir`,
   `species_reuse_clusters.py --busco-genome-dir`, `collect_busco_stats.py
   --genome-busco-dir/--pep-busco-dir`, `collect_asm_stats.py --genomedir` all need
   `os.path.join(dir, hash_bucket(key), filename)`.
5. Do not touch `tables/*.csv.gz` layout or `db/BFD.duckdb` build mechanics in this pass beyond the
   glob fix in (4) — see §D for the separate Parquet follow-on.
6. **Review note (Fable) — this IS a breaking cutover, call it out explicitly, don't gloss over
   it**: changing `storeDir` changes the *declared output path* Nextflow's `-resume` cache keys on.
   Every already-cached BUSCO/asm_stats/gene_stats/etc. task will look "not done" to Nextflow after
   this change regardless of whether the legacy migration (§C) has already physically moved the
   files to the new paths — `-resume` cache identity and "does this file exist at this path" are
   different mechanisms. **Decision required before merging step 2**: either (a) accept one full
   recompute pass across all cached BFD stats processes after cutover (expensive — these are
   HPC-hour-costly BUSCO/stats jobs across ~8-55k units), or (b) run the legacy migration (§C)
   *before* deploying the new code so `storeDir`'s existence check (not the resume cache) finds the
   file already present at the new path and skips recompute — `storeDir` semantics: if the declared
   output file already exists, Nextflow does not resubmit that task, independent of `-resume`. This
   plan recommends (b): migrate first (§C), deploy code second, so no recompute is triggered.
7. **Review note (Fable) — check for in-flight runs at cutover**: confirm no Nextflow BFD run is
   mid-flight against the old flat `storeDir` paths when the new code is deployed; a run started
   under old paths and interrupted mid-cutover risks orphaned or duplicated outputs. Quiesce/finish
   or explicitly kill and restart any in-flight run before deploying step 2's code.

## C. One-time legacy migration

New script `scripts/one-off/reorg_genome_stats_hash_buckets.py` (per repo convention: one-off →
`scripts/one-off/`, defaults to dry-run, `--apply` required to write):

1. Load `samples.csv` → `{ASMID: LOCUSTAG, GENUS, SPECIES, STRAIN}` and the reverse
   `Genus_species_strain → ASMID/LOCUSTAG` map, reusing (factored out, not reimplemented) the same
   fallback cascade already coded in `collect_asm_stats.py`.
2. For each of the 17 flat source dirs: for every existing file, resolve its key (ASMID for
   asm_stats/asm_reports/BUSCO_genome; LOCUSTAG for the rest), compute `hash_bucket(key, width)`
   per §A's per-type width table, record `(old_path, new_path)`.
3. **Review note (Fable) — check target-path uniqueness BEFORE any rename executes**, not just via
   a post-hoc total-count check: two different keys landing on the same bucket is fine, but if two
   *different* keys somehow produce the same target filename (shouldn't happen given keys are
   unique, but verify), fail loudly before writing anything, not after.
4. Files whose key can't be resolved (no samples.csv match) go to `results/_UNMATCHED/<type>/` —
   never silently dropped, never renamed. Given 329 NOGENUS rows and known historic species-name
   drift, expect a non-zero unmatched count; this is a required, tested output of the dry run
   (**review note: the test-profile dry run must explicitly cover the `_UNMATCHED` path**, not just
   the happy-path bucket assignment).
5. Move, don't copy (`os.rename` within same filesystem — instant metadata op). Write the full
   manifest (`old_path,new_path,sha256,size,mtime`) to a CSV *before* any rename, so the migration
   is auditable and reversible from the manifest alone.
6. **Review note (Fable) — make the migration idempotent/resumable**: an sbatch job over ~130k
   renames that dies partway (OOM, walltime, preemption — all documented risks on this cluster)
   must be safely re-runnable. Before acting on each planned rename, check "is old_path already gone
   and new_path already present matching the manifest row" — if so, skip as already-done, don't
   error or reprocess.
7. Verification pass required before declaring done (robust-analysis convention: fail loudly, assert
   counts): total file count before == after (old+new+unmatched), spot-check N random files' sha256
   before/after match, per-bucket file counts roughly uniform (flag any bucket >3x the mean as a
   hash-function bug, not "just how it fell").
8. **Review note (Fable) — check disk/inode headroom before running**: even metadata-only renames
   can transiently pressure inode-count quotas (symlinks count too, relevant once step (9) runs).
   Confirm quota headroom before starting.
9. Execution: single-threaded `os.rename` is metadata-only and fast even at 130k files, but run via
   `sbatch -p short` (not the login node — matches existing learning about not running heavy
   filesystem/process work on shared login nodes), with periodic progress logging; do one
   `os.walk` up front rather than repeated globbing (NFS listing here is documented as slow).
10. After migration: run the human-facing symlink-tree generator (§A) once to produce
    `results/genome_stats_by_name/` and `results/function_by_name/`.
11. Rollback: the manifest CSV is sufficient to reverse (swap old_path/new_path and re-run the
    idempotent rename step) if something is found wrong after the fact.

## D. Database recommendation: DuckDB stays; Parquet replaces CSV.gz, unpartitioned by default

Keep DuckDB as the query/serving layer — proven, has an MCP server built on it, does SQL joins
fluently across genome_stats + `samples.csv`-derived `species` table. Nothing here argues for
replacing it.

Change the staging format from gzipped CSV (`tables/*.csv.gz`) to **Parquet, one file per data
type, unpartitioned**:
- DuckDB reads Parquet natively, faster than CSV (typed columns, no re-inference, better
  compression for numeric-heavy tables like asm_stats).
- **Review note (Fable) — do NOT partition Parquet by hash bucket by default**: partitioning a
  table-of-rows dataset into 256 small Parquet files per type reintroduces the exact small-file NFS
  metadata problem this whole plan exists to solve, one layer up. The current DB is 3.3GB total;
  these tables are small enough that a single unpartitioned Parquet file per type already beats
  gzipped CSV on read speed and size with none of the small-file cost. Reserve bucket-partitioned
  Parquet only if/when a single type's table genuinely grows large enough that DuckDB's partition
  pruning becomes worth the file-count cost — not by default, and not in this pass.
- One Parquet file per data type (asm_stats, busco, aa_freq, codon_freq, gene_stats, intergenic,
  function/*), each loaded as a DuckDB table/view — matches "one table per type, co-located in one
  DB, linked via samples.csv" exactly; `species` (from `samples.csv`) stays the join spine.
- Independent, separately-shippable from §B/§C — the existing CSV.gz → full-rebuild approach keeps
  working throughout the file-layout migration.

### D.1 Update strategy, non-duplication, and staleness (added 2026-07-31, pre-#11 design discussion)

**Staleness relative to a regenerated annotation is already solved upstream, at the file layer —
the Parquet swap does not need to solve it again.** When a genome is re-annotated (new GFF/genome
mtime), `clearIfStale()` (`BFD_GENOME_STATS.nf`) deletes the stale bucketed stat file; the next BFD
run recomputes it fresh at the same bucket+key path. `MERGE_*` reads its inputs via
`gatedGlobIn()` → `toManifest()`, which bakes each source file's `(path, mtime, size)` into the
manifest *content* — so a changed stat file changes the manifest hash, which is Nextflow's own
`-resume` cache key. `MERGE_*` is therefore already incrementally cache-aware end-to-end: it only
re-runs when something underneath it actually changed, and when it does re-run it writes a
*complete, fresh* table for its type from every genome currently known — never an append.
Switching its output format from `tables/*/type.csv.gz` to `tables/*/type.parquet` is a pure format
change at the tail of an already-correct pipeline; no new staleness logic needed at this layer.

**Non-duplication is structural, not validated after the fact.** `collect_busco_stats.py` /
`summarize_asm_stats.py` / the `MERGE_*` scripts all accumulate into a dict keyed by ASMID or
LOCUSTAG (`records[asmid] = stats`) before writing output — a re-annotated genome's fresh stats
*overwrite* the same key. Each `MERGE_*` run produces exactly one row per genome. As long as the
rule stays "each `MERGE_*` run writes a complete, fresh Parquet file for its type, never appends,"
duplication is impossible by construction.

**Decision: full rebuild from Parquet at both layers, not incremental upsert.**
- `MERGE_*` → `tables/*/type.parquet`: already full-rewrite-per-run (see above); only the output
  format changes.
- `tables/*.parquet` → `db/BFD.duckdb`: keep `build_BFD_duckDB.sh`'s existing `rm -f "$DB"` +
  full `CREATE TABLE AS SELECT` shape, just point it at Parquet instead of CSV.gz. DuckDB reading
  typed, columnar Parquet is materially faster than re-parsing gzipped CSV, so this step gets
  *faster*, not more complex, for free.

**Rejected for now: `INSERT ... ON CONFLICT (key) DO UPDATE` incremental upsert.** DuckDB supports
it, but Parquet itself is immutable/columnar, not designed for row-level updates — doing this
properly needs a staging/delta Parquet file for changed rows plus a periodic compaction step that
rewrites the main Parquet file (`COPY (SELECT * FROM main WHERE key NOT IN (SELECT key FROM delta)
UNION ALL SELECT * FROM delta) TO 'main.parquet'`), a mini Delta-Lake/Iceberg pattern. Real
complexity not justified at current scale (tens of thousands of rows per table, 3.3GB DB total) —
revisit only if `build_BFD_duckDB.sh`'s full-rebuild time becomes an actual bottleneck in practice.

**Considered and rejected: DuckDB `VIEW`s over Parquet instead of materialized `TABLE`s**
(`CREATE VIEW aa_freq AS SELECT * FROM read_parquet('tables/*/aa_freq.parquet')`). This eliminates
the rebuild step entirely — every query re-reads current Parquet content, zero staleness concept
at the DB layer at all. Rejected as the default because it has no indexes, so repeated/interactive
queries (the MCP server's use case) are slower per-query than a materialized, indexed table; the
MCP server is presumably latency-sensitive, so pay the (now cheap, Parquet-sourced) full-rebuild
cost instead. Worth keeping in mind as a lower-maintenance escape hatch if the rebuild step is ever
skipped/forgotten in practice — a set of views could serve as a safety-net "always current" layer
underneath, independent of whether anyone remembers to run the table rebuild.

### D.2 Retire per-taxon merge re-runs; master-only build + separate extraction tool (added 2026-07-31)

**Decision**: stop producing per-taxon `tables/<Taxon>/*` subset folders by re-running `MERGE_*`
scoped to `--taxon`. `taxonRowFilter()` is applied at the *base* genome channel (`BFD.nf`, before
any per-genome computation), which means the hash-bucketed per-genome store (T-014) is already
universal/cumulative — a genome computed under any `--taxon` scope or a full run lands in the same
shared bucketed directories regardless. The current `tables/<Taxon>/` re-merge only exists as a
side effect of scoping the manifest glob to avoid pulling in genomes outside the requested clade
(`BFD_MERGE.nf`'s own comment); it was never a deliberately designed "small scoped artifact"
feature. Once the master merge is always current (glob mode, covers every genome ever computed),
a SQL filter over master data gives the exact same correctness guarantee at far lower cost — no
re-glob/re-cat/re-dedup of the whole per-genome store for every taxon someone wants to look at.

- **`MERGE_SAMPLES` simplifies**: builds one always-complete master `samples`/`species` table
  (`subset_samples.py`/`build_species_table.py` still used, just always unscoped — no more
  `--taxon`-filtered subset calls).
- **New tool, its own issue, not folded into #11**: `EXTRACT_BFD_TAXONOMIC_SUBSET` — a Python CLI
  wrapping DuckDB (`ATTACH`/`COPY`), e.g. `extract_bfd_taxonomic_subset.py --taxon GENUS:Aspergillus
  --master db/BFD.duckdb --out db/BFD.Aspergillus.duckdb`. Not a Nextflow process (pure SQL
  filtering over already-materialized data, no per-genome parallelism to gain, would only add
  scheduler overhead) and not a bare `.sql` script (a CLI gives argument validation, logging,
  testability that `duckdb -c` alone doesn't).
- **Correctness requirement, not optional**: `species`/`samples` is the *authoritative* key list.
  Every other table (`aa_freq`, `gene_stats`, `asm_stats`, all 9 function tables) must be filtered
  by `INNER JOIN` against that key set, never by independently re-evaluating the taxon predicate
  per table — otherwise a schema quirk or missing taxonomy field in one table can silently produce
  a subset DB with orphaned rows (a function-table row for a genome that got filtered out of
  `species`). Add a post-extract row-count assertion (matching-key row count per table vs.
  `species`) before writing the output file, and fail loudly per this repo's robust-analysis
  convention.
- **No new staleness logic needed**: since the master build is already "full rebuild every time"
  (§D.1), the extractor just runs its filter fresh against whatever `db/BFD.duckdb`/master Parquet
  currently is — no separate caching or version-stamping for the extract itself.
- **Reviewed by Fable (2026-07-31)**: flagged one question requiring explicit confirmation before
  implementing — does the current per-taxon build ever serve as an access boundary (a taxon-scoped
  DB handed to a collaborator who should specifically not see other taxa's data), not just a
  performance artifact? **Confirmed by the user: purely performance/convenience, never an access
  boundary** — the extractor can be built as a straightforward correctness-focused (not
  security-focused) tool.
- **Verified before deciding**: grepped for other consumers of the per-taxon pattern.
  `MCP/BFD_mcp_server.py` already defaults to the single master `db/BFD.duckdb`
  (`PROTEIN_DB_PATH` env var, no hardcoded per-taxon path) — the primary query interface is
  unaffected by retiring per-taxon builds. `scripts/Rscripts/*.R` reference table paths and should
  get a quick spot-check before implementation, but nothing found suggests a hardcoded
  taxon-subset dependency.
- **Sequencing**: keep `EXTRACT_BFD_TAXONOMIC_SUBSET` out of #11. Land the Parquet conversion
  first (master-only, glob-mode build, simplified `MERGE_SAMPLES`), then build the extractor as
  its own issue against a stable Parquet/DuckDB target rather than a moving one.

## E. Stepwise delivery (issues/branches/PRs)

1. Issue + branch: `hashBucket()`/`genome_stats_paths.py` + Groovy/Python parity unit test — no
   behavior change yet.
2. Issue + branch: `scripts/one-off/reorg_genome_stats_hash_buckets.py` — dry-run reviewed against
   production tree (including a specific test of the `_UNMATCHED` path and idempotent re-run), then
   a separately-tracked `--apply` run. **This runs before step 3**, per the §B.6 sequencing decision
   (migrate first, so the storeDir cutover finds files already in place and triggers no recompute).
3. Issue + branch: update BFD module `storeDir`s + `BFD_GENOME_STATS.nf` planning logic + all
   consumer scripts to bucket-aware paths, validated against `tests/output/genome_stats` (test
   profile) with `-stub-run` first, then a real small run against the now-migrated production tree
   to confirm zero unexpected recomputation.
4. Issue + branch: symlink-tree generator for `*_by_name/` human views, with the staleness-hash
   trigger and an explicit backup/rsync exclusion note added to `DATA_MANIFEST.md`.
5. Issue + branch (independent, can run in parallel with 1-4): Parquet staging conversion
   (unpartitioned) + incremental-vs-full-rebuild decision for `db/BFD.duckdb`.
6. Update `.living/decisions.md` (ASMID/LOCUSTAG + per-type hash width over genus grouping, and the
   migrate-before-cutover sequencing decision), `.living/learnings.md` if anything unexpected turns
   up during migration, and the relevant `*_MANIFEST.md` files per CLAUDE.md's post-action protocol.

## Open items before starting implementation

- Confirm `results/function/*` producer modules (BFD_FUNCTIONAL) actually follow the same flat
  `storeDir` pattern assumed here — not yet directly inspected.
- Decide who/when runs the `--apply` migration (§C) — recommend scheduling during a window with no
  in-flight BFD runs, given §B.7's in-flight-run risk.
- Confirm whether any external tooling (MCP server, R scripts, notebooks) reads
  `results/genome_stats/*` paths directly rather than through `tables/`/`db/` — if so, those need
  updating too and aren't currently enumerated above.
