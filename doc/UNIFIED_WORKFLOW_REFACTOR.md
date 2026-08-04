# Unified Workflow Refactor Plan

## Goal

Go from `samples.csv` to annotated genomes in a single command:

```
samples.csv → CLEAN → MASK → ANI → ASM_STATS → BUSCO → BFD_MERGE
          → PICK_REPRESENTATIVE → RNASEQ → TRAIN → PREDICT → ANNOTATE
```

## Current Architecture: 3 Manually-Sequenced Pipelines

```
┌─────────────────────────────────────────────────────────────────────┐
│ Pipeline 1: BFD (--pipeline BFD)                                     │
│   samples.csv → GENOME_STATS(ASM_STATS, BUSCO, TELOMERES)           │
│              → BFD_MERGE → tables/{asm_stats,busco_genome}.parquet  │
└─────────────────────────────────────────────────────────────────────┘
                                  ↓ (manual: user runs 3 commands)
┌─────────────────────────────────────────────────────────────────────┐
│ Pipeline 2: compare_ANI (--pipeline compare_ani)                    │
│   samples.csv → ANI_COMPARE → CONCAT_ANI_TSVS                       │
│              → PICK_REPRESENTATIVE_STRAIN                            │
│              → abinitio_reuse_assignments.csv                        │
└─────────────────────────────────────────────────────────────────────┘
                                  ↓
┌─────────────────────────────────────────────────────────────────────┐
│ Pipeline 3: funannotate (--pipeline funannotate)                    │
│   samples.csv → GENOME_PREP(CLEAN, MASK)                            │
│              → RNASEQ → PREDICTION(TRAIN, PREDICT)                  │
│              → ANNOTATION                                           │
└─────────────────────────────────────────────────────────────────────┘
```

### Implicit file-on-disk dependencies (the user must sequence these)

| Producer pipeline | File on disk | Consumer |
|---|---|---|
| BFD (`MERGE_ASM_STATS`) | `tables/asm_stats.parquet` | `pick_representative_strain.py` (PICK_REPRESENTATIVE_STRAIN) |
| BFD (`MERGE_BUSCO_GENOME`) | `tables/busco_genome.parquet` | `pick_representative_strain.py` (PICK_REPRESENTATIVE_STRAIN) |
| compare_ANI (`PICK_REPRESENTATIVE_STRAIN`) | `genome_annotation/_reuse_assignments/abinitio_reuse_assignments.csv` | `loadAbinitioReuseMap()` in funannotate.nf |

### Why 3 pipelines exist separately

The current split reflects a historical ordering: BFD was built first (functional annotation + genome stats), compare_ANI was added for representative selection, and funannotate was ported from 1KFG. Each has its own profile config (`conf/profile_BFD.config`, `conf/profile_ANI.config`, `conf/profile_funannotate.config`) with different `workDir`, `singularity` settings, and param defaults.

## Dependency Analysis: The Critical Path

The unified pipeline's DAG is:

```
samples.csv
    │
    ▼
GENOME_PREP (CLEAN + MASK)
    │
    ▼ (masked + cleaned genomes)
    ├──→ ANI (ANI_COMPARE_METHOD)
    │        → CONCAT_ANI_TSVS → merged ANI TSV ─────────────────┐
    │                                                             │
    ├──→ ASM_STATS (CALC_ASM_STATS) → per-genome stats ──┐       │
    │                                                     ├→ BFD_MERGE
    ├──→ BUSCO (BUSCO_GENOME) → per-genome BUSCO ────────┘       │
    │                                                     ↓       │
    │                                              parquet tables │
    │                                                     ↓       │
    │                                          PICK_REP_STRAIN ◄──┘
    │                                                     ↓
    │                                          abinitio_reuse_assignments.csv
    │                                                     ↓
    │                                          READ_REUSE_MAP → abinitioReuseMap (channel)
    │                                                     ↓
    ├──→ RNASEQ (SRA fetch + Trinity) ◄──────────────────────────
    │           ↓
    │      TRAIN (funannotate train)
    │           ↓
    │      PREDICT (funannotate predict)
    │           ↓
    └──→ ANNOTATE (funannotate annotate)
```

### Key findings

1. **ANI, ASM_STATS, and BUSCO can all fork in parallel** after CLEAN+MASK. No dependency between them.
2. **BFD_MERGE must wait for ASM_STATS + BUSCO** to produce per-genome outputs.
3. **PICK_REPRESENTATIVE_STRAIN is the barrier** — it needs both the merged ANI TSV and the parquet tables from BFD_MERGE. Only after it finishes can RNASEQ start.
4. **RNASEQ needs the reuse map** to pick the right assembly for the shared Trinity-GG build. This means RNASEQ is gated on PICK_REPRESENTATIVE_STRAIN.
5. **FUNANNOTATE_ANNOTATION is the terminal sink** — no emits, writes annotation outputs to disk.

## The One Hard Problem: `loadAbinitioReuseMap()` Is Synchronous

```groovy
// funannotate.nf line 141
def abinitioReuseMap = loadAbinitioReuseMap()  // reads CSV from disk NOW
FUNANNOTATE_RNASEQ(FUNANNOTATE_GENOME_PREP.out.predict_genome, abinitioReuseMap)
```

This function reads `abinitio_reuse_assignments.csv` from disk at **workflow construction time** — before any process runs. In the current 3-pipeline flow, this works because Pipeline 2 (compare_ANI) already produced the CSV on disk.

In a unified pipeline, the CSV wouldn't exist until `PICK_REPRESENTATIVE_STRAIN` runs mid-pipeline. So `loadAbinitioReuseMap()` would return an empty map, and all strains would train independently — defeating the purpose of ANI-driven reuse.

**This is the single barrier to unification.** Everything else is just channel wiring.

## Two Options

### Option A — Temporary: Launcher Script (no code changes)

A `run_unified.sh` that runs the 3 pipelines in sequence with `-resume`. Each phase is independently resumable.

**Pros:** Zero code changes, works today.
**Cons:** Three Nextflow driver processes, three sets of work dirs, no cross-phase channel wiring, manual error recovery between phases.

### Option B — Refactor: Unified Workflow (`--pipeline unified`)

The real refactor. The changes needed are:

#### B.1. Convert `loadAbinitioReuseMap()` from a function to a process

The core change. Instead of reading the CSV synchronously at workflow construction time, emit it as a channel:

```groovy
// New process: READ_REUSE_ASSIGNMENTS
// Input:  path assignments_csv  (from PICK_REPRESENTATIVE_STRAIN.out.outCSV)
// Output: val reuse_map_json    (JSON string of the map)
process READ_REUSE_ASSIGNMENTS {
    input:  path assignments_csv
    output: val reuse_map_json
    script:
    """
    python3 - << 'PYEOF'
    import csv, json
    m = {}
    with open('${assignments_csv}') as f:
        for row in csv.DictReader(f):
            m[row['out']] = {
                'species': row['species'],
                'reuse_eligible': row['reuse_eligible'].lower() == 'true',
                'is_representative': row.get('is_representative', '').lower() == 'true',
            }
    print(json.dumps(m))
    PYEOF
    """
}
```

Then in the workflow:
```groovy
// PICK_REPRESENTATIVE_STRAIN emits assignments_csv
// READ_REUSE_ASSIGNMENTS converts it to a channel value
// FUNANNOTATE_RNASEQ and FUNANNOTATE_PREDICTION consume the channel
```

This is the deep module change — it moves the reuse map from a synchronous disk read to an asynchronous channel, allowing it to be produced mid-pipeline and consumed downstream.

#### B.2. Create `workflows/unified.nf`

```groovy
workflow UNIFIED {
    // ── Build jobs channel from samples.csv ────────────────────────────
    def jobs = buildJobsChannel()

    // ── Phase 1: GENOME_PREP (CLEAN + MASK) ────────────────────────────
    FUNANNOTATE_GENOME_PREP(jobs)
    def predict_genome = FUNANNOTATE_GENOME_PREP.out.predict_genome

    // ── Phase 2: Fork into ANI, ASM_STATS, BUSCO ───────────────────────
    //    All need the cleaned/masked genome from Phase 1
    ANI_REPRESENTATIVE_SELECT(predict_genome)
    def merged_ani_tsv = ANI_REPRESENTATIVE_SELECT.out.merged_tsv

    BFD_GENOME_STATS(
        aa_freq_ch, codon_freq_ch, intergenic_ch, gene_stats_ch,
        asm_stats_ch, telomeres_ch, busco_genome_ch, busco_pep_ch
    )

    // ── Phase 3: BFD_MERGE (after ASM_STATS + BUSCO) ───────────────────
    BFD_MERGE(
        BFD_GENOME_STATS.out.aa_cached,
        BFD_GENOME_STATS.out.aa_csv,
        // ... other channels ...
        BFD_GENOME_STATS.out.busco_genome,
        use_glob
    )

    // ── Phase 4: PICK_REPRESENTATIVE_STRAIN (after ANI + BFD_MERGE) ────
    //    Needs: merged ANI TSV, parquet tables, predict_input
    //    Produces: abinitio_reuse_assignments.csv
    PICK_REPRESENTATIVE_STRAIN(
        merged_ani_tsv.ifEmpty(file('/dev/null')),
        writePredictInput(predict_genome),
        file(params.samples),
        file(params.tables)
    )

    // ── Phase 5: Load reuse map from CSV (now it exists!) ──────────────
    def reuse_map_ch = READ_REUSE_ASSIGNMENTS(
        PICK_REPRESENTATIVE_STRAIN.out.outCSV
    ).out.reuse_map_json

    // ── Phase 6: RNASEQ (gated on reuse map) ───────────────────────────
    FUNANNOTATE_RNASEQ(predict_genome, reuse_map_ch)

    // ── Phase 7: PREDICTION (gated on reuse map + RNASEQ output) ───────
    FUNANNOTATE_PREDICTION(
        FUNANNOTATE_RNASEQ.out.predict_input,
        reuse_map_ch,
        forceIndependentSet
    )

    // ── Phase 8: ANNOTATION ────────────────────────────────────────────
    FUNANNOTATE_ANNOTATION(
        FUNANNOTATE_PREDICTION.out.metadata,
        FUNANNOTATE_RNASEQ.out.reads,
        taxonFilter,
        asmidFilter
    )
}
```

#### B.3. Interface design (codebase-design vocabulary)

The **deep module** here is the `UNIFIED` workflow itself. Its **interface** is:

| Interface element | What the caller provides |
|---|---|
| `--pipeline unified` | One flag activates the entire pipeline |
| `samples.csv` | One input file |
| `--run_ani_reuse true` | One flag controls ANI-driven reuse |
| All other params | Same params as the 3 separate pipelines |

The **implementation** hides:
- The 3-pipeline orchestration (CLEAN→MASK→ANI→BUSCO→MERGE→PICK_REP→RNASEQ→PREDICT→ANNOTATE)
- The file-on-disk dependencies between pipelines
- The `loadAbinitioReuseMap()` synchronous-read problem (converted to a channel)
- The parallel forking of ANI/ASM_STATS/BUSCO after CLEAN+MASK

The **seam** is at `PICK_REPRESENTATIVE_STRAIN.out.outCSV → READ_REUSE_ASSIGNMENTS`. This is where the reuse map transitions from a file-on-disk (produced by the ANI phase) to a channel value (consumed by the prediction phase). This seam is necessary because the reuse map is produced mid-pipeline and consumed downstream — it can't be a synchronous function call.

The **depth** comes from:
- **Leverage:** The caller gets the entire 9-step pipeline from one `--pipeline unified` command
- **Locality:** All the orchestration complexity is concentrated in `unified.nf`, not spread across 3 separate entry points

#### B.4. Backward compatibility

The 3 separate pipelines (BFD, compare_ANI, funannotate) would continue to work unchanged. The unified workflow would be a new entry point that chains them together.

The only shared code change is converting `loadAbinitioReuseMap()` from a function to a process. This change would be backward compatible because:
- The existing `funannotate.nf` calls `loadAbinitioReuseMap()` at workflow construction time (reads from disk)
- The new `unified.nf` would use `READ_REUSE_ASSIGNMENTS` (reads from a channel)
- Both can coexist — the function stays for the standalone `funannotate` pipeline, the process is new for the `unified` pipeline

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| `loadAbinitioReuseMap()` conversion breaks existing funannotate pipeline | Low | High | Keep the function for standalone funannotate; add the process only for unified |
| ANI compute takes too long before RNASEQ can start | Medium | Medium | ANI and BUSCO run in parallel; RNASEQ is gated but not unnecessarily |
| Profile config conflicts (BFD vs ANI vs funannotate) | Medium | High | Create a new `conf/profile_unified.config` that merges all three |
| Work dir conflicts between phases | Low | Medium | Use separate work dirs per phase (already the pattern in existing profiles) |

## Implementation Roadmap

### Phase 0: Temporary launcher (Option A) — now

1. Write `nextflow/run_unified.sh` — a launcher script that runs the 3 pipelines in sequence
2. Each phase uses its own profile and `-resume`
3. Error recovery: just re-run the failed phase

### Phase 1: Unified workflow (Option B) — future

1. Create `conf/profile_unified.config` that merges params from all 3 profiles
2. Create `workflows/unified.nf` that chains all the steps together
3. Add `--pipeline unified` to `nextflow/main.nf`
4. Convert `loadAbinitioReuseMap()` to `READ_REUSE_ASSIGNMENTS` process
5. Test with `-stub-run` first, then real data

### Phase 2: Profile consolidation — future

1. Merge `conf/profile_BFD.config`, `conf/profile_ANI.config`, `conf/profile_funannotate.config` into a single `conf/profile_unified.config`
2. Keep the 3 separate profiles for backward compatibility
3. The unified profile should be a superset of all 3

## Assessment: Strategy to Go from samples.csv → Annotated Genomes

**Recommended approach:** Start with Option A (launcher script) as a temporary measure, then refactor to Option B (unified workflow) when ready.

**Option A is the right temporary approach because:**
1. Zero code changes — works today
2. Each phase is independently resumable with `-resume`
3. Error recovery is simple — just re-run the failed phase
4. The launcher script can be enhanced to check pre-requisites (e.g., verify parquet tables exist before running funannotate)

**Option B is the right long-term refactor because:**
1. Single Nextflow driver process — simpler monitoring
2. Cross-phase channel wiring — no file-on-disk dependencies
3. Single `-resume` — one command recovers from any failure point
4. The `loadAbinitioReuseMap()` → `READ_REUSE_ASSIGNMENTS` conversion is the only deep change

**The key risk in Option B** is the `loadAbinitioReuseMap()` conversion. If the reuse map channel is not properly gated, RNASEQ could start before `PICK_REPRESENTATIVE_STRAIN` finishes, leading to incorrect reuse assignments. The fix is to ensure `READ_REUSE_ASSIGNMENTS` is a process (not a function) that emits the reuse map as a channel value, which naturally gates downstream consumers.

**My recommendation:** Start with Option A — it's a ~50-line script that works today. Then, when ready to invest in the refactor, Option B is a well-scoped change (one new workflow file, one new process, one function-to-process conversion) that delivers the unified pipeline.
