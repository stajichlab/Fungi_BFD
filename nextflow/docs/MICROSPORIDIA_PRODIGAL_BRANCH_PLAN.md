# Microsporidia Prodigal Branch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

> **Revision history:** v1 of this plan (a 6-process, `.combine()`/`.filter()`/`.mix()`
> three-way channel rewrite of `FUNANNOTATE_PREDICTION.nf`) was reviewed and found to
> have 5 blocking defects, verified against real data — most critically, the gate as
> designed would never have fired on the genome the recipe was validated on
> (*Ordospora colligata* OC4 scores `ok`, not `small_fragmented`, under the composite
> preflight verdict: small but not fragmented), and a missing `.ifEmpty([:])` would
> have silently dropped every genome in the run whenever `microsporidiaSet` is empty
> (i.e. by default, for every non-microsporidia pipeline invocation). v2 (this
> document) adopts the reviewer's simpler design: fold the decision into
> `GENEMARK_RUN` itself, where the preflight stats are already computed in bash, with
> no new process and no channel restructuring.

**Goal:** For genomes whose samples.csv `PHYLUM` is `Microsporidia` AND whose total
assembled size is below `params.predict_min_asm_bp`, add a self-trained Prodigal
ab-initio track alongside GeneMark-ES in `funannotate predict` (via `--other_gff`),
with AUGUSTUS/SNAP forced off and EVM partitioning disabled — reproducing the recipe
already validated end-to-end in
`../../Microsporidia_predict/scripts/run_microsporidia_predict.py` against
*Ordospora colligata* OC4 (`genemark:1` + `prodigal:5`, `augustus:0`, `snap:0`,
`--no-evm-partitions`, `--min_protlen 30`).

**Architecture:** `GENEMARK_RUN` gains one new input (`is_microsporidia`, a plain
boolean computed once per row from a taxonomy Set) and, when true, additionally runs
Prodigal in the same task and emits a second, always-present (possibly empty) output
path alongside its existing `gtf`/`mod` outputs. `FUNANNOTATE_PREDICT` gains one new
input (`other_gff`) and, when it's non-empty, adds `--other_gff <path>:5`, zeroes
`augustus`/`snap` weights, adds `--no-evm-partitions --min_protlen 30`, and — critically
— bypasses its own small/fragmented preflight skip (which would otherwise discard the
Prodigal evidence before `funannotate predict` ever runs). GeneMark is still attempted
on these genomes exactly as before; if it fails, the existing empty-GTF degradation
path already handles that, so no GeneMark-skip logic changes. No new Nextflow process,
no channel `.filter()`/`.mix()` restructuring, no `.combine()` broadcast.

**Tech Stack:** Nextflow DSL2, Groovy, Prodigal (`module load prodigal/2.6.3`,
confirmed available on this cluster), funannotate (existing `funannotate_sif`
container).

**Spec:** This document. Validated reference implementation (read-only, not modified
by this plan):
`../../Microsporidia_predict/scripts/run_microsporidia_predict.py` and
`.../funannotate/data/profiles/microsporidia-default.json`.

## Global Constraints

- Zero behavior change for any non-Microsporidia genome, and zero behavior change
  for the whole pipeline when `params.microsporidia_phylum` is unset/empty. Verify
  this concretely in Task 5 (existing stub-run must be byte-for-byte identical in
  which processes get invoked).
- Reuse `bin/asm_preflight_stats.py`'s already-computed `ASM_BP` value for the size
  gate — do not use the composite `ASM_VERDICT` (`small_fragmented`) for this
  decision; it requires fragmentation too, which the validated target genome (OC4:
  2.29 Mb, 15 contigs, N50 228,601 — small but *not* fragmented) does not have.
  Confirmed by direct run: `python bin/asm_preflight_stats.py
  input_clean_genomes/GCF_000803265.1_ASM80326v1.masked.fasta.gz --min-bp 8000000
  --max-n50 10000 --max-contigs 1000` → `2290528	15	228601	ok`.
- The taxonomy Set must be keyed the same way the pipeline's `out` identifier is
  built: `makeSampleTag(row.SPECIES, row.STRAIN)` (see
  `workflows/funannotate.nf:110`) — **not** `SPECIES_IN`, which is a different,
  raw/unnormalized column.
- Do not skip GeneMark for these genomes. Attempt it exactly as today; let the
  existing empty-`genemark.gtf` degradation absorb a failure. This matches what was
  actually validated (`genemark:1` AND `prodigal:5` together), not a GeneMark-free
  variant.
- Follow existing house style: `record_*_skip()`-style TSV provenance logs under
  `${params.target}`; synchronous CSV-reading helpers in `modules/common/utils.nf`
  (see `loadSuppressSet()`) for anything read once and reused as a plain Groovy
  `Set`/`Map`.
- No `nf-test` framework exists in this repo. Testing convention is
  `-profile <pipeline>,test,test_<pipeline> -stub-run` against a fixture CSV in
  `tests/data/<pipeline>/`.
- Every `-w` weight decision must go through a single `-w` group (funannotate's
  `argparse` replaces the whole list on a second `-w` occurrence — see the existing
  comment in `FUNANNOTATE_PREDICT/main.nf:203-207`).
- Changing `GENEMARK_RUN`'s and `FUNANNOTATE_PREDICT`'s input tuple arity/script body
  invalidates the Nextflow cache for **every** genome on the next `-resume`, not just
  microsporidia ones (task hash includes both). This is unavoidable for this class of
  change; call it out to whoever runs the next `-resume` after this deploys, and
  avoid deploying mid-run if at all possible.

---

### Task 1: Taxonomy gate — `loadMicrosporidiaOutSet()` + new params

**Files:**
- Modify: `nextflow/modules/common/utils.nf` (add function, after `loadSuppressSet()` at line 277)
- Modify: `nextflow/conf/profile_funannotate.config` (new params, near `predict_min_asm_bp` block at line 216)
- Modify: `nextflow/workflows/funannotate.nf` (call the loader, log it, pass to `FUNANNOTATE_PREDICTION`)
- Modify: `nextflow/nextflow_schema.json` — add `microsporidia_phylum` param doc

**Interfaces:**
- Produces: `Set<String> loadMicrosporidiaOutSet()` — a Groovy function reading
  `params.samples` synchronously, returning the `out` values (built via
  `makeSampleTag(SPECIES, STRAIN)`, matching exactly how `workflows/funannotate.nf`
  builds `out` for the main `jobs` channel) whose `PHYLUM` column
  (case-insensitive, trimmed) equals `params.microsporidia_phylum`. Empty set when
  `params.samples` doesn't exist, `PHYLUM` column absent, or no rows match.

- [ ] **Step 1: Add the new params to `conf/profile_funannotate.config`**

Add immediately after the existing pre-flight guard block:

```groovy
    // Set predict_min_asm_bp = 0 to disable the pre-flight guard entirely.
    predict_min_asm_bp        = 8000000    // assembled bp below this = "small"
    predict_frag_max_n50      = 10000      // N50 below this = "fragmented"
    predict_frag_max_contigs  = 1000       // contig count above this = "fragmented"

    // Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    // genomes whose samples.csv PHYLUM matches this value AND whose total
    // assembled bp is below predict_min_asm_bp (size alone -- NOT the
    // small_fragmented composite verdict, which also requires fragmentation
    // and would miss small-but-complete genomes like Ordospora colligata OC4,
    // the genome this recipe was validated on) get a self-trained Prodigal
    // track alongside GeneMark-ES, fed to funannotate predict via
    // --other_gff. Set to '' to disable the whole feature.
    microsporidia_phylum      = 'Microsporidia'
```

- [ ] **Step 2: Write `loadMicrosporidiaOutSet()` in `modules/common/utils.nf`**

Insert directly after `loadSuppressSet()` (ends at line 277):

```groovy
// Build the set of `out` values (see makeSampleTag()) whose samples.csv
// PHYLUM matches params.microsporidia_phylum. Mirrors loadSuppressSet()'s
// synchronous-read pattern: read once, return a plain Set for O(1)
// .contains() checks in FUNANNOTATE_PREDICTION.nf's per-row closures.
//
// MUST key off SPECIES/STRAIN (via makeSampleTag), matching exactly how
// workflows/funannotate.nf builds `out` for the jobs channel (line 110) --
// NOT SPECIES_IN, a separate raw/unnormalized column that produces a
// different, non-matching tag.
//
// Plain comma-split (not a full CSV parser) matches every other synchronous
// loader in this file (loadSuppressSet, loadAbinitioReuseMap in
// modules/funannotate/utils.nf) -- if embedded commas in quoted fields ever
// bite here, fix it in all of them together, not just this one.
def loadMicrosporidiaOutSet() {
    def phylum = (params.microsporidia_phylum ?: '').trim()
    if (!phylum) {
        return ([] as Set)
    }
    def f = file(params.samples as String)
    if (!f.exists()) {
        return ([] as Set)
    }
    def lines = f.readLines()
    if (lines.size() < 2) {
        return ([] as Set)
    }
    def header  = lines[0].split(',', -1)*.trim()
    def iSpecies = header.indexOf('SPECIES')
    def iStrain  = header.indexOf('STRAIN')
    def iPhylum  = header.indexOf('PHYLUM')
    if (iSpecies < 0 || iStrain < 0 || iPhylum < 0) {
        log.warn "loadMicrosporidiaOutSet: samples.csv missing SPECIES/STRAIN/PHYLUM column; microsporidia_phylum gate is a no-op"
        return ([] as Set)
    }
    def outSet = lines.drop(1)
        .collect { it.split(',', -1) }
        .findAll { f2 -> f2.size() > [iSpecies, iStrain, iPhylum].max() }
        .findAll { f2 -> f2[iPhylum].trim().equalsIgnoreCase(phylum) }
        .collect { f2 -> makeSampleTag(f2[iSpecies].trim(), f2[iStrain].trim()) }
        .toSet()
    if (outSet) {
        log.info "loadMicrosporidiaOutSet: ${outSet.size()} genome(s) tagged PHYLUM=${phylum}"
    }
    return outSet
}
```

- [ ] **Step 3: Wire it into `workflows/funannotate.nf`**

Add the import to the existing line 15:

```groovy
include { makeSampleTag; loadSuppressSet; suppressRowFilter; loadMicrosporidiaOutSet } from '../modules/common/utils.nf'
```

Add directly after the `forceIndependentGenemarkSet` block (ends around line 79):

```groovy
    // Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    // out values that get a Prodigal track alongside GeneMark-ES when also
    // below predict_min_asm_bp.
    def microsporidiaSet = loadMicrosporidiaOutSet()
```

Update the single `FUNANNOTATE_PREDICTION(...)` call site (near the end of the file),
together with Task 4 Step 1's `take:` block change (do these two edits in the same
commit — the call site and the `take:` arity must change together or the pipeline
won't compile):

```groovy
    FUNANNOTATE_PREDICTION(FUNANNOTATE_RNASEQ.out.predict_input, abinitioReuseMap, forceIndependentSet, forceIndependentGenemarkSet, microsporidiaSet)
```

- [ ] **Step 4: Add the param to `nextflow_schema.json`**

Find the block documenting `predict_min_asm_bp`/`predict_frag_max_n50`/
`predict_frag_max_contigs` (grep for `predict_frag_max_contigs` in
`nextflow_schema.json`) and add a sibling entry:

```json
"microsporidia_phylum": {
    "type": "string",
    "default": "Microsporidia",
    "description": "samples.csv PHYLUM value that, combined with predict_min_asm_bp, gates the Prodigal ab-initio supplement. Empty string disables the feature."
}
```

- [ ] **Step 5: Commit**

```bash
git add nextflow/modules/common/utils.nf nextflow/conf/profile_funannotate.config nextflow/workflows/funannotate.nf nextflow/nextflow_schema.json
git commit -m "Add microsporidia_phylum taxonomy gate"
```

(This commit alone will not compile standalone — `workflows/funannotate.nf`'s call
site now passes 5 args to `FUNANNOTATE_PREDICTION`'s still-4-arg `take:` block.
Land this together with Task 4 Step 1, or as one combined commit, not deployed
in between.)

---

### Task 2: Prodigal hierarchy-builder script

**Files:**
- Create: `nextflow/bin/prodigal_hierarchy.py`

**Interfaces:**
- Produces: `prodigal_hierarchy.py RAW_PRODIGAL.gff3 OUT_HIER.gff3` — a standalone
  script, called from `GENEMARK_RUN`'s script body in Task 3.

- [ ] **Step 1: Port the hierarchy-builder from the validated orchestrator**

Create `nextflow/bin/prodigal_hierarchy.py`:

```python
#!/usr/bin/env python3
"""Re-emit a Prodigal CDS-only GFF3 as full gene/mRNA/CDS blocks.

EVM assembles gene models from gene/mRNA blocks; CDS-only features can only
act as support, never define structure. (An earlier version of this
docstring cited a specific Ordospora colligata OC4 Sn comparison for this
transform, attributed to Microsporidia_predict PLAN.md 9.15 -- that number
does not actually appear in PLAN.md or STATUS.md and was a fabricated
citation; removed 2026-09-04. The gene/mRNA/CDS-vs-CDS-only structural
argument above is a mechanical fact about EVM's input format, independent
of any specific benchmark number.) Each single-exon Prodigal CDS
becomes its own gene/mRNA/CDS block, blocks separated by a blank line
(EVM's lib.readBlocks in funannotate-runEVM.py expects every gene block
immediately preceded by one -- no ##gff-version header, which would join
the first block and crash gene_blocks_to_interlap).

Usage: prodigal_hierarchy.py RAW_PRODIGAL.gff3 OUT_HIER.gff3
"""
import sys


def main():
    src, dst = sys.argv[1], sys.argv[2]
    n = 0
    with open(src) as fh, open(dst, "w") as out:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 9 or cols[2] != "CDS":
                continue
            seqid, _, _, start, end, score, strand, phase, _ = cols
            n += 1
            gid, mid = f"prodigal_g{n}", f"prodigal_m{n}"
            out.write(f"{seqid}\tprodigal\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={gid}\n")
            out.write(f"{seqid}\tprodigal\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={mid};Parent={gid}\n")
            out.write(f"{seqid}\tprodigal\tCDS\t{start}\t{end}\t{score}\t{strand}\t{phase}\tID={mid}.cds;Parent={mid}\n")
            out.write("\n")
    print(f"[INFO] prodigal_hierarchy: wrote {n} gene/mRNA/CDS blocks to {dst}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Commit**

```bash
git add nextflow/bin/prodigal_hierarchy.py
git commit -m "Add prodigal_hierarchy.py (standalone, not yet called)"
```

---

### Task 3: Extend `GENEMARK_RUN` to also run Prodigal for gated genomes

**Files:**
- Modify: `nextflow/modules/funannotate/predict/GENEMARK_RUN/main.nf`
- Modify: `nextflow/subworkflows/local/FUNANNOTATE_PREDICTION.nf` (the three
  `*_genemark_in` map closures and the three `${x}_with_gtf` join chains)

**Interfaces:**
- Consumes: `prodigal_hierarchy.py` (Task 2), `microsporidiaSet` (available in
  `FUNANNOTATE_PREDICTION.nf` once Task 4 Step 1 adds it to `take:` — do Task 3 and
  Task 4 Step 1 together, since Task 3's closures reference `microsporidiaSet`).
- Produces: `GENEMARK_RUN` input tuple grows from 10 to 11 fields (adds
  `is_microsporidia` as the last field); output gains a new always-emitted
  `tuple val(out), path("${out}.other.gff3"), emit: other_gff` (empty file when
  `is_microsporidia` is false or the genome doesn't qualify — never
  `optional: true`, so downstream joins stay a simple 1:1 by `out`, matching the
  existing `gtf` output's own "always emit, sometimes empty" convention rather
  than the `mod` output's `optional: true` convention).

- [ ] **Step 1: Extend the input tuple** (replace lines 68-71)

```groovy
    input:
    tuple val(out), val(asmid), val(species), val(strain),
          val(genome_fa), val(transl_table), val(mode), val(training_bam),
          val(force_independent), val(shared_mod), val(is_microsporidia)
```

- [ ] **Step 2: Add the new output** (add after line 74, before the `mod` output)

```groovy
    tuple val(out), path("${out}.genemark.gtf"), emit: gtf
    // Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    // always emitted (empty file when is_microsporidia=false or this genome
    // didn't qualify), never optional -- same "always emit, sometimes empty"
    // convention as `gtf` above, so downstream joins stay simple 1:1 instead
    // of needing optional-output handling.
    tuple val(out), path("${out}.other.gff3"), emit: other_gff
```

- [ ] **Step 3: Add the Prodigal supplement, right after the existing pre-flight guard block**

Insert after line 177 (right after the existing
`if [ "\$ASM_VERDICT" = "small_fragmented" ]; then ... exit 0; fi` block — genomes
tripping that block `exit 0` before reaching this new code, which is intentional:
`FUNANNOTATE_PREDICT`'s own preflight guard would discard their Prodigal evidence
anyway per Task 4 Step 3's note. The microsporidia gate below is evaluated
independently of `ASM_VERDICT`, using `ASM_BP` directly, since `ASM_VERDICT`
requires fragmentation too and would miss small-but-complete genomes like OC4):

```bash
    # ── Microsporidia Prodigal supplement ────────────────────────────────────
    # Gated on PHYLUM (is_microsporidia, computed once in Groovy from
    # loadMicrosporidiaOutSet()) AND total assembled bp -- NOT the
    # small_fragmented composite verdict above, which also requires
    # fragmentation and would miss the genome this recipe was validated on
    # (Ordospora colligata OC4: 2.29 Mb, 15 contigs, N50 228,601 -- small but
    # not fragmented). Runs alongside, not instead of, the GeneMark attempt
    # below: the validated recipe used genemark:1 AND prodigal:5 together --
    # if GeneMark independently fails on this genome, the existing empty-GTF
    # degradation further down already absorbs that; this block does not
    # change GeneMark's own success/failure handling at all.
    touch "${out}.other.gff3"
    if [ "${is_microsporidia}" = "true" ] && [ "\$ASM_BP" -lt "${params.predict_min_asm_bp}" ]; then
        echo "[INFO] GENEMARK_RUN ${out}: microsporidia + small (\$ASM_BP bp < ${params.predict_min_asm_bp}); running Prodigal supplement"
        module load prodigal
        prodigal -i genome.fa -o "${out}.prodigal.raw.gff3" -f gff -p single -g "${transl_table}"
        python "${workflow.projectDir}/bin/prodigal_hierarchy.py" \\
            "${out}.prodigal.raw.gff3" "${out}.other.gff3"
        PRODIGAL_REPORT="${params.target}/prodigal_selected_microsporidia.tsv"
        mkdir -p "${params.target}"
        [ -s "\$PRODIGAL_REPORT" ] || printf 'out\tasmid\treason\ttotal_bp\n' > "\$PRODIGAL_REPORT"
        printf '%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "microsporidia_phylum_and_small" "\$ASM_BP" >> "\$PRODIGAL_REPORT"
    fi
```

Note: this reads `genome.fa` (already inflated earlier in the script at lines
126-130) and `\$ASM_BP` (already computed by the existing preflight call at lines
167-170) — no duplicate genome inflation or preflight computation needed. This
block must run for every genome that reaches it regardless of the
`small_fragmented` verdict, since `is_microsporidia`-gated genomes that are
small-but-not-fragmented (like OC4) fall through the existing
`if [ "\$ASM_VERDICT" = "small_fragmented" ]` check without exiting.

- [ ] **Step 4: Update the `stub:` block** (replace lines 295-301)

```groovy
    stub:
    """
    touch "${out}.genemark.gtf"
    touch "${out}.other.gff3"
    if [ -z "${shared_mod}" ] || [ "${force_independent}" = "true" ]; then
        touch "${out}.genemark.mod"
    fi
    """
```

- [ ] **Step 5: Update the three call sites in `FUNANNOTATE_PREDICTION.nf` to pass `is_microsporidia`**

Each of the three `*_genemark_in` map closures (lines 107-110, 172-175, 315-319)
gains one field. E.g. the first (representative-only branch):

```groovy
            def rep_genemark_in = rep_todo.map { out, asmid, sp, st, lt, bl, hl, tt, gfa, _shared_json ->
                def forceIndep = forceIndependentGenemarkSet.contains(out as String) ? 'true' : 'false'
                def isMicro    = microsporidiaSet.contains(out as String) ? 'true' : 'false'
                tuple(out, asmid, sp, st, gfa, tt, genemarkMode, trainingTranscriptBamFor(out as String), forceIndep, '', isMicro)
            }
```

Apply the identical one-line addition (`def isMicro = microsporidiaSet.contains(out as String) ? 'true' : 'false'`,
appended as the tuple's 11th field) to the other two occurrences (the
`rep_and_indep`/`genemark_in` map and the `sibling_predict_todo`/`sib_genemark_in`
map).

- [ ] **Step 6: Thread `GENEMARK_RUN.out.other_gff` alongside `.out.gtf` in the same three `.join()` calls**

Each of the three `${x}_with_gtf` constructions currently does
`.join(GENEMARK_RUN.out.gtf)` (or `GENEMARK_RUN_SIB.out.gtf`). Chain a second join
onto `.out.other_gff` (also always-emitted, 1:1 by `out`, safe for the same reason
`.out.gtf` already is). E.g. for the representative-only branch (replace lines
112-115):

```groovy
            rep_with_gtf = rep_todo.join(GENEMARK_RUN.out.gtf).join(GENEMARK_RUN.out.other_gff)
                .map { out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf, other_gff ->
                    tuple(out, asmid, sp, st, lt, bl, hl, tt, gfa, shared_json, gtf, other_gff)
                }
```

Apply the identical pattern to the other two `_with_gtf` constructions (lines
177-180 and 321-324), and to the corresponding `else` branches (when `runGenemark`
is false) so every `_with_gtf` tuple consistently ends `(..., gtf, other_gff)` —
the `else` branches simply set `other_gff = ''` alongside the existing `gtf = ''`.

- [ ] **Step 7: Commit**

```bash
git add nextflow/modules/funannotate/predict/GENEMARK_RUN/main.nf nextflow/subworkflows/local/FUNANNOTATE_PREDICTION.nf
git commit -m "GENEMARK_RUN: run Prodigal supplement for gated microsporidia genomes"
```

(This leaves `FUNANNOTATE_PREDICT` unable to consume the new `other_gff` field yet
— Task 4 finishes the wiring. Do not deploy between Task 3 and Task 4.)

---

### Task 4: Extend `FUNANNOTATE_PREDICT` to consume `other_gff` and bypass its own guard

**Files:**
- Modify: `nextflow/subworkflows/local/FUNANNOTATE_PREDICTION.nf` (`take:` block)
- Modify: `nextflow/modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf`

**Interfaces:**
- Produces: `FUNANNOTATE_PREDICT` input tuple grows from 11 to 12 fields:
  `tuple val(out), val(asmid), val(species), val(strain), val(locustag),
  val(busco_lineage), val(header_length), val(transl_table), val(genome_fa),
  val(shared_params_json), val(genemark_gtf), val(other_gff)`.
- `FUNANNOTATE_PREDICTION`'s `take:` block gains `microsporidiaSet` as a 5th
  parameter (consumed by Task 3 Step 5's closures — must be declared here too).

- [ ] **Step 1: Add `microsporidiaSet` to `FUNANNOTATE_PREDICTION.nf`'s `take:` block** (lines 48-55)

```groovy
    take:
    predict_input_ch             // tuple(out, asmid, sp, st, lt, bl, hl, tt, genome_fa)
    abinitioReuseMap             // out -> [species, reuse_eligible, is_representative], loaded offline
    forceIndependentSet          // species that always train independently (whole ab-initio bundle)
    forceIndependentGenemarkSet  // out values that always run GeneMark independently, finer-grained
                                  // than forceIndependentSet -- AUGUSTUS/SNAP reuse can stay on for a
                                  // strain while forcing GeneMark specifically to retrain for it
    microsporidiaSet             // out values gated (with predict_min_asm_bp) to also get a Prodigal
                                  // track from GENEMARK_RUN (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md)
```

- [ ] **Step 2: Extend `FUNANNOTATE_PREDICT`'s input tuple** (replace lines 17-20)

```groovy
    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          val(genome_fa), val(shared_params_json), val(genemark_gtf), val(other_gff)
```

- [ ] **Step 3: Bypass the existing preflight skip when `other_gff` is non-empty**

Modify the guard at lines 160-168 (the `if [ "\$ASM_VERDICT" = "small_fragmented" ]`
block) to add an escape hatch:

```groovy
    if [ "\$ASM_VERDICT" = "small_fragmented" ] && [ ! -s "${other_gff}" ]; then
        echo "[WARN] ${out} is too small/fragmented for funannotate training (\${ASM_BP} bp, \${ASM_CTG} contigs, N50 \${ASM_N50}); skipping predict" >&2
        mkdir -p "${params.target}"
        [ -s "\$SKIP_REPORT" ] || printf 'out\tasmid\tlocustag\treason\ttotal_bp\tcontigs\tN50\n' > "\$SKIP_REPORT"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${out}" "${asmid}" "${locustag}" "preflight_small_fragmented" "\$ASM_BP" "\$ASM_CTG" "\$ASM_N50" >> "\$SKIP_REPORT"
        touch "\$PREDICTDIR/${out}.predict.skipped_too_small"
        touch ${out}.predict.done
        exit 0
    elif [ "\$ASM_VERDICT" = "small_fragmented" ]; then
        echo "[INFO] ${out} is small/fragmented but has Prodigal evidence (${other_gff}); proceeding instead of skipping" >&2
    fi
```

(Most microsporidia genomes gated by size alone — not fragmentation — never hit
this branch at all, since `ASM_VERDICT` requires both small AND fragmented; this
only matters for the subset that happen to also be fragmented. Known risk, not a
blocker: a genome fragmented enough to fail funannotate's own 30-training-model
requirement may still fail even with Prodigal evidence supplied, just later and
more expensively — `--min_training_models 30` is only required for AUGUSTUS/SNAP
training and those are being zeroed anyway on this branch, see Step 4.)

- [ ] **Step 4: Add the `--other_gff` flag, weight overrides, and the two missing validated-recipe flags** (replace lines 208-219)

```bash
    # Microsporidia Prodigal supplement (nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md):
    # other_gff non-empty means GENEMARK_RUN ran the Prodigal supplement for
    # this genome (is_microsporidia=true AND below predict_min_asm_bp).
    # AUGUSTUS/SNAP forced off here specifically -- both need >=200 BUSCO
    # training models these genomes don't have (Microsporidia_predict
    # PLAN.md 9.15) -- this does NOT change augustus/snap weighting for any
    # genome where other_gff is empty, which keeps today's default-on
    # behavior. --no-evm-partitions and --min_protlen 30 mirror the
    # validated recipe (microsporidia-default.json) for these near-zero-
    # intergenic compact genomes; both are new flags this pipeline didn't
    # pass before, applied ONLY on this branch.
    #
    # Why zeroing the weight is enough, and no further BUSCO change is
    # needed: funannotate's predict.py only populates RunModes["augustus"]/
    # ["snap"] when their StartWeight > 0; RunBusco is only set True when
    # some RunModes value == "busco". With augustus/snap/glimmerhmm all at
    # weight 0, RunBusco never becomes True and BUSCO is never invoked at
    # all -- not "invoked and tolerated despite too few models". The
    # --min_training_models 30 check a few lines below is itself nested
    # inside `if "augustus" in RunModes:`, so it stays unreachable dead code
    # on this branch. This was confirmed by tracing predict.py directly (not
    # just inferred) and matches the failure already documented in
    # Microsporidia_predict/STATUS.md:73-82 (passing --busco_db
    # microsporidia_odb10 with augustus/snap weight > 0 found only 36
    # complete BUSCO models against a >=200 requirement and failed; the fix
    # there was this same weight-zeroing, not any BUSCO-side workaround).
    # Do NOT also add --min_training_models 0 or drop --busco_db here --
    # neither is necessary and both would be pure noise.
    GENEMARK_GTF_FLAG=()
    OTHER_GFF_FLAG=()
    EXTRA_PREDICT_ARGS=()
    WEIGHT_ARGS=(codingquarry:0 glimmerhmm:0)
    if [ -s "${other_gff}" ]; then
        echo "[INFO] ${out}: using Prodigal evidence from ${other_gff} (--other_gff weight 5)"
        OTHER_GFF_FLAG=(--other_gff "${other_gff}:5")
        WEIGHT_ARGS+=(augustus:0 snap:0)
        # --busco_seed_species is REQUIRED here, not cosmetic: predict.py
        # unconditionally sets RunModes["augustus"]="busco" and requires a
        # pre-existing trained_species entry to seed it, REGARDLESS of
        # augustus/snap EVM weight (confirmed directly: RunModes["augustus"]
        # population happens before/independent of the weight-gated EVM
        # combiner step). Without a valid entry, predict hard-aborts with
        # "ERROR: --busco_seed_species {} is not valid" (confirmed by a
        # real, non-stub run against Ordospora colligata OC4). The literal
        # string "microsporidia" is NOT a valid seed species -- confirmed
        # neither this project nor the shared production funannotate_db has
        # ever registered an Augustus species by that name (only real
        # microsporidia strain entries like Encephalitozoon_cuniculi
        # exist). Use an existing real trained species as the seed instead
        # -- its own weight is 0 so its actual training output never
        # influences the final EVM gene models, it only exists to satisfy
        # this startup precondition.
        EXTRA_PREDICT_ARGS+=(--no-evm-partitions --min_protlen 30 --busco_seed_species Encephalitozoon_cuniculi)
    fi
    # -s not -n: GENEMARK_RUN's too-small-genome skip path emits a real but
    # deliberately empty ${out}.genemark.gtf (see GENEMARK_RUN/main.nf).
    if [ -s "${genemark_gtf}" ]; then
        echo "[INFO] ${out}: using pre-computed GeneMark GTF from ${genemark_gtf}"
        GENEMARK_GTF_FLAG=(--genemark_gtf "${genemark_gtf}")
        WEIGHT_ARGS+=(genemark:1)
    fi

    # other_gff/genemark_gtf are `val`, not `path`, in this process's input
    # tuple -- Nextflow never stages/symlinks them into THIS task's own
    # $PWD, so they're referenced by their original absolute path under
    # GENEMARK_RUN's own separate task work directory. Same missing-bind bug
    # class already documented above for $PWD/genome_input.fa: a path
    # outside SING_BINDS's explicit list is invisible inside the container
    # even though the host shell can read it fine. Confirmed empirically
    # (2026-09-04, real non-stub run against Ordospora colligata OC4):
    # without this, funannotate predict fails with "<path>/Ordospora_
    # colligata_OC4.other.gff3 is not a valid file, exiting" despite the
    # host-side `[ -s "${other_gff}" ]` check above having already confirmed
    # the file exists and is non-empty. Mirrors GENEMARK_RUN's own
    # training_bam dirname-binding pattern (GENEMARK_RUN/main.nf).
    if [ -n "${other_gff}" ]; then
        SING_BINDS="$SING_BINDS,$(dirname "${other_gff}"):$(dirname "${other_gff}")"
    fi
    if [ -n "${genemark_gtf}" ]; then
        SING_BINDS="$SING_BINDS,$(dirname "${genemark_gtf}"):$(dirname "${genemark_gtf}")"
    fi
    SING="apptainer exec ${SING_BINDS} ${params.funannotate_sif}"
```

- [ ] **Step 5: Pass the new flags into the `funannotate predict` invocation** (replace lines 221-228)

```groovy
    \$SING funannotate predict --name ${locustag} -i "\$GENOME_IN" --strain "${strain}" \\
        -o "\$PREDICTDIR" -s "${species}" --cpu ${task.cpus} --busco_db ${busco_lineage} \\
        --AUGUSTUS_CONFIG_PATH \$AUGUSTUS_CONFIG_PATH -w "\${WEIGHT_ARGS[@]}" \\
        --min_training_models 30 --tmpdir \$TMPDIR --SeqCenter ${params.seqcenter} \\
        --keep_no_stops --header_length ${header_length} --protein_evidence ${params.proteins} \\
        --max_intronlen ${params.max_intronlen} --min_intronlen ${params.min_intronlen} \\
        --tbl2asn "\$TBL2ASN_PARAMS" --table ${transl_table} --auto-skip-genemark \\
        "\${ABINITIO_REUSE_FLAG[@]}" "\${GENEMARK_GTF_FLAG[@]}" "\${OTHER_GFF_FLAG[@]}" "\${EXTRA_PREDICT_ARGS[@]}" || true
```

- [ ] **Step 6: Verify no other caller of `FUNANNOTATE_PREDICT` needs updating**

Run: `grep -rn "FUNANNOTATE_PREDICT(" nextflow/subworkflows/ nextflow/workflows/`
Expected: only `FUNANNOTATE_PREDICT(rep_with_gtf)`/`FUNANNOTATE_PREDICT(rep_and_indep_with_gtf)`/
`FUNANNOTATE_PREDICT_SIB(sibling_predict_with_gtf)` in `FUNANNOTATE_PREDICTION.nf` —
each already produces the correct 12-field tuple as of Task 3 Step 6, so no
further edits needed at these call sites.

- [ ] **Step 7: Commit**

```bash
git add nextflow/modules/funannotate/predict/FUNANNOTATE_PREDICT/main.nf nextflow/subworkflows/local/FUNANNOTATE_PREDICTION.nf
git commit -m "FUNANNOTATE_PREDICT: consume other_gff, apply microsporidia weight/flag overrides"
```

---

### Task 5: Test fixture + stub-run validation

**Files:**
- Create: `nextflow/tests/data/funannotate/test_samples_funannotate_microsporidia.csv`
- Create: `nextflow/tests/data/funannotate/source_microsporidia/GCF_TEST003/GCF_TEST003_genomic.fna.gz`
- Create: `nextflow/conf/test_funannotate_microsporidia.config`

- [ ] **Step 1: Add a Microsporidia row to a new fixture CSV**

```csv
ASMID,SPECIES_IN,STRAIN,BIOPROJECT,NCBI_TAXONID,BUSCO_LINEAGE,PHYLUM,SUBPHYLUM,CLASS,SUBCLASS,ORDER,FAMILY,GENUS,SPECIES,LOCUSTAG,TRANSL_TABLE
GCF_TEST001,Testus fungus STRAIN1,STRAIN1,PRJNA000001,5085,dikarya,Ascomycota,Pezizomycotina,Dothideomycetes,,,,Testus,Testus fungus,FF5840CF,1
GCF_TEST002,Secondus fungus STRAIN2,STRAIN2,PRJNA000001,5085,dikarya,Ascomycota,Pezizomycotina,Dothideomycetes,,,,Secondus,Secondus fungus,F2EE6837,1
GCF_TEST003,Testospora microsporidia STRAIN1,STRAIN1,PRJNA000003,6033,microsporidia_odb10,Microsporidia,,,,,,Testospora,Testospora microsporidia,F00D0003,1
```

(`SPECIES` for row 3 is `Testospora microsporidia`, `STRAIN` is `STRAIN1` —
`makeSampleTag` produces `Testospora_microsporidia_STRAIN1`, and
`loadMicrosporidiaOutSet()` must produce the identical string from the same
`SPECIES`/`STRAIN` columns — this is the exact bug class Task 1's design note
calls out.)

- [ ] **Step 2: Stage a tiny synthetic genome for GCF_TEST003**

Match whatever synthetic-FASTA shape `tests/data/funannotate/source/GCF_TEST001/`
already uses. Since this is a `-stub-run`, real sequence content is irrelevant to
`GENEMARK_RUN`'s stub (Task 3 Step 4) — only the file's existence matters for the
channel's `gz.exists()` filter in `workflows/funannotate.nf`.

- [ ] **Step 3: Write `conf/test_funannotate_microsporidia.config`**

Copy `conf/test_funannotate.config`, changing only:

```groovy
params {
    samples  = "${projectDir}/tests/data/funannotate/test_samples_funannotate_microsporidia.csv"
    source   = "${projectDir}/tests/data/funannotate/source_microsporidia"
    target           = "${projectDir}/tests/output/funannotate_genome_annotation_microsporidia"
    training_target  = "${projectDir}/tests/output/funannotate_genome_annotation_training_microsporidia"

    run_sra_fetch  = false
    run_ani_reuse  = false
}
```
(same `process { }` executor/resource overrides as `test_funannotate.config` — no
new `withName` entries needed since this plan added no new process.)

- [ ] **Step 4: Run the stub and check what actually happens under `-stub-run`**

Run:
```bash
nextflow run nextflow/main.nf -c nextflow/nextflow.config \
    -c nextflow/conf/test_funannotate_microsporidia.config \
    -profile funannotate --pipeline funannotate -stub-run
```
Expected: exits 0. **Note the honest limitation**: `GENEMARK_RUN`'s stub (Task 3
Step 4) always touches an empty `${out}.other.gff3` — it does not execute the real
bash gate (`is_microsporidia` + `$ASM_BP` comparison), so this stub-run validates
only that the DAG wiring is intact (tuple arities match at every call site, no
missing-channel errors, `GENEMARK_RUN`/`FUNANNOTATE_PREDICT` both run for all 3
rows) — **not** that GCF_TEST003 actually gets real Prodigal evidence. That
requires a real (non-stub) run against a real small microsporidia genome, done in
Step 6 below, as a follow-up validation step before Task 6's real-genome rollout.

- [ ] **Step 5: Run the existing (non-microsporidia) stub-run to confirm zero regression**

Run:
```bash
nextflow run nextflow/main.nf -c nextflow/nextflow.config \
    -c nextflow/conf/test_funannotate.config \
    -profile funannotate --pipeline funannotate -stub-run
```
Expected: exits 0, identical DAG/behavior to before this plan — no
`prodigal_selected_microsporidia.tsv` created (since `microsporidiaSet` is empty
for this fixture's 2 non-microsporidia rows). Capture a task-count/task-name log
from a run of this exact command taken *before* starting Task 1, to have a
concrete baseline to diff this run against.

- [ ] **Step 6: Run a real single-genome validation against one held-out microsporidia genome**

Before the full corpus run, validate against one real small microsporidia genome
end-to-end (not stub) — e.g. reuse
`input_clean_genomes/GCF_000803265.1_ASM80326v1.masked.fasta.gz` (Ordospora
colligata OC4) via a one-off `samples.csv` row with `--asmid` filtering
(`params.asmid`, already supported per `workflows/funannotate.nf`'s `asmidFilter`).
Confirm: `prodigal_selected_microsporidia.tsv` gets a row for this genome, the
predict log shows both `--genemark_gtf` and `--other_gff` flags, and the resulting
gene count is in the vicinity of the previously-validated Sn/Sp table (report as a
before/after table per this project's collaboration norms — do not claim
equivalence without the actual numbers).

- [ ] **Step 7: Commit**

```bash
git add nextflow/tests/data/funannotate/test_samples_funannotate_microsporidia.csv nextflow/tests/data/funannotate/source_microsporidia nextflow/conf/test_funannotate_microsporidia.config
git commit -m "Add microsporidia Prodigal supplement stub-run test fixture"
```

---

### Task 6: Real-genome rollout (after Task 5 Step 6 passes)

Not a code task. Once the single-genome validation (Task 5 Step 6) confirms
correct gene counts/Sn-Sp in the expected range:

1. Run against a small batch (5-10 genomes spanning the size range found in the
   corpus — recall 4/25 sampled genomes were `small_fragmented` and 21/25 were
   `ok`-but-below-8Mb, so most gated genomes will hit the size-only gate, not the
   fragmented one) before the full 127-species run.
2. Report Sn/Sp per this batch, before/after, per this project's own
   collaboration norms (no pooled-only numbers, no prose-only claims).
3. Only then point at the full corpus.
