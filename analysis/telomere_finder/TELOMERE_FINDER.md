# Telomere Finder for BFD Assemblies

## Goal
Add a flexible, reusable telomere-feature detection step to the BFD Nextflow
pipeline. The step should improve on the reference `find_telomeres.py` by
supporting multiple fungal monomer patterns, exact and fuzzy matching, per-tract
repeat count and total telomeric length, and ~500 bp of inward flanking
sequence.

## Inputs
- Per-genome genome FASTA (`.scaffolds.fa` or raw `<asmid>.fa.gz`).
- Fungal telomere monomers, including canonical `TTAGGG`/`CCCTAA` plus
  TeloBase-derived fungal variants.

## Outputs
### Per-genome
- `results/genome_stats/telomeres/<bucket>/<asmid>.telomeres.tsv.gz`
  - `scaffold`, `end` (5prime/3prime), `strand`, `monomer`
  - `repeat_count`, `tract_length` (total telomeric length for the tract)
  - `start`, `end_coord`, `terminal`
  - `tract_seq`, `flank_seq` (~500 bp inward)

### Merged table
- `tables/telomeres.parquet` — per-genome aggregates:
  - `telomere_scaffolds`, `telomere_tracts`
  - `telomere_total_length_bp`, `telomere_total_repeats`
  - `telomere_5prime_count`, `telomere_3prime_count`
  - `telomere_plus_count`, `telomere_minus_count`
  - `telomere_terminal_count`, `telomere_internal_count`
  - `telomere_top_monomer`, `telomere_monomers`

## Implementation
- `nextflow/bin/find_telomeres.py` — standalone Python finder.
  - Exact mode: regex-based tandem-repeat scan in terminal windows.
  - Fuzzy mode: `finditer_approx` using Biopython-like approximate string
    matching to tolerate mismatches/indels.
  - Reverse-complements patterns for canonical 5'/3' orientation;
    `--both-ends` searches all combinations.
  - IUPAC ambiguities expanded in regex mode.
- `nextflow/bin/summarize_telomeres.py` — aggregate per-genome TSVs into the
  merged table.
- `nextflow/modules/BFD/FIND_TELOMERES/main.nf` — per-genome process.
- `nextflow/modules/BFD/MERGE_TELOMERES/main.nf` — merge process.
- Wired into `nextflow/subworkflows/local/BFD_GENOME_STATS.nf`,
  `nextflow/subworkflows/local/BFD_MERGE.nf`, and `nextflow/workflows/BFD.nf`.
- Parameters added to `nextflow/conf/profile_BFD.config` and
  `nextflow/nextflow_schema.json`.

## Test
Initial validation on *Neurospora crassa* OR74A (ASMID
`GCA_009805915.1_ASM980591v1`) only found 2 terminal tracts, both at the 5'
end — this later turned out to be a symptom of a bug, not confirmation of
correctness (see Bug fix below).

- Nextflow `-stub-run` wiring check passed.
- Local `-process.executor=local` run of `FIND_TELOMERES` + `CALC_ASM_STATS`
  succeeded and published to `results/genome_stats/telomeres/`.
- `summarize_telomeres.py` → DuckDB Parquet conversion verified.

## Bug fix (2026-08-04)
Review + testing on 5 real genomes (*N. crassa* OR74A, 2x *A. nidulans*, 2x
*F. graminearum* incl. the chromosome-level `GCA_900073075.1_CML3066.v2`)
found that 3'-end detection was structurally broken and the `terminal` check
had zero tolerance, together dropping the large majority of real telomeres.
Fix plan was independently reviewed by a bioinformatics-expert agent and by
the Fable model before implementation. Full diagnosis, review notes, and fix
details in `.living/learnings.md` and `.living/decisions.md` ("H", 2026-08-04).
Post-fix, `GCA_900073075.1_CML3066.v2` finds 7/8 real telomeric ends (up from
1/8); *N. crassa* 1→13 hits; *A. nidulans* 1→9 hits.

**Known remaining gap (not fixed, tracked as follow-up)**: the "canonical
monomer at 5', RC at 3'" convention assumes a fixed strand orientation per
chromosome, which is actually arbitrary per assembly/contig. Only the
`TTAGGG`/`CCCTAA` pair in the default pattern set is RC-symmetric; other
default monomers have no RC twin, so a contig assembled in the "wrong" global
orientation can still miss real telomeres at both ends unless `--both-ends`
is passed. `--both-ends` defaults to `false` pipeline-wide.

## Reproduce
```bash
bash analysis/telomere_finder/run.sh
```

## Decisions
- Report `tract_length` as the actual matched span in bp and `repeat_count` as
  the inferred number of monomer copies; both are summed in the merged table.
- Include `tract_seq` and `flank_seq` by default in per-genome outputs because
  the user requested ~500 bp inward flanking sequence for manual inspection.
- Keep the merged table as aggregates only to avoid bloating the Parquet file.
- Store per-genome outputs in the existing hash-bucket layout under
  `genome_stats/telomeres/` for consistency with `asm_stats`, `BUSCO_genome`,
  etc.

## Tags
- telomeres
- genome-statistics
- nextflow
- bfd-pipeline
- fungi
