<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-06-25

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 0 sections | 2026-06-18 | — |
| decisions.md | 3 entries | 2026-06-25 | Persona-based ideation over BFD genome-size/composition data (Broader-9), Trinity-GG usability tiers keyed on NUM_TRANSCRIPTS, not N50, Initialized mycelium with non-destructive scaffold, not restructure |
| learnings.md | 6 entries | 2026-06-25 | MERGE_* steps: pass a manifest file, not thousands of staged inputs; bake mtime+size in for staleness, Nextflow `combine` spreads list-valued left items — use a scalar barrier, Any BFD genome-size analysis must control two confounds + needs the All_Taxa merge, scilintr needs a python3.12 venv here; flags "old" inside "threshold", TrinityGG_failed.tsv LOW_COUNT rows double-list successes, not failures |
| log/ | 15 sessions | 2026-06-25 | fungi-bfd (15) |
| findings/ | 2 findings across 4 topics | 2026-06-25 | annotation-failure-modes, ncbi-assembly-curation, genome-size-architecture, rnaseq-annotation-evidence |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-06-25 (heuristic)

## Tag clusters

- **mycelium** (3 entries) — L-4, L-6, D-3
- **BFD.nf** (2 entries) — L-1, L-2
- **caching** (2 entries) — L-1, L-2
- **comparative-genomics** (2 entries) — L-3, D-1
- **genome-size** (2 entries) — L-3, D-1
- **nextflow** (2 entries) — L-1, L-2

## Most recent (10)

- [2026-06-25] L-1: MERGE_* steps: pass a manifest file, not thousands of staged inputs; bake mtime+size in for staleness
- [2026-06-25] L-2: Nextflow `combine` spreads list-valued left items — use a scalar barrier
- [2026-06-19] L-3: Any BFD genome-size analysis must control two confounds + needs the All_Taxa merge
- [2026-06-18] L-4: scilintr needs a python3.12 venv here; flags "old" inside "threshold"
- [2026-06-18] L-5: TrinityGG_failed.tsv LOW_COUNT rows double-list successes, not failures
- [2026-06-18] L-6: Mycelium scripts need python3.12, not the default python3
- [2026-06-18] D-1: Persona-based ideation over BFD genome-size/composition data (Broader-9)
- [2026-06-18] D-2: Trinity-GG usability tiers keyed on NUM_TRANSCRIPTS, not N50
- [2026-06-18] D-3: Initialized mycelium with non-destructive scaffold, not restructure

## By tag

- `mycelium`: L-4, L-6, D-3
- `BFD.nf`: L-1, L-2
- `caching`: L-1, L-2
- `comparative-genomics`: L-3, D-1
- `genome-size`: L-3, D-1
- `nextflow`: L-1, L-2
- `python-version`: L-4, L-6
- `qc`: L-5, D-2
- `rnaseq`: L-5, D-2
- `tooling`: L-4, L-6
- `trinity`: L-5, D-2
- `arg-max`: L-1
- `asm-stats`: L-1
- `asm_stats`: L-3
- `assembly-quality`: L-3
- `barrier-channel`: L-2
- `channel-semantics`: L-2
- `collectFile`: L-1
- `combine-operator`: L-2
- `confounds`: L-3
- `data-reconciliation`: L-5
- `denominator`: L-5
- `double-counting`: L-5
- `environment`: L-6
- `false-positive`: L-4
- `hooks`: L-6
- `hpcc`: L-6
- `idea-generator`: D-1
- `ideation`: D-1
- `init`: D-3
- `linting`: L-4
- `manifest`: L-1
- `merge`: L-1
- `non-destructive`: D-3
- `phylogenetic-non-independence`: L-3
- `preview-lint`: L-2
- `repo-structure`: D-3
- `resume`: L-1
- `scilintr`: L-4
- `sensitivity-analysis`: D-2
- `staging`: L-1
- `staleness`: L-1
- `subagents`: D-1
- `thresholds`: D-2

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 skills/core/scripts/recall_lessons.py --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
