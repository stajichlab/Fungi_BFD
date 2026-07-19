# Strategy: Placing Taxonomically-Unassigned ("Orphan") Genomes via ANI

## Problem

`samples.csv` has 22,929 genomes; **329 have no GENUS assigned**:

| Bucket | Count | Deepest known rank |
|---|---|---|
| A | 152 | FAMILY assigned |
| B | 121 | CLASS only |
| C | 38  | PHYLUM only |
| D | 18  | nothing |

Goal: for buckets A/B, search each orphan against sketches of named relatives
in its known FAMILY/CLASS to see if ANI reveals a close match that suggests a
genus. For orphans with no useful relatives (bucket C/D), check whether they
cluster tightly with *each other* — a sign of an undescribed shared lineage.

## Constraint that shaped the approach

`nextflow/compare_ANI.nf` groups genomes by a single taxonomic rank
(`--compare RANK`) and computes a full **all-vs-all** ANI triangle within each
group. That's fine for bucket A: family sizes for the affected families range
from 24 to 2268 members, well within what `skani triangle` is designed for.

It breaks down for bucket B: the affected classes are enormous —
Sordariomycetes (5827 members), Eurotiomycetes (3123), Dothideomycetes (2670),
Saccharomycetes (2621) — and triangulating an entire class just to place a
handful of orphans (25, 16, 16, 26 respectively) is disproportionate: cost is
O(N²) when only O(orphans × N) is actually needed.

## Solution

**Phase 1 — Bucket A (has FAMILY): use `compare_ANI.nf` unmodified.**

A row is only dropped from a `--compare FAMILY` group if FAMILY itself is
empty, so orphans with FAMILY assigned are already included in the same
triangle as their named relatives — no code change needed.

```bash
nextflow run nextflow/compare_ANI.nf -c nextflow/nextflow.config -profile ani \
    -params-file params_clean_genomes_ani.yaml \
    --compare FAMILY --taxon FAMILY:Boletaceae -resume
# repeat --taxon per affected family (~20 families; see counts below)
```

Read the result from `results/ANI/skani/FAMILY/<Family>/<Family>_ANI_report.txt`
— each orphan will show up either inside a named cluster (≥95% ANI) or in the
outliers list (best ANI <90%) with its best match named.

**Phase 2 — Bucket B/C (CLASS or PHYLUM only): new `query_ANI.nf` workflow.**

Added `nextflow/query_ANI.nf`, a companion to `compare_ANI.nf` for this
specific case. Instead of a full triangle, it splits each `--compare` group
into **query** genomes (missing `--query_rank`, default GENUS) and
**reference** genomes (everyone else in the group), then runs
`skani dist --ql <queries> --rl <references>` — cost is O(queries × refs)
instead of O(refs²). Groups with zero orphans are skipped automatically, so no
`--taxon` restriction is needed to scope a run. Sketches are cached in the
same `work/ANI/sketch_cache/` as `compare_ANI.nf`, so a genome sketched by
either workflow is reused by the other.

```bash
nextflow run nextflow/query_ANI.nf -c nextflow/nextflow.config -profile ani \
    --samples samples.csv \
    --genome_name_style asmid --genome_dir input_clean_genomes --genome_suffix .fa.gz \
    --compare CLASS --query_rank GENUS -resume
```

Output per group: `results/ANI/skani_query/CLASS/<Class>/<Class>_query_report.txt`
(human-readable, one block per orphan with its top reference hits) and
`<Class>_query_calls.tsv`. All groups' calls are merged into
`results/ANI/skani_query/CLASS/orphan_calls.csv`, sorted strongest match first.

Also use `query_ANI.nf` for orphans that have FAMILY but you additionally want
checked against the full CLASS (not just their family) by passing
`--compare CLASS --query_rank FAMILY`.

**Phase 3 — Orphan-vs-orphan clustering (buckets C/D, or all 329 as a sanity
check): `compare_ANI.nf` again, on a reduced sample sheet.**

For orphans with no useful named relatives, build a samples CSV containing
only the orphan rows (any grouping column works — even a constant sentinel —
since there's nothing to compare them against except each other) and run the
normal triangle:

```bash
awk -F',' 'NR==1 || $13==""' samples.csv > samples.orphans_only.csv
nextflow run nextflow/compare_ANI.nf -c nextflow/nextflow.config -profile ani \
    -params-file params_clean_genomes_ani.yaml \
    --samples samples.orphans_only.csv --compare PHYLUM -resume
```

Tight orphan-orphan clusters (≥95% ANI) with no named relative in Phase 1/2
are candidates for an undescribed shared lineage.

## Classification rule (Phase 4)

Reuses the existing `--ani_cluster_threshold` (95%) / `--ani_outlier_threshold`
(90%) semantics from `report_ani.py`, and is built directly into
`report_query_ani.py`'s per-orphan tiers:

| Tier | Condition | Action |
|---|---|---|
| `same_genus_high_confidence` | best ANI ≥ 95% to one named genus | propose adopting that genus |
| `closely_related_review` | 90% ≤ best ANI < 95% | flag "closely related to Genus X, possibly novel species" for manual review |
| `no_close_match` | best ANI < 90% (but some alignment) | leave unresolved |
| `no_alignment` | nothing passed the aligned-fraction floor | leave unresolved; isolate needing formal taxonomic description |
| *(orphan-orphan, Phase 3)* | tight cluster, no named relative ≥90% | candidate undescribed genus-level clade — flag with a placeholder tag, don't assign a real genus |

## New files added

- `nextflow/query_ANI.nf` — asymmetric orphan-vs-reference ANI search (skani dist).
- `nextflow/bin/report_query_ani.py` — per-group tiered classification report + `_query_calls.tsv`.
- `nextflow/bin/combine_query_calls.py` — merges all groups' calls into one `orphan_calls.csv`.

## Validation performed

- `query_ANI.nf` stub-run against the real Boletaceae subset correctly split the
  407-genome family into 22 query genomes vs the rest as references — exactly
  matching the known orphan count — with correct per-role sketch batching.
- `report_query_ani.py` / `combine_query_calls.py` unit-tested against synthetic
  ANI pairs covering all four tiers; classification and sort order confirmed.

## Not yet done

Real (non-stub) runs haven't been launched — both phases involve real cluster
time (sketching + comparing potentially thousands of genomes for the larger
classes). Launch Phase 1 and Phase 2 via `sbatch run_ANI.sh` (or an equivalent
launcher for `query_ANI.nf`) when ready.
