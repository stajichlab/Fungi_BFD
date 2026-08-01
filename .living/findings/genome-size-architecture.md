---
topic: genome-size-architecture
description: Determinants and cross-taxonomic patterns of genome size and genomic composition (TEs/repeats, codon/AA usage, gene structure, domain repertoire) in fungi.
created: 2026-06-18
last_updated: 2026-06-18
status: active
---

# Genome Size & Genomic Architecture

> Topic opened 2026-06-18 alongside the genome-size/composition ideation session
> (`analysis/ideas/2026-06-18-genome-size-composition-framework/`). No empirical
> findings recorded yet — entries will accrue as the proposed analyses run.
> Candidate first analyses: variance partitioning of size across taxonomy (#1/#13),
> order-parameter test of TE expansion (#8), Zipf lexicon shifts (#17).

## Data-state notes (not yet findings)

- **Corpus:** `tables/asm_stats.tsv.gz` — 22,412 assemblies / 8,061 species, with
  `total_length_bp` (size), `masked_pct` (TE/repeat proxy), `gc_pct`, contiguity,
  telomere completeness. Mean size ≈ 35 Mb.
- **Quality tails flagged:** genome size spans ~0 Mb to ~1711 Mb — both tails are
  likely data-quality artifacts (partial/contaminant assemblies; extreme outliers),
  not biology. Any size analysis must filter/inspect these and control assembly
  quality (contig_count/N50/BUSCO).
- **Prerequisite:** the kingdom-wide `All_Taxa` compositional merge (aa_freq,
  codon_freq, intergenic, pfam/cazy/merops, idp, etc.) is **not yet materialized**;
  most composition-vs-size analyses need one BFD.nf merge run first.

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| (none yet) | | | | | |

### Open Questions
- How much genome-size variance sits at each taxonomic level (species→phylum)?
- Is `masked_pct`→size linear, threshold-like, or regime-dependent (TE runaway)?
- Which indicators *cause* size vs. correlate vs. are quality artifacts?
