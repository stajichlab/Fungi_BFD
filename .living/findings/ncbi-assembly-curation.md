# Finding topic: NCBI assembly curation (dedup & sanitization)

## F-002 — Reference-genome dedup tracks annotation quality, not GCF/GCA prefix or recency

**Status**: established (empirical, full-dataset)
**Date**: 2026-06-19
**Tags**: ncbi, dedup, gcf-gca, samples-csv, curation

### Claim
When the BFD fungal genome set is deduplicated to one assembly per
species+strain isolate, the genome that was kept is the best-annotated /
reference assembly — which is **not** predicted by a simple "prefer RefSeq GCF"
or "newest accession" rule.

### Evidence
Comparing the full pre-curation set (`misc/samples.csv.prechange`, 22,412 rows)
against the hand-curated `samples.csv`:

- Of **39** isolates present as **both** a GCF and a GCA assembly:
  - GCA kept / GCF dropped: **20**
  - GCF kept / GCA dropped: **16**
  - both kept (intentional): **3**
  → a blanket "prefer GCF" rule would reverse 20 human decisions.
- "Keep newest accession" reproduces only **71 / 110** clean single-pick
  dedup groups (39 mismatches).
- Mismatches consistently kept curated/annotated assemblies (JGI MycoCosm
  genomes `Babin1`, `Cybja1`, `Lipst1_1`, …; classic references H99 `CNA3`,
  *P. oryzae* 70-15) over newer un-annotated re-sequences.

### Implications
- Encode dedup as an explicit curated exclude list
  (`data/curation/exclude_asmids.txt`), not a mechanical rule.
- Apply an automatic tie-breaker (prefer GCF, else newest GCA) **only** to new
  uncurated collisions; protect intentional multi-keeps via
  `data/curation/keep_dupes.csv` (6 isolates: PH-1, S288C R64, ME14, …).

### Provenance
Analysis in this session via `scripts/sample_sanitize.py` keys; see
`.living/decisions.md` (samples.csv reproducibility) and
`data/curation/README.md`.
