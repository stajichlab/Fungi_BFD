---
topic: intraspecific-ani-diversity
description: Within-species ANI distributions across strain collections, and what they imply for treating a SPECIES label as a proxy for "safe to pool" (ab-initio parameter reuse, strain-level analyses in general).
created: 2026-07-24
last_updated: 2026-07-24
status: active
---

# Intraspecific ANI Diversity

> Opened 2026-07-24 while backfilling `_shared_abinitio/` stores for
> `todo/species_level_abinitio_reuse.md` (T-004/T-006). Tracks within-species ANI
> distributions across strain collections in this dataset — relevant beyond the
> ab-initio-reuse use case to anything that assumes a shared SPECIES label implies
> low genomic divergence.

## Findings

### F-005 — Beauveria bassiana's within-species ANI distribution is consistent with an undifferentiated species complex, not simple strain-level polymorphism

**Status:** preliminary (2026-07-24, n=223 strains vs. one representative)

**Claim:** Computed via `nextflow/bin/species_reuse_clusters.py --species "Beauveria
bassiana"` against the production `results/ANI/skani/GENUS/ani.db` (206,885
within-species pairs loaded across the dataset). Against the highest-BUSCO
ANI-covered representative (`Beauveria_bassiana_KW2`), the 223 other *B. bassiana*
strains show: **median ANI 97.16%, p25 96.99%, p75 97.38%, min 90.89%, max 99.92%**.
Only 32/223 (14%) strains clear a 99.0% ANI threshold to the representative; 190/223
(85%) sit in the 90.89–98.97% range — well below what's typically expected for
conspecific fungal strains.

Contrast: the same computation for *Aspergillus fumigatus* (n=375 vs. representative
`Aspergillus_fumigatus_Z5`) gives **median ANI 99.62%**, with 373/375 (99.5%) strains
clearing the same 99.0% threshold. *A. fumigatus* behaves like a tight, low-diversity
species under this measure; *B. bassiana* does not.

**Why it matters:** *Beauveria bassiana sensu lato* is independently documented in the
literature as a species complex with substantial cryptic diversity not fully resolved
by current taxonomy — this ANI distribution is consistent with that, not with a
labeling or pipeline error. It directly validates the ab-initio-reuse plan's design
decision to gate parameter sharing on measured ANI rather than trusting the SPECIES
column alone (`todo/species_level_abinitio_reuse.md` S2, S6 decision 1): a naive
"pool everyone with the same SPECIES label" policy would have wrongly shared one
strain's trained AUGUSTUS/SNAP/GeneMark-ES parameters across a genuinely
heterogeneous set for 190/223 strains here.

**Implications beyond ab-initio reuse:** any future analysis in this dataset that
pools strains by SPECIES label (comparative genomics, pangenome construction,
phylogenetic sampling assumptions, etc.) should check the within-species ANI
distribution first for genera/species known or suspected to be complexes — *Beauveria*
is a concrete example already in hand; there are likely others in this ~23,000-genome
collection (worth a follow-up scan: `SELECT query_species, MEDIAN ANI ... GROUP BY
query_species` style query against `ani.db`, flagging species whose median
within-species ANI falls well below the ~99% seen for well-resolved species like
*A. fumigatus*).

**Evidence:** `genome_annotation/_reuse_assignments/abinitio_reuse_assignments.Beauveria_bassiana.csv`
(224 rows, ANI-to-representative per strain); `results/ANI/skani/GENUS/ani.db`
(production, 385,795 pairs / 204 genera as of 2026-07-24).

**Not yet done:** a systematic scan across all species in `ani.db` for similarly wide
within-species ANI distributions (this finding is from investigating one species while
implementing T-006, not a deliberate survey).
