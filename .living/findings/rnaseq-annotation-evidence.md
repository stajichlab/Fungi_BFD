---
topic: rnaseq-annotation-evidence
description: Availability and quality of RNA-seq-derived transcript evidence (Trinity assemblies) used to support genome annotation.
created: 2026-06-18
last_updated: 2026-06-18
status: active
---

# RNA-seq Transcript Evidence for Annotation

## F-001: The transcript-evidence gap is dominated by missing assemblies, not poor-quality ones
**Status:** preliminary
**Claim:** Across 1178 fungal species for which genome-guided Trinity assembly
was attempted, only 239 (~20%) currently have a transcript FASTA on disk; 939
are no-assembly failures (874 MISSING + 65 EMPTY). Among the 239 that *do*
exist, assembly quality is overwhelmingly adequate: 234 PASS / 3 BORDERLINE / 2
FAIL using a NUM_TRANSCRIPTS tier (FAIL<100, PASS≥1000), and this split is
stable across swept thresholds. The annotation-evidence bottleneck is therefore
*production* (getting an assembly to exist) rather than *quality filtering*.
**Implications:** Effort to improve RNA-seq-supported annotation should target
the 939 missing assemblies (read fetch / re-assembly, cross-ref
`rnaseq_rerun_report.tsv`) rather than re-filtering the existing 239. Only 2
existing assemblies (Alternaria_alstroemeriae=12, Zychaea_mexicana=36 transcripts)
should be treated as no-evidence despite having a FASTA.
**Tags:** rnaseq, trinity, annotation, funannotate, qc, coverage

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-06-18 | analysis/rnaseq_trinity_qc | TrinityGG_summary.tsv (239) + TrinityGG_failed.tsv (944) | Fungi_BFD | 239/1178 have FASTA; 234 PASS / 3 BORDERLINE / 2 FAIL; robust to threshold sweep | supports |

### Open Questions
- Of the 939 missing assemblies, how many have usable SRA accessions available
  (vs. genuinely lacking public RNA-seq)? Cross-ref `samples.rnaseq_sra.csv` and
  `rnaseq_rerun_report.tsv`.
- Should transcript N50 (contiguity), not just count, gate usability? Two
  BORDERLINE/FAIL species have moderate counts but low N50.
