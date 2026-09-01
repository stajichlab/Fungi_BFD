# INVESTIGATE TELOMERES

## Purpose

Identify how telomeric repeats vary across fungi, are there phyletic patterns? Are there correlations with length of telomere and repeat type? Are there correlations to specific genes / proteins in the genome and telomere length (eg telomerase) or other genes.

## Status

**Status**: draft

## Datasets

- 'tables/telomeres.parquet' - will have all the telomere data loaded for all genomes we could predict them in. Develop a data mining approach to look at how the number of repeats, the length of repeats, sequence composition all vary and if these have taxonomic or phylogenetic patterns. Examine GC content differences between telomeres and rest of genome. Start to look at pfam domain composition to see if specific domains are present (probably) but perhaps dig deep into know telomere sequence interacting proteins and see how these are evolving, their length, maybe a structure prediction of the protein?

## Algorithms

<!-- Reference entries from algorithms/ALGORITHM_MANIFEST.md -->
- `[algorithm-name]` — [brief description of how it's applied]

## Parent Analysis

<!-- If this builds on prior work, reference it. Otherwise delete this section. -->
- Parent: `[parent-analysis-name]`
- What's different: [what this analysis adds or changes]

## Key Findings

<!-- Update as work progresses. Use bullet points. -->
- [Finding 1]
- [Finding 2]

## Open Questions

<!-- What remains unresolved or needs follow-up? -->
- [Question 1]
- [Question 2]

## Reproducibility

To reproduce all outputs:

```bash
cd analysis/[analysis-name]
bash run.sh
```

## Outputs

| File | Description |
|------|-------------|
| `outputs/[filename]` | [description] |
