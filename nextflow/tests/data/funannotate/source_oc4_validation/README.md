# source_oc4_validation — real-genome fixture (not tracked in git)

This directory backs `nextflow/conf/test_funannotate_oc4_validation.config`, a
real (non-synthetic) validation run against *Ordospora colligata* OC4. Unlike
the synthetic fixtures under `nextflow/tests/data/funannotate/source/`, this
is an actual NCBI RefSeq assembly, so it is excluded from the repo
(`nextflow/.gitignore`'s `tests/` rule) rather than committed.

## Regenerate

```bash
mkdir -p GCF_000803265.1_ASM80326v1
curl -o GCF_000803265.1_ASM80326v1/GCF_000803265.1_ASM80326v1_genomic.fna.gz \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/803/265/GCF_000803265.1_ASM80326v1/GCF_000803265.1_ASM80326v1_genomic.fna.gz"
```

Accession: GCF_000803265.1 (ASM80326v1), *Ordospora colligata* OC4,
BioProject PRJNA210314. See `samples_oc4_validation.csv` for the sample-sheet
metadata and `nextflow/docs/MICROSPORIDIA_PRODIGAL_BRANCH_PLAN.md` for why
this genome was chosen (small-but-complete: 2.29 Mb, 15 contigs, N50
228,601 — trips the size gate but not the fragmentation gate).
