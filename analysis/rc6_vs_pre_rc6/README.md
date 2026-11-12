# rc.6 vs pre_rc6 gene-model comparison, and scoring against a reference annotation

`compare_gff_pairs.py` compares funannotate predictions in `Fungi_BFD_runs/genome_annotation.pre_rc6/`
(old) and `Fungi_BFD_runs/genome_annotation/` (rc.6, new). Run it with `python3.12`. It needs
`gffcompare` on PATH (`module load gffcompare/0.12.10`). It reads plain, `.gz` and `.zst` files.

| Subcommand | What it does |
|---|---|
| `pairs` | Finds genomes with a GFF3 in both trees (key `out` from `samples.csv`, same rule as `makeSampleTag`). Optional phylum-stratified pilot sample. |
| `compare` | Per genome: gene, mRNA, exon, CDS and intron counts for old and new, the same counts on scaffolds that are identical in both FASTA files, tRNA counts, and gffcompare sensitivity and precision (old = reference). `--shard I/N` splits the work. |
| `summarize` | Per-phylum medians and the distribution of the new/old gene-count ratio. |
| `reference` | Scores old and new against a reference GFF3 (for example NCBI RefSeq) in two modes: `exon` (exons as annotated, UTRs included) and `cds` (CDS segments as exons, the fair mode when the prediction has no UTRs). Optional `--introns-dir` (one per-intron outcome table per genome and set) and `--save-dir` (raw gffcompare `.tmap`, `.refmap`, `.stats`). |

Only protein-coding genes (genes with an mRNA child) are scored. tRNA genes are counted separately.

## Example

```bash
module load gffcompare/0.12.10
S=Fungi_BFD/analysis/rc6_vs_pre_rc6/compare_gff_pairs.py
python3.12 $S pairs --samples Fungi_BFD_runs/samples.csv --runs-root Fungi_BFD_runs --out pairs_all.tsv
python3.12 $S compare --pairs pairs_all.tsv --shard 0/2 --out compare_shard0of2.tsv --procs 8
python3.12 $S summarize --compare 'compare_shard*of2.tsv' --out-dir summary
# table columns: out, ref_gff, old_gff, new_gff
python3.12 $S reference --table reference_pairs.tsv --out reference_compare.tsv --procs 4 \
    --introns-dir reference_introns --save-dir reference_gffcompare
```

Outputs of the 2026-11-11 runs are in `Fungi_BFD_runs/analysis/rc6_vs_pre_rc6/`.

## Things to know

- Scaffolds are "shared" when name and length match in both FASTA files. If a FASTA is missing, the
  seqids with features in both GFFs are used (`scaffold_source` says which).
- gffcompare calls its second number precision, not specificity.
- In the gffcompare stats file, "Missed introns" and "Novel introns" are not 100 minus Sn or Pr.
  Use the status column of the per-intron table (`reference_introns`), which agrees with Sn and Pr.
  What the stats lines count was not established.
- The per-intron `phase` is the number of coding bases 5' of the intron, mod 3. The CDS phase column
  of the first segment is not used, so models with a partial 5' end can be off.
- Reference scores are agreement with the reference, not accuracy. RefSeq annotations carry
  alternative isoforms while predictions have one mRNA per gene, so transcript-level sensitivity is
  capped; the locus level is the better figure. Protein evidence (`swissprot_fungi.faa`) probably
  contains proteins of the model species (not checked), which can inflate their scores.
- `make_sample_tag()` copies `nextflow/bin/pick_representative_strain.py`. Keep them in sync.
