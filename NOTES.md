# Notes

## Issues

### 1. Improperly coded paired-end data — separate SRA submissions

These two datasets had paired-end reads split across separate SRA accessions instead of a single paired submission:

```
Pseudocercospora_eumusae,321146,SRR2093558,26033251,ILLUMINA
Pseudocercospora_eumusae,321146,SRR2093556,26033251,ILLUMINA
```

### 2. Same problem

```
Pseudocercospora_musae,113226,SRR2093568,23314101,ILLUMINA
Pseudocercospora_musae,113226,SRR2093567,23314101,ILLUMINA
```

**Fix:** Download FASTQ directly from EBI, renumber reads, truncate to equal length, normalize coverage, then trim with fastp.

```bash
curl -O -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/008/SRR2093568/SRR2093568.fastq.gz
curl -O -L ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR209/008/SRR2093567/SRR2093567.fastq.gz

seqkit replace -j 8 -p '.+' -r '{nr}' SRR2093567.fastq.gz \
  | ./scripts/fix_fastq_headers --read 1 \
  | pigz -c > $SCRATCH/Pseudocercospora_musae_R1.fq.gz

seqkit replace -j 8 -p '.+' -r '{nr}' SRR2093568.fastq.gz \
  | ./scripts/fix_fastq_headers --read 1 \
  | pigz -c > $SCRATCH/Pseudocercospora_musae_R2.fq.gz

./scripts/enforce_seqpair_readlen \
  in=$SCRATCH/Pseudocercospora_musae_R1.fq.gz \
  in2=$SCRATCH/Pseudocercospora_musae_R2.fq.gz \
  out=$SCRATCH/Pseudocercospora_musae_trunc_R1.fq.gz \
  out2=$SCRATCH/Pseudocercospora_musae_trunc_R2.fq.gz \
  minlen=75

bbnorm.sh \
  in=$SCRATCH/Pseudocercospora_musae_trunc_R1.fq.gz \
  in2=$SCRATCH/Pseudocercospora_musae_trunc_R2.fq.gz \
  out=$SCRATCH/Pseudocercospora_musae_norm_R1.fq.gz \
  out2=$SCRATCH/Pseudocercospora_musae_norm_R2.fq.gz \
  target=30 ecc=t
fastp \
  --in1 $SCRATCH/Pseudocercospora_musae_norm_R1.fq.gz \
  --in2 $SCRATCH/Pseudocercospora_musae_norm_R2.fq.gz \
  --out1 rnaseq_reads/Pseudocercospora_musae_norm_R1.fastq.gz \
  --out2 rnaseq_reads/Pseudocercospora_musae_norm_R2.fastq.gz \
  --thread 16 --detect_adapter_for_pe \
  --cut_front --cut_front_window_size 1 --cut_front_mean_quality 5 \
  --cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 5 \
  --cut_right --cut_right_window_size 4 --cut_right_mean_quality 5 \
  --length_required 25
```

### 3. Nextflow thinks jobs moved to different queue are failed

Change the request of queue out of the nextflow.conf and conf/profile_funannotate.conf for now so that restriction is not in there

### 4. syncing finished but not fully copied results

python sync_predict_results.py --dry-run
python sync_predict_results.py 


