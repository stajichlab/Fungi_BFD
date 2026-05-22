#!/usr/bin/bash -l

mkdir -p dna pep cds gff3
INDIR=$1
if [ -z "$INDIR" ]; then 
  INDIR=/bigdata/stajichlab/shared/projects/1KFG/common_annotate
fi

