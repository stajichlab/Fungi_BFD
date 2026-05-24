#!/usr/bin/bash -l

mkdir -p dna pep cds gff3
SAMPLES=samples.csv
INDIR=$1
if [ -z "$INDIR" ]; then 
  INDIR=/bigdata/stajichlab/shared/projects/1KFG/common_annotate
fi
IFS=,
tail -n +2 $SAMPLES | while read ASMID SPECIES_IN STRAIN BIOPROJECT NCBI_TAXONID BUSCO_LINEAGE PHYLUM SUBPHYLUM CLASS SUBCLASS ORDER FAMILY GENUS SPECIES TRANSL_TABLE LOCUSTAG
do
	TAG=$(echo -n "$SPECIES $STRAIN" | perl -p -e 's/ /_/g;')
	if [ ! -d "$INDIR/annotate/$TAG" ]; then
		echo "cannot find input files $SPECIES $STRAIN $ASMID"
	fi
	ln -sf $INDIR/annotate/$TAG/predict_results/${TAG}.gff3 gff3/
	ln -sf $INDIR/annotate/$TAG/predict_results/${TAG}.proteins.fa pep/
	ln -sf $INDIR/annotate/$TAG/predict_results/${TAG}.cds-transcripts.fa cds/
	ln -sf $INDIR/annotate/$TAG/predict_results/${TAG}.scaffolds.fa dna/
done
