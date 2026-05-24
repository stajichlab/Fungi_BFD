#!/usr/bin/bash -l

mkdir -p data
pushd data
mkdir -p dna pep cds gff3 trna
popd
SAMPLES=samples.csv
INDIR=$1
if [ -z "$INDIR" ]; then 
	INDIR=genome_annotation
fi
INDIR=$(realpath $INDIR)
IFS=,
tail -n +2 $SAMPLES | while read ASMID SPECIES_IN STRAIN BIOPROJECT NCBI_TAXONID BUSCO_LINEAGE PHYLUM SUBPHYLUM CLASS SUBCLASS ORDER FAMILY GENUS SPECIES TRANSL_TABLE LOCUSTAG
do
	TAG=$(echo -n "$SPECIES $STRAIN" | perl -p -e 's/ /_/g;')
	if [ ! -d "$INDIR/$TAG" ]; then
		echo "cannot find input files $SPECIES $STRAIN $ASMID"
	fi
	ln -sf $INDIR/$TAG/predict_results/${TAG}.gff3 data/gff3/
	ln -sf $INDIR/$TAG/predict_results/${TAG}.proteins.fa data/pep/
	ln -sf $INDIR/$TAG/predict_results/${TAG}.cds-transcripts.fa data/cds/
	ln -sf $INDIR/$TAG/predict_results/${TAG}.scaffolds.fa data/dna/
	ln -sf $INDIR/$TAG/predict_misc/trnascan.no-overlaps.gff3 data/trna/${TAG}.tRNA.gff3
done
