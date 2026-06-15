find genome_annotation -maxdepth 2 -type d -name training | while read src;  
do   
	out=$(basename $(dirname "$src"));
	dest="genome_annotation_training/${out}";
	mkdir -p "genome_annotation_training/${out}";
	rsync -a "$src" "$dest";
	echo "Moved: $out";
done
