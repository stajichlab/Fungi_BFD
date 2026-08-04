process GENOME_CLEAN {
    tag "$asmid"

    // container '/rhome/jstajich/projects/AAFTF/AAFTF_v0.6.1-signed.sif'

    // Nextflow skips this task when input_clean_genomes/<asmid>.fa.gz already exists.
    storeDir "${launchDir}/input_clean_genomes"

    cpus   16
    memory '450 GB'
    time   '6h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path(genome_gz), val(taxonid), val(taxondb)

    output:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table),
          path("${asmid}.fa.gz"), val(taxonid), emit: genome

    script:
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    # Ensure /dev/shm/gxdb is present on this node; register for cleanup when done.
    source ${projectDir}/bin/setup_fcs_shm.sh
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    TAXONKIT_DB=${taxondb}
    module load taxonkit
    phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
    if [ -z "\$phylum" ]; then
    	phylum=\$(echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
	# weird we are getting 2 lines from name2taxid when input is Fungi add the uniq/head -n 1 to ensure only one line
    fi
    if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
        echo "ERROR: could not resolve a numeric taxid for ${asmid} (input taxonid=${taxonid}); got '\$phylum' instead." >&2
        echo "[DEBUG] taxonkit lineage output:" >&2
        echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage >&2
        echo "[DEBUG] taxonkit reformat -f {p} output:" >&2
        echo ${taxonid} | taxonkit --data-dir \$TAXONKIT_DB lineage | taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result >&2
        echo "[DEBUG] TAXONKIT_DB=\$TAXONKIT_DB" >&2
        echo "This usually means the phylum/kingdom name is missing from the taxonkit names.dmp in TAXONKIT_DB (outdated or mismatched taxonomy dump). Update TAXONKIT_DB and retry." >&2
        module unload taxonkit
        exit 1
    fi
    module unload taxonkit
    echo "[INFO] Phylum for ${asmid} (taxonid=${taxonid}): \$phylum"
    echo "[INFO] Decompressing and cleaning genome for ${asmid}..."
    module load AAFTF
    pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
        -i \$SCRATCH/${asmid}.raw.fa --cpus ${task.cpus} \
        -o \$SCRATCH/${asmid}.purge.fasta \
        -t "\$phylum" -w \$SCRATCH/${asmid}.fcs_report
    mkdir -p ${launchDir}/input_clean_genomes/clean
    cat \$SCRATCH/${asmid}.purge.fasta | \
        ${params.clean_script} --len ${params.min_contig_len} > \$SCRATCH/${asmid}.clean.fa
    echo "[INFO] Clean genome written: ${asmid}.fa (\$(du -sh \$SCRATCH/${asmid}.clean.fa | cut -f1)); compressing to ${asmid}.fa.gz"
    pigz -c \$SCRATCH/${asmid}.clean.fa > ${asmid}.fa.gz
    rm -f \$SCRATCH/${asmid}.clean.fa
    pigz \$SCRATCH/${asmid}.purge.fasta
    [ -f \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv ] && pigz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv
    mv \$SCRATCH/${asmid}.purge.fasta.gz ${launchDir}/input_clean_genomes/clean/
    [ -f \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ] && \
        mv \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ${launchDir}/input_clean_genomes/clean/
    """

    stub:
    """
    echo ">stub_${asmid}" | pigz -c > ${asmid}.fa.gz
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fcs_gx-taxonomy.tsv
    """
}
