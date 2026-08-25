process GENOME_CLEAN {
    tag "$asmid"

    // AAFTF (+ in-image taxonkit) is invoked inside the script via
    // `apptainer exec ${params.aaftf_sif} ...` (see SING_BINDS/SING below), not
    // as a Nextflow `process.container` -- this process also needs host-side tools
    // (pigz, Rscript/conda) and the host-staged FCS-GX /dev/shm db.

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
    path "TO_ADD_TO_SUPRESS.csv", emit: suppress

    script:
    """
    if [ ! -f "${genome_gz}" ]; then
        echo "ERROR: genome_gz not found at path: ${genome_gz}" >&2
        exit 1
    fi

    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load miniconda3
    eval "\$(conda shell.bash hook)"
    module load apptainer
    AAFTF_SIF=${params.aaftf_sif}
    SCRATCH=\$(printf '%s' "\${SCRATCH}" | tr -d '\\n\\r')
    # AAFTF + taxonkit now come from the AAFTF SIF (replaces `module load AAFTF`
    # and `module load taxonkit`). The image hardcodes AAFTF_DB=/opt/aaaftf_db and
    # AAFTF fcs_gx_purge needs the host-staged FCS-GX db (/dev/shm/gxdb) visible,
    # so bind the host DB + /dev/shm for every in-container call. \$SCRATCH (node-
    # local, NOT under /bigdata) also needs an explicit bind -- fcs_gx_purge reads/
    # writes \$SCRATCH/${asmid}.*.fa[sta] directly, and manually-built \$SING
    # commands get none of Nextflow's automatic task-workdir binding (same missing-
    # bind bug class as GENEMARK_RUN/FUNANNOTATE_PREDICT, confirmed 2026-08-24).
    SING_BINDS="--bind \${SCRATCH}:\${SCRATCH},${taxondb}:${taxondb},/srv/projects/db/AAFTF_DB:/opt/aaaftf_db,/dev/shm:/dev/shm"
    SING="apptainer exec \${SING_BINDS} \${AAFTF_SIF}"
    # Ensure /dev/shm/gxdb is present on this node; register for cleanup when done.
    # Persist per-task FCS-GX /dev/shm sync timing outside work/ (cleanup=true).
    export FCS_GX_TIMING_LOG="${launchDir}/logs/nextflow/fcs_gx_shm_timing.tsv"
    source ${projectDir}/bin/setup_fcs_shm.sh
    TAXONKIT_DB=${taxondb}
    phylum=\$(echo ${taxonid} | \$SING taxonkit --data-dir \$TAXONKIT_DB lineage | \$SING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | \$SING taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
    if [ -z "\$phylum" ]; then
    	phylum=\$(echo ${taxonid} | \$SING taxonkit --data-dir \$TAXONKIT_DB lineage | \$SING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | \$SING taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
	# weird we are getting 2 lines from name2taxid when input is Fungi add the uniq/head -n 1 to ensure only one line
    fi
    if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
        echo "ERROR: could not resolve a numeric taxid for ${asmid} (input taxonid=${taxonid}); got '\$phylum' instead." >&2
        echo "[DEBUG] taxonkit lineage output:" >&2
        echo ${taxonid} | \$SING taxonkit --data-dir \$TAXONKIT_DB lineage >&2
        echo "[DEBUG] taxonkit reformat -f {p} output:" >&2
        echo ${taxonid} | \$SING taxonkit --data-dir \$TAXONKIT_DB lineage | \$SING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result >&2
        echo "[DEBUG] TAXONKIT_DB=\$TAXONKIT_DB" >&2
        echo "This usually means the phylum/kingdom name is missing from the taxonkit names.dmp in TAXONKIT_DB (outdated or mismatched taxonomy dump). Update TAXONKIT_DB and retry." >&2
        exit 1
    fi
    echo "[INFO] Phylum for ${asmid} (taxonid=${taxonid}): \$phylum"
    echo "[INFO] Decompressing and cleaning genome for ${asmid}..."
    pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    \$SING AAFTF fcs_gx_purge --db /dev/shm/gxdb/all \
        -i \$SCRATCH/${asmid}.raw.fa --cpus ${task.cpus} \
        -o \$SCRATCH/${asmid}.purge.fasta \
        -t "\$phylum" -w \$SCRATCH/${asmid}.fcs_report
    mkdir -p ${launchDir}/input_clean_genomes/clean
    cat \$SCRATCH/${asmid}.purge.fasta | \
        ${params.clean_script} --len ${params.min_contig_len} > \$SCRATCH/${asmid}.clean.fa
    echo "[INFO] Clean genome written: ${asmid}.fa (\$(du -sh \$SCRATCH/${asmid}.clean.fa | cut -f1)); compressing to ${asmid}.fa.gz"
    pigz -c \$SCRATCH/${asmid}.clean.fa > ${asmid}.fa.gz
    rm -f \$SCRATCH/${asmid}.clean.fa

    # Assemblies below params.min_assembly_len bp are almost certainly junk/
    # contamination-only or a failed FCS-GX purge, not a usable genome; flag for suppression.
    : > TO_ADD_TO_SUPRESS.csv
    total_len=\$(pigz -dc ${asmid}.fa.gz | grep -v '^>' | tr -d '\\n\\r' | wc -c)
    if [ "\$total_len" -lt "${params.min_assembly_len}" ]; then
        echo "[WARN] ${asmid} cleaned assembly is only \$total_len bp (< ${params.min_assembly_len}); flagging for suppression" >&2
        printf '%s,Assembly too short %s bp in length\\n' "${asmid}" "\$total_len" >> TO_ADD_TO_SUPRESS.csv
    fi

    pigz \$SCRATCH/${asmid}.purge.fasta
    [ -f \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv ] && pigz \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv
    mv \$SCRATCH/${asmid}.purge.fasta.gz ${launchDir}/input_clean_genomes/clean/
    [ -f \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ] && \
        mv \$SCRATCH/${asmid}.purge.fcs_gx-taxonomy.tsv.gz ${launchDir}/input_clean_genomes/clean/
    """

    stub:
    """
    echo ">stub_${asmid}" | pigz -c > ${asmid}.fa.gz
    : > TO_ADD_TO_SUPRESS.csv
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fcs_gx-taxonomy.tsv
    """
}
