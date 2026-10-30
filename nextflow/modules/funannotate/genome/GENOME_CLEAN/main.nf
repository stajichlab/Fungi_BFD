process GENOME_CLEAN {
    tag "$asmid"

    // AAFTF is invoked inside the script via
    // `apptainer exec ${params.aaftf_sif} ...` (see SING_BINDS/SING below), not
    // as a Nextflow `process.container` -- this process also needs host-side tools
    // (pigz, Rscript/conda) and the FCS-GX db staged under params.fcsgx_stage_dir.

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
    # AAFTF comes from the AAFTF SIF (replaces `module load AAFTF`); taxonkit
    # comes from params.taxonkit_container via \$TSING. fcs_gx_purge takes its DB via an explicit
    # --db <stage dir>/gxdb/all flag (not the image's AAFTF_DB env default), so
    # only the stage dir needs binding for that (\$SCRATCH is bound below) -- no separate AAFTF_DB bind.
    # \$SCRATCH (node-local, NOT under /bigdata) also needs an explicit bind --
    # fcs_gx_purge reads/writes \$SCRATCH/${asmid}.*.fa[sta] directly, and
    # manually-built \$SING commands get none of Nextflow's automatic
    # task-workdir binding (same missing-bind bug class as GENEMARK_RUN/
    # FUNANNOTATE_PREDICT, confirmed 2026-08-24).
    # FCS-GX DB stage dir = params.fcsgx_stage_dir ('\$SCRATCH' = node-local NVMe; it expands on the task node).
    # Each task gets its own mktemp dir, so tasks on one node never share or delete each other's copy.
    FCS_GX_STAGE_ROOT="${params.fcsgx_stage_dir}"
    mkdir -p "\$FCS_GX_STAGE_ROOT"
    FCS_GX_STAGE="\$(mktemp -d "\$FCS_GX_STAGE_ROOT/fcsgx_XXXXXX")" || { echo "[ERROR] cannot create a stage dir under \$FCS_GX_STAGE_ROOT" >&2; exit 1; }
    trap 'rm -rf "\$FCS_GX_STAGE"' EXIT
    export FCS_GX_SHM_DIR="\$FCS_GX_STAGE/gxdb"
    STAGE_BIND=""
    [ "\$FCS_GX_STAGE_ROOT" != "\$SCRATCH" ] && STAGE_BIND=",\$FCS_GX_STAGE_ROOT:\$FCS_GX_STAGE_ROOT"
    SING_BINDS="--bind \${SCRATCH}:\${SCRATCH},${taxondb}:${taxondb},/dev/shm:/dev/shm\${STAGE_BIND}"
    SING="apptainer exec \${SING_BINDS} \${AAFTF_SIF}"
    # taxonkit is not in the AAFTF >= 0.7 image; it runs from params.taxonkit_container.
    TSING="apptainer exec --bind ${taxondb}:${taxondb}:ro ${params.taxonkit_container}"
    # Stage the FCS-GX DB into \$FCS_GX_STAGE (this job's own copy; the EXIT trap above removes it).
    # Persist per-task FCS-GX sync timing outside work/ (cleanup=true).
    export FCS_GX_KEEP_SHM=1
    export FCS_GX_TIMING_LOG="${launchDir}/logs/nextflow/fcs_gx_shm_timing.tsv"
    source ${projectDir}/bin/setup_fcs_shm.sh || { echo "[ERROR] FCS-GX DB staging failed" >&2; exit 1; }
    TAXONKIT_DB=${taxondb}
    phylum=\$(echo ${taxonid} | \$TSING taxonkit --data-dir \$TAXONKIT_DB lineage | \$TSING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result | cut -f3 | \$TSING taxonkit --data-dir \$TAXONKIT_DB name2taxid | cut -f2 | uniq | head -n 1)
    if [ -z "\$phylum" ]; then
    	phylum=\$(echo ${taxonid} | \$TSING taxonkit --data-dir \$TAXONKIT_DB lineage | \$TSING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{K}" --output-ambiguous-result | cut -f3 | \$TSING taxonkit --data-dir \$TAXONKIT_DB name2taxid | uniq | cut -f2 | head -n 1)
	# weird we are getting 2 lines from name2taxid when input is Fungi add the uniq/head -n 1 to ensure only one line
    fi
    if ! [[ "\$phylum" =~ ^[0-9]+\$ ]]; then
        echo "ERROR: could not resolve a numeric taxid for ${asmid} (input taxonid=${taxonid}); got '\$phylum' instead." >&2
        echo "[DEBUG] taxonkit lineage output:" >&2
        echo ${taxonid} | \$TSING taxonkit --data-dir \$TAXONKIT_DB lineage >&2
        echo "[DEBUG] taxonkit reformat -f {p} output:" >&2
        echo ${taxonid} | \$TSING taxonkit --data-dir \$TAXONKIT_DB lineage | \$TSING taxonkit --data-dir \$TAXONKIT_DB reformat -f "{p}" --output-ambiguous-result >&2
        echo "[DEBUG] TAXONKIT_DB=\$TAXONKIT_DB" >&2
        echo "This usually means the phylum/kingdom name is missing from the taxonkit names.dmp in TAXONKIT_DB (outdated or mismatched taxonomy dump). Update TAXONKIT_DB and retry." >&2
        exit 1
    fi
    echo "[INFO] Phylum for ${asmid} (taxonid=${taxonid}): \$phylum"
    echo "[INFO] Decompressing and cleaning genome for ${asmid}..."
    pigz -dc ${genome_gz} > \$SCRATCH/${asmid}.raw.fa
    # AAFTF >= 0.7: fcs_gx_purge -o is an output DIRECTORY. It gets <input prefix>.fcsgx_cleaned.fasta and
    # <input prefix>.fcs_gx_report.txt; the input prefix here is ${asmid}.raw.
    FCS_OUT=\$SCRATCH/${asmid}.fcs_out
    mkdir -p \$FCS_OUT
    \$SING AAFTF fcs_gx_purge -d "\$FCS_GX_SHM_DIR/all" \
        -i \$SCRATCH/${asmid}.raw.fa -c ${task.cpus} \
        -o \$FCS_OUT \
        -t "\$phylum" -w \$SCRATCH/${asmid}.fcs_report
    PURGED=\$FCS_OUT/${asmid}.raw.fcsgx_cleaned.fasta
    # fcs_gx_purge's exit code is not a reliable success signal; require the cleaned file.
    if [ ! -s "\$PURGED" ]; then
        echo "ERROR: fcs_gx_purge produced no \$PURGED for ${asmid}" >&2
        exit 1
    fi
    mkdir -p ${launchDir}/input_clean_genomes/clean
    ${params.clean_script} --len ${params.min_contig_len} < \$PURGED > \$SCRATCH/${asmid}.clean.fa
    if [ ! -s \$SCRATCH/${asmid}.clean.fa ]; then
        echo "ERROR: no contigs left for ${asmid} after the FCS-GX purge and length filter" >&2
        exit 1
    fi
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

    pigz -c \$PURGED > ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta.gz
    [ -f \$FCS_OUT/${asmid}.raw.fcs_gx_report.txt ] && \
        pigz -c \$FCS_OUT/${asmid}.raw.fcs_gx_report.txt > ${launchDir}/input_clean_genomes/clean/${asmid}.fcs_gx_report.txt.gz
    rm -rf \$FCS_OUT
    """

    stub:
    """
    echo ">stub_${asmid}" | pigz -c > ${asmid}.fa.gz
    : > TO_ADD_TO_SUPRESS.csv
    mkdir -p ${launchDir}/input_clean_genomes/clean
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.purge.fasta
    touch ${launchDir}/input_clean_genomes/clean/${asmid}.fcs_gx_report.txt
    """
}
