process ANTISMASH_RUN {
    tag "$out"

    cpus   8
    memory '16 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table), val(antismash_db_dir)

    output:
    tuple val(out), path("${out}/antismash_local/**")

    script:
    def gbk = "${params.target}/${out}/predict_results/${out}.gbk"
    """
    # Accept a compressed prediction (.gbk.gz); antismash needs it uncompressed, so
    # inflate a local copy in the work dir when only the gzipped form is present.
    GBK="${gbk}"
    if [ ! -f "\$GBK" ] && [ -f "${gbk}.gz" ]; then
        zcat "${gbk}.gz" > ${out}.predict.gbk
        GBK=${out}.predict.gbk
    fi
    if [ ! -f "\$GBK" ]; then
        echo "ERROR: predict GBK not found: ${gbk}[.gz]" >&2
        exit 1
    fi
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load apptainer
    mkdir -p ${out}/antismash_local
    # ${antismash_db_dir} is bind-mounted read-only and pointed at via --databases --
    # see SETUP_ANTISMASH_DB/main.nf for the version-match assumption this relies on
    # when antismash_databases is pointed at a pre-existing (not pipeline-downloaded)
    # database directory.
    export TMPDIR=\${SCRATCH:-/tmp}
    SING_BINDS="--bind ${antismash_db_dir}:${antismash_db_dir},\${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    apptainer exec \\
        \${SING_BINDS} \\
        ${params.antismash_sif} \\
        antismash --taxon ${params.antismash_taxon} \\
        --databases ${antismash_db_dir} \\
        --output-dir ${out}/antismash_local \\
        --genefinding-tool none \\
        --fullhmmer --clusterhmmer --cb-general --pfam2go \\
        -c ${task.cpus} \\
        \$GBK
    pigz ${out}/antismash_local/*.json
    """

    stub:
    """
    mkdir -p ${out}/antismash_local
    touch ${out}/antismash_local/${out}.json.gz
    touch ${out}/antismash_local/index.html
    """
}
