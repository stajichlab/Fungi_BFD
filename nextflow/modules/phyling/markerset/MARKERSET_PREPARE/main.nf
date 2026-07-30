process MARKERSET_PREPARE {
    tag   "${markerset}"
    label 'phyling_markerset'
    storeDir "${params.phylo_outdir}/markersets"

    input:
        val markerset

    output:
        tuple val(markerset), path("${markerset}"), emit: ready

    script:
    """
    src="${params.markerset_db}/${markerset}"
    if [ -d "\$src" ]; then
        ln -sfn "\$src" "${markerset}"
        echo "[INFO] Using pre-extracted markerset: \$src"
    else
        TARBALL=\$(ls "${params.markerset_db}/${markerset}".*.tar.gz 2>/dev/null | head -1)
        if [ -z "\$TARBALL" ]; then
            echo "ERROR: markerset '${markerset}' not found in ${params.markerset_db}" >&2
            exit 1
        fi
        echo "[INFO] Extracting \$TARBALL ..."
        mkdir -p "${markerset}"
        tar -xzf "\$TARBALL" -C "${markerset}" --strip-components=1
    fi
    """

    stub:
    """
    mkdir -p "${markerset}/hmms"
    touch "${markerset}/hmms/dummy.hmm"
    touch "${markerset}/links_to_ODB10.txt"
    """
}
