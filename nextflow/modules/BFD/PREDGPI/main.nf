process RUN_PREDGPI {
    tag        "${meta.locustag}"
    label      'predgpi'
    storeDir   "${params.outdir}/predgpi"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.predgpi.gff3.gz"), emit: gff3

    script:
    """
    module load predgpi
    predgpi.py -f ${proteins} -m gff3 -o ${meta.id}.predgpi.gff3
    pigz ${meta.id}.predgpi.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${meta.id}.predgpi.gff3.gz
    """
}
