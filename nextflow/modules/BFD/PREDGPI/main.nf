process RUN_PREDGPI {
    tag        "${locustag}"
    label      'predgpi'
    storeDir   "${params.outdir}/predgpi"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.predgpi.gff3.gz"), emit: gff3

    script:
    """
    module load predgpi
    predgpi.py -f ${proteins} -m gff3 -o ${basename}.predgpi.gff3
    pigz ${basename}.predgpi.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${basename}.predgpi.gff3.gz
    """
}
