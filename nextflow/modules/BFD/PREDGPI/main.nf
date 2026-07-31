include { hashBucketForType } from '../../common/utils.nf'

process RUN_PREDGPI {
    tag        "${meta.locustag}"
    label      'predgpi'
    storeDir   { "${params.outdir}/predgpi/${hashBucketForType('predgpi', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.predgpi.gff3.gz"), emit: gff3

    script:
    """
    module load predgpi
    predgpi.py -f ${proteins} -m gff3 -o ${meta.locustag}.predgpi.gff3
    pigz ${meta.locustag}.predgpi.gff3
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${meta.locustag}.predgpi.gff3.gz
    """
}
