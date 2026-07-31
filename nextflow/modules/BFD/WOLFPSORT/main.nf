include { hashBucketForType } from '../../common/utils.nf'

process RUN_WOLFPSORT {
    tag        "${meta.locustag}"
    label      'wolfpsort'
    storeDir   { "${params.outdir}/wolfpsort/${hashBucketForType('wolfpsort', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.wolfpsort.results.txt.gz"), emit: results

    script:
    """
    module load wolfpsort
    cat ${proteins} | runWolfPsortSummary fungi > ${meta.locustag}.wolfpsort.results.txt
    pigz ${meta.locustag}.wolfpsort.results.txt
    """

    stub:
    """
    printf '# WoLF PSORT\\n' | gzip > ${meta.locustag}.wolfpsort.results.txt.gz
    """
}
