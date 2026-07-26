process RUN_WOLFPSORT {
    tag        "${meta.locustag}"
    label      'wolfpsort'
    storeDir   "${params.outdir}/wolfpsort"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.wolfpsort.results.txt.gz"), emit: results

    script:
    """
    module load wolfpsort
    cat ${proteins} | runWolfPsortSummary fungi > ${meta.id}.wolfpsort.results.txt
    pigz ${meta.id}.wolfpsort.results.txt
    """

    stub:
    """
    printf '# WoLF PSORT\\n' | gzip > ${meta.id}.wolfpsort.results.txt.gz
    """
}
