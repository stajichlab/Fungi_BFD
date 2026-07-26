process RUN_WOLFPSORT {
    tag        "${locustag}"
    label      'wolfpsort'
    storeDir   "${params.outdir}/wolfpsort"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.wolfpsort.results.txt.gz"), emit: results

    script:
    """
    module load wolfpsort
    cat ${proteins} | runWolfPsortSummary fungi > ${basename}.wolfpsort.results.txt
    pigz ${basename}.wolfpsort.results.txt
    """

    stub:
    """
    printf '# WoLF PSORT\\n' | gzip > ${basename}.wolfpsort.results.txt.gz
    """
}
