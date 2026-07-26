process RUN_TMHMM {
    tag        "${locustag}"
    label      'tmhmm'
    storeDir   "${params.outdir}/tmhmm"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.tmhmm_short.tsv.gz"),   emit: short_tsv
        path("${basename}.tmhmm_results.tsv.gz"), emit: full_tsv

    script:
    """
    module load tmhmm
    tmhmm --noplot         < ${proteins} > ${basename}.tmhmm_results.tsv
    tmhmm --short --noplot < ${proteins} > ${basename}.tmhmm_short.tsv
    pigz ${basename}.tmhmm_results.tsv ${basename}.tmhmm_short.tsv
    """

    stub:
    """
    printf '# TMHMM\\n' | gzip > ${basename}.tmhmm_results.tsv.gz
    printf '# TMHMM\\n' | gzip > ${basename}.tmhmm_short.tsv.gz
    """
}
