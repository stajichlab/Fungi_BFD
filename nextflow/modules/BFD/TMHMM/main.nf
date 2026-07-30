process RUN_TMHMM {
    tag        "${meta.locustag}"
    label      'tmhmm'
    storeDir   "${params.outdir}/tmhmm"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.tmhmm_short.tsv.gz"),   emit: short_tsv
        path("${meta.id}.tmhmm_results.tsv.gz"), emit: full_tsv

    script:
    """
    module load tmhmm
    tmhmm --noplot         < ${proteins} > ${meta.id}.tmhmm_results.tsv
    tmhmm --short --noplot < ${proteins} > ${meta.id}.tmhmm_short.tsv
    pigz ${meta.id}.tmhmm_results.tsv ${meta.id}.tmhmm_short.tsv
    """

    stub:
    """
    printf '# TMHMM\\n' | gzip > ${meta.id}.tmhmm_results.tsv.gz
    printf '# TMHMM\\n' | gzip > ${meta.id}.tmhmm_short.tsv.gz
    """
}
