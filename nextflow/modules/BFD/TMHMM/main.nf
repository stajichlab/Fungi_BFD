include { hashBucketForType; plainName; gzInflate } from '../../common/utils.nf'

process RUN_TMHMM {
    tag        "${meta.locustag}"
    label      'tmhmm'
    storeDir   { "${params.outdir}/tmhmm/${hashBucketForType('tmhmm', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.tmhmm_short.tsv.gz"),   emit: short_tsv
        path("${meta.locustag}.tmhmm_results.tsv.gz"), emit: full_tsv

    script:
    def proteins_plain = plainName(proteins)
    """
    ${gzInflate(proteins)}
    module load tmhmm
    tmhmm --noplot         < ${proteins_plain} > ${meta.locustag}.tmhmm_results.tsv
    tmhmm --short --noplot < ${proteins_plain} > ${meta.locustag}.tmhmm_short.tsv
    pigz ${meta.locustag}.tmhmm_results.tsv ${meta.locustag}.tmhmm_short.tsv
    """

    stub:
    """
    printf '# TMHMM\\n' | gzip > ${meta.locustag}.tmhmm_results.tsv.gz
    printf '# TMHMM\\n' | gzip > ${meta.locustag}.tmhmm_short.tsv.gz
    """
}
