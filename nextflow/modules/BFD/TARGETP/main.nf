include { hashBucketForType } from '../../common/utils.nf'

process RUN_TARGETP {
    tag        "${meta.locustag}"
    label      'targetp'
    storeDir   { "${params.outdir}/targetP/${hashBucketForType('targetP', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}_summary.targetp2.gz"), emit: summary

    script:
    """
    TMPD=\$(mktemp -d)
    module load targetp
    targetp -batch 50 -tmp \$TMPD -format short \\
        -fasta ${proteins} -org non-pl -prefix ${meta.locustag}
    pigz -f ${meta.locustag}_summary.targetp2
    rm -rf \$TMPD
    """

    stub:
    """
    printf '# TargetP-2.0\\n' | gzip > ${meta.locustag}_summary.targetp2.gz
    """
}
