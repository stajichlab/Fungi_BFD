process RUN_TARGETP {
    tag        "${locustag}"
    label      'targetp'
    storeDir   "${params.outdir}/targetP"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}_summary.targetp2.gz"), emit: summary

    script:
    """
    TMPD=\$(mktemp -d)
    module load targetp
    targetp -batch 50 -tmp \$TMPD -format short \\
        -fasta ${proteins} -org non-pl -prefix ${basename}
    pigz -f ${basename}_summary.targetp2
    rm -rf \$TMPD
    """

    stub:
    """
    printf '# TargetP-2.0\\n' | gzip > ${basename}_summary.targetp2.gz
    """
}
