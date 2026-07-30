process RUN_SIGNALP {
    tag        "${meta.locustag}"
    label      'signalp'
    storeDir   "${params.outdir}/signalp"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.signalp.gff3.gz"),         emit: gff3
        path("${meta.id}.signalp.results.txt.gz"),  emit: results

    script:
    """
    module load signalp/6-gpu
    OUTD=\$(mktemp -d)
    signalp6 -od \$OUTD -org euk --mode fast -format txt \\
        -fasta ${proteins} --write_procs ${task.cpus} -bs 100
    pigz -c \$OUTD/output.gff3             > ${meta.id}.signalp.gff3.gz
    pigz -c \$OUTD/prediction_results.txt  > ${meta.id}.signalp.results.txt.gz
    rm -rf \$OUTD
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${meta.id}.signalp.gff3.gz
    printf '# SignalP-6.0\\n'   | gzip > ${meta.id}.signalp.results.txt.gz
    """
}
