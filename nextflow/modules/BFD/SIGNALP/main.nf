process RUN_SIGNALP {
    tag        "${locustag}"
    label      'signalp'
    storeDir   "${params.outdir}/signalp"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.signalp.gff3.gz"),         emit: gff3
        path("${basename}.signalp.results.txt.gz"),  emit: results

    script:
    """
    module load signalp/6-gpu
    OUTD=\$(mktemp -d)
    signalp6 -od \$OUTD -org euk --mode fast -format txt \\
        -fasta ${proteins} --write_procs ${task.cpus} -bs 100
    pigz -c \$OUTD/output.gff3             > ${basename}.signalp.gff3.gz
    pigz -c \$OUTD/prediction_results.txt  > ${basename}.signalp.results.txt.gz
    rm -rf \$OUTD
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${basename}.signalp.gff3.gz
    printf '# SignalP-6.0\\n'   | gzip > ${basename}.signalp.results.txt.gz
    """
}
