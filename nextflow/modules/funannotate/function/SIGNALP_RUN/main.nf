process SIGNALP_RUN {
    tag "$out"

    cpus   8
    memory '16 GB'
    time   '12h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/annotate_misc/signalp.results.txt")

    script:
    def proteins = "${params.target}/${out}/predict_results/${out}.proteins.fa"
    """
    PROTEINS="${proteins}"
    if [ ! -f "\$PROTEINS" ] && [ -f "${proteins}.gz" ]; then
        # predict compresses proteins.fa; inflate into node-local scratch for the tool
        PROTEINS="\${SCRATCH:?SCRATCH not set}/${out}.proteins.fa"
        pigz -dc "${proteins}.gz" > "\$PROTEINS"
    fi
    if [ ! -f "\$PROTEINS" ]; then
        echo "ERROR: protein FASTA not found: ${proteins}[.gz]" >&2
        exit 1
    fi
    module load signalp/6-gpu
    TMPDIR=\${SCRATCH:-/tmp}
    signalp6 -od \$TMPDIR/${out}_signalp \\
        -org euk --mode fast -format txt \\
        -fasta \$PROTEINS \\
        --write_procs ${task.cpus} -bs 16
    mkdir -p ${out}/annotate_misc
    cp \$TMPDIR/${out}_signalp/prediction_results.txt ${out}/annotate_misc/signalp.results.txt
    rm -rf \$TMPDIR/${out}_signalp
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/signalp.results.txt
    """
}
