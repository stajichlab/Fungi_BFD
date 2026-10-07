// IPRSCAN5
process INTERPROSCAN_RUN {
    tag "$out"

    cpus   8
    memory '32 GB'
    time   '60h'

    publishDir "${params.target}", mode: 'copy', overwrite: true

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}/annotate_misc/iprscan.xml")

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
    mkdir -p ${out}/annotate_misc
    module load interproscan
    interproscan.sh -i \$PROTEINS -f XML -o ${out}/annotate_misc/iprscan.xml \\
        -dp -goterms -pa -t p -cpu ${task.cpus}
    """

    stub:
    """
    mkdir -p ${out}/annotate_misc
    touch ${out}/annotate_misc/iprscan.xml
    """
}
