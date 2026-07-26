process RUN_PFAM {
    tag        "${locustag}"
    label      'pfam'
    storeDir   "${params.outdir}/pfam_hmmscan"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.domtblout.gz"),    emit: domtbl
        path("${basename}.tblout.gz"),  emit: tblout

    script:
    def mpi_launch = params.pfam_tasks > 1 ? "srun -N ${params.pfam_nodes} -n ${params.pfam_tasks}" : ""
    def mpi_flag   = params.pfam_tasks > 1 ? "--mpi" : ""
    """
    if [ ! -z "${mpi_flag}" ]; then
        module load hmmer/3.4-mpi
    else
        module load hmmer/3.4
    fi
    module load db-pfam
    ${mpi_launch} hmmsearch ${mpi_flag} --cut_ga --noali --cpu ${task.cpus} \\
        --domtbl    ${basename}.domtblout \\
        --tblout    ${basename}.tblout \\
        \$PFAM_DB/Pfam-A.hmm ${proteins} > /dev/null
    pigz ${basename}.domtblout ${basename}.tblout
    """

    stub:
    """
    printf '#\\n' | gzip > ${basename}.domtblout.gz
    printf '' | gzip     > ${basename}.tblout.gz
    """
}
