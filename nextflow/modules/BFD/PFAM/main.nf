process RUN_PFAM {
    tag        "${meta.locustag}"
    label      'pfam'
    storeDir   "${params.outdir}/pfam_hmmscan"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.domtblout.gz"),    emit: domtbl
        path("${meta.id}.tblout.gz"),  emit: tblout

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
        --domtbl    ${meta.id}.domtblout \\
        --tblout    ${meta.id}.tblout \\
        \$PFAM_DB/Pfam-A.hmm ${proteins} > /dev/null
    pigz ${meta.id}.domtblout ${meta.id}.tblout
    """

    stub:
    """
    printf '#\\n' | gzip > ${meta.id}.domtblout.gz
    printf '' | gzip     > ${meta.id}.tblout.gz
    """
}
