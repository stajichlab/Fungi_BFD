include { hashBucketForType } from '../../common/utils.nf'

process RUN_PFAM {
    tag        "${meta.locustag}"
    label      'pfam'
    storeDir   { "${params.outdir}/pfam_hmmscan/${hashBucketForType('pfam_hmmscan', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.domtblout.gz"),    emit: domtbl
        path("${meta.locustag}.tblout.gz"),  emit: tblout

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
        --domtbl    ${meta.locustag}.domtblout \\
        --tblout    ${meta.locustag}.tblout \\
        \$PFAM_DB/Pfam-A.hmm ${proteins} > /dev/null
    pigz ${meta.locustag}.domtblout ${meta.locustag}.tblout
    """

    stub:
    """
    printf '#\\n' | gzip > ${meta.locustag}.domtblout.gz
    printf '' | gzip     > ${meta.locustag}.tblout.gz
    """
}
