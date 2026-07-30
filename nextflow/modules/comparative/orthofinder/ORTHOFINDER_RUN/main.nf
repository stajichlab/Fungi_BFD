process ORTHOFINDER_RUN {
    label    'comparative_orthofinder'
    tag      "orthofinder"
    storeDir "${params.outdir}/${params.project}/orthofinder/run"

    input:
    val proteins_dir

    output:
    path "orthofinder_out", emit: out_dir

    script:
    def msa_flag = params.orthofinder_msa.toBoolean() ? '-M msa' : ''
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load orthofinder
    orthofinder -f "${proteins_dir}" \\
        -t ${task.cpus} -a ${task.cpus} \\
        ${msa_flag} \\
        -o orthofinder_out
    """

    stub:
    """
    mkdir -p orthofinder_out/Orthogroups orthofinder_out/Comparative_Genomics_Statistics
    printf 'Orthogroup\n' > orthofinder_out/Orthogroups/Orthogroups.tsv
    printf 'Number of species\t1\n' > orthofinder_out/Comparative_Genomics_Statistics/Statistics_Overall.tsv
    """
}
