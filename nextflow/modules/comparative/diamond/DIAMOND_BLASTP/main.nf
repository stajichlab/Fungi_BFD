process DIAMOND_BLASTP {
    label    'comparative_diamond'
    tag      "blastp"
    storeDir "${params.outdir}/${params.project}/mcl/search"

    input:
    path combined_faa
    path db

    output:
    path "blastp.tsv", emit: blastp_tsv

    script:
    def sens_flag = params.diamond_more_sensitive.toBoolean() ? '--more-sensitive' : '--sensitive'
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load diamond
    diamond blastp \\
        -d combined \\
        -q ${combined_faa} \\
        ${sens_flag} \\
        --max-target-seqs 500 \\
        -e ${params.mcl_evalue} \\
        --outfmt 6 \\
        --threads ${task.cpus} \\
        -o blastp.tsv
    """

    stub:
    """
    printf 'A\tB\t100.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\n' > blastp.tsv
    """
}
