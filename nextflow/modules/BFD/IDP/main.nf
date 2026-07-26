process RUN_IDP {
    tag        "${locustag}"
    label      'idp'
    storeDir   "${params.outdir}/aiupred"

    input:
        tuple val(locustag), val(basename), val(species), val(strain), path(proteins)

    output:
        path("${basename}.aiupred.txt.gz"),     emit: raw
        path("${basename}.idp.csv.gz"),         emit: idp_csv
        path("${basename}.idp_summary.csv.gz"), emit: idp_summary_csv

    script:
    """
    module load aiupred
    aiupred.py -i ${proteins} -o ${basename}.aiupred.txt
    pigz ${basename}.aiupred.txt
    python3 ${params.scripts}/gather_AIUPred.py ${basename}.aiupred.txt.gz \\
        --outfile      ${basename}.idp.csv \\
        --outfilesum   ${basename}.idp_summary.csv
    pigz ${basename}.idp.csv ${basename}.idp_summary.csv
    """

    stub:
    """
    printf '' | gzip > ${basename}.aiupred.txt.gz
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > ${basename}.idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > ${basename}.idp_summary.csv.gz
    """
}
