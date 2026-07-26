process RUN_IDP {
    tag        "${meta.locustag}"
    label      'idp'
    storeDir   "${params.outdir}/aiupred"

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.id}.aiupred.txt.gz"),     emit: raw
        path("${meta.id}.idp.csv.gz"),         emit: idp_csv
        path("${meta.id}.idp_summary.csv.gz"), emit: idp_summary_csv

    script:
    """
    module load aiupred
    aiupred.py -i ${proteins} -o ${meta.id}.aiupred.txt
    pigz ${meta.id}.aiupred.txt
    python3 ${params.scripts}/gather_AIUPred.py ${meta.id}.aiupred.txt.gz \\
        --outfile      ${meta.id}.idp.csv \\
        --outfilesum   ${meta.id}.idp_summary.csv
    pigz ${meta.id}.idp.csv ${meta.id}.idp_summary.csv
    """

    stub:
    """
    printf '' | gzip > ${meta.id}.aiupred.txt.gz
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > ${meta.id}.idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > ${meta.id}.idp_summary.csv.gz
    """
}
