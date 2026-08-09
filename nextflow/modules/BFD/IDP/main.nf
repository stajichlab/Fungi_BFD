include { hashBucketForType } from '../../common/utils.nf'

process RUN_IDP {
    tag        "${meta.locustag}"
    label      'idp'
    storeDir   { "${params.outdir}/aiupred/${hashBucketForType('aiupred', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.aiupred.txt.gz"),     emit: raw
        path("${meta.locustag}.idp.csv.gz"),         emit: idp_csv
        path("${meta.locustag}.idp_summary.csv.gz"), emit: idp_summary_csv

    script:
    """
    module load aiupred
    aiupred.py -i ${proteins} -o ${meta.locustag}.aiupred.txt
    pigz ${meta.locustag}.aiupred.txt
    python3 ${projectDir}/bin/gather_AIUPred.py ${meta.locustag}.aiupred.txt.gz \\
        --outfile      ${meta.locustag}.idp.csv \\
        --outfilesum   ${meta.locustag}.idp_summary.csv
    pigz ${meta.locustag}.idp.csv ${meta.locustag}.idp_summary.csv
    """

    stub:
    """
    printf '' | gzip > ${meta.locustag}.aiupred.txt.gz
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > ${meta.locustag}.idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > ${meta.locustag}.idp_summary.csv.gz
    """
}
