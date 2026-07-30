include { tablesDir } from '../../common/utils.nf'

process MERGE_IDP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path 'idp/*'
        path 'sum/*'

    output:
        path("idp.csv.gz"),         emit: idp
        path("idp_summary.csv.gz"), emit: summary

    script:
    """
    first=1
    for f in idp/*.idp.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > idp.csv.gz

    first=1
    for f in sum/*.idp_summary.csv.gz; do
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done | gzip > idp_summary.csv.gz
    """

    stub:
    """
    printf 'protein_id,idp_status,disordered_residues,total_residues\\n' | gzip > idp.csv.gz
    printf 'protein_id,idp_status\\n'                                     | gzip > idp_summary.csv.gz
    """
}
