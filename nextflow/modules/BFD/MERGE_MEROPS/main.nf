include { tablesDir } from '../../common/utils.nf'

process MERGE_MEROPS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(blasttabs)

    output:
        path("merops.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_merops.py -o merops.csv ${blasttabs}
    pigz merops.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,merops_id,percent_identity,aln_length,mismatches,gap_openings,q_start,q_end,s_start,s_end,evalue,bitscore\\n' | gzip > merops.csv.gz
    """
}
