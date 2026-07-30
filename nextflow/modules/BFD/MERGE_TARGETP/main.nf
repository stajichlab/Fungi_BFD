include { tablesDir } from '../../common/utils.nf'

process MERGE_TARGETP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(summaries)

    output:
        path("targetP.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_targetp.py -o targetP.csv ${summaries}
    pigz targetP.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,prediction,probability,cleavage_position_start,cleavage_position_end,cleavage_probability,motif\\n' | gzip > targetP.csv.gz
    """
}
