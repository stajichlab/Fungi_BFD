include { tablesDir } from '../../common/utils.nf'

process MERGE_TMHMM {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(tsvs)

    output:
        path("tmhmm.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_tmhmm.py -o tmhmm.csv ${tsvs}
    pigz tmhmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,len,ExpAA,First60,PredHel,Topology\\n' | gzip > tmhmm.csv.gz
    """
}
