include { tablesDir } from '../../common/utils.nf'

process MERGE_PREDGPI {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("predgpi.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_predgpi.py -o predgpi.csv ${gff3s}
    pigz predgpi.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,source,feature,start,end,score,strand,phase,attributes\\n' | gzip > predgpi.csv.gz
    """
}
