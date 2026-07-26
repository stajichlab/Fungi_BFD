include { tablesDir } from '../../common/utils.nf'

process MERGE_WOLFPSORT {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(results)

    output:
        path("wolfpsort.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_wolfpsort.py -o wolfpsort.csv ${results}
    pigz wolfpsort.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,localization,score\\n' | gzip > wolfpsort.csv.gz
    """
}
