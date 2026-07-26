include { tablesDir } from '../../common/utils.nf'

process MERGE_AA_FREQ {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "aa_freq.csv.gz", emit: csv

    script:
    """
    first=1
    while IFS=\$'\\t' read -r f _mtime _size; do
        [ -n "\$f" ] || continue
        if [ "\$first" = "1" ]; then zcat "\$f"; first=0
        else zcat "\$f" | tail -n +2; fi
    done < ${manifest} | gzip > aa_freq.csv.gz
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' | gzip > aa_freq.csv.gz
    """
}
