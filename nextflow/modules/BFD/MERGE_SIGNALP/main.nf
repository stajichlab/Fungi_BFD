include { tablesDir } from '../../common/utils.nf'

process MERGE_SIGNALP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(gff3s)

    output:
        path("signalp.signal_peptide.csv.gz"), emit: csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_signalp.py -o signalp.signal_peptide.csv ${gff3s}
    pigz signalp.signal_peptide.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,peptide_start,peptide_end,probability\\n' | gzip > signalp.signal_peptide.csv.gz
    """
}
