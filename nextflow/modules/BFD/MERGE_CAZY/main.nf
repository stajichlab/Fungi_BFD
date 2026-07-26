include { tablesDir } from '../../common/utils.nf'

process MERGE_CAZY {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path(overviews)
        path(cazymes)

    output:
        path("cazy.overview.csv.gz"),    emit: overview_csv
        path("cazy.cazymes_hmm.csv.gz"), emit: hmm_csv

    script:
    """
    export PATH="${projectDir}/bin:\$PATH"
    merge_cazy.py \\
        --overviews ${overviews} \\
        --cazymes   ${cazymes} \\
        --out-overview cazy.overview.csv \\
        --out-hmm      cazy.cazymes_hmm.csv
    pigz cazy.overview.csv cazy.cazymes_hmm.csv
    """

    stub:
    """
    printf 'species_prefix,protein_id,EC,cazyme_fam,sub_fam,diamond_fam,substrate,toolcount\\n' | gzip > cazy.overview.csv.gz
    printf 'species_prefix,HMM_id,profile_length,protein_id,protein_length,evalue,q_start,q_end,s_start,s_end,coverage\\n' | gzip > cazy.cazymes_hmm.csv.gz
    """
}
