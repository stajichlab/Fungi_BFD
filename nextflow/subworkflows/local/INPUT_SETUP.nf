//
// INPUT_SETUP — symlink genome_annotation/predict_results files into input/ subdirs
//

include { SETUP_SYMLINKS } from '../../modules/common/SETUP_SYMLINKS/main.nf'

workflow INPUT_SETUP {
    take:
    ch   // tuple(locustag, basename, species, strain)

    main:
    rows_file = ch
        .map { locustag, basename, species, strain -> "${basename}\t${locustag}" }
        .collectFile(name: 'setup_rows.tsv', newLine: true)
    SETUP_SYMLINKS(rows_file)
    done_ch = SETUP_SYMLINKS.out.manifest
        .combine(ch)
        .map { _manifest, locustag, basename, species, strain ->
            tuple(locustag, basename, species, strain)
        }

    emit:
    done = done_ch
}
