//
// INPUT_SETUP — symlink genome_annotation/predict_results files into input/ subdirs
//

include { SETUP_SYMLINKS } from '../../modules/common/SETUP_SYMLINKS/main.nf'

workflow INPUT_SETUP {
    take:
    ch   // meta maps

    main:
    rows_file = ch
        .map { meta -> "${meta.id}\t${meta.locustag}" }
        .collectFile(name: 'setup_rows.tsv', newLine: true)
    SETUP_SYMLINKS(rows_file)
    done_ch = SETUP_SYMLINKS.out.manifest
        .combine(ch)
        .map { _manifest, meta -> meta }

    emit:
    done = done_ch
}
