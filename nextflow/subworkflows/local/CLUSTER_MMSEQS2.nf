include { MMSEQS_CREATEDB }   from '../../modules/comparative/mmseqs2/MMSEQS_CREATEDB/main.nf'
include { MMSEQS_SEARCH }     from '../../modules/comparative/mmseqs2/MMSEQS_SEARCH/main.nf'
include { MMSEQS_CLUST }      from '../../modules/comparative/mmseqs2/MMSEQS_CLUST/main.nf'
include { MMSEQS_CREATETSV }  from '../../modules/comparative/mmseqs2/MMSEQS_CREATETSV/main.nf'

workflow CLUSTER_MMSEQS2 {
    take:
    proteins_dir  // val: path to directory containing {LOCUSTAG}.faa files

    main:
    MMSEQS_CREATEDB(proteins_dir)
    MMSEQS_SEARCH(MMSEQS_CREATEDB.out.db_files)
    MMSEQS_CLUST(MMSEQS_CREATEDB.out.db_files, MMSEQS_SEARCH.out.result_files)
    MMSEQS_CREATETSV(MMSEQS_CREATEDB.out.db_files, MMSEQS_CLUST.out.cluster_files)

    emit:
    tsv = MMSEQS_CREATETSV.out.tsv
}
