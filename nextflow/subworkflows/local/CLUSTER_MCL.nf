include { DIAMOND_MAKEDB }  from '../../modules/comparative/diamond/DIAMOND_MAKEDB/main.nf'
include { DIAMOND_BLASTP }  from '../../modules/comparative/diamond/DIAMOND_BLASTP/main.nf'
include { MCL_PREPARE }     from '../../modules/comparative/mcl/MCL_PREPARE/main.nf'
include { MCL_RUN }         from '../../modules/comparative/mcl/MCL_RUN/main.nf'

workflow CLUSTER_MCL {
    take:
    proteins_dir  // val: path to directory containing {LOCUSTAG}.faa files

    main:
    DIAMOND_MAKEDB(proteins_dir)
    DIAMOND_BLASTP(DIAMOND_MAKEDB.out.combined_faa, DIAMOND_MAKEDB.out.db)
    MCL_PREPARE(DIAMOND_BLASTP.out.blastp_tsv)
    MCL_RUN(MCL_PREPARE.out.abc)

    emit:
    tsv = MCL_RUN.out.tsv
}
