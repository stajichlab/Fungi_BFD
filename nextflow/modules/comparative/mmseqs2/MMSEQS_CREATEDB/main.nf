process MMSEQS_CREATEDB {
    label    'comparative_mmseqs'
    tag      "createdb"
    storeDir "${params.outdir}/${params.project}/mmseqs2/db"

    input:
    val proteins_dir

    output:
    path "combined.faa", emit: combined_faa
    path "combined_db*", emit: db_files

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load mmseqs2
    cat "${proteins_dir}"/*.faa > combined.faa
    mmseqs createdb combined.faa combined_db
    """

    stub:
    """
    touch combined.faa combined_db combined_db.index combined_db.dbtype
    touch combined_db_h combined_db_h.index combined_db_h.dbtype
    """
}
