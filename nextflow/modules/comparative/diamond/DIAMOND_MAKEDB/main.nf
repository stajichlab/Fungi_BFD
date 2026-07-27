process DIAMOND_MAKEDB {
    label    'comparative_diamond'
    tag      "makedb"
    storeDir "${params.outdir}/${params.project}/mcl/db"

    input:
    val proteins_dir

    output:
    path "combined.faa",  emit: combined_faa
    path "combined.dmnd", emit: db

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load diamond
    cat "${proteins_dir}"/*.faa > combined.faa
    diamond makedb --in combined.faa -d combined --threads ${task.cpus}
    """

    stub:
    """
    touch combined.faa combined.dmnd
    """
}
