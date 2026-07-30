// Download and extract NCBI taxdump once; storeDir caches it at params.taxondb so
// subsequent runs skip this entirely.
process SETUP_TAXONDB {
    storeDir params.taxondb

    cpus   1
    memory '4 GB'
    time   '1h'

    output:
    path "names.dmp",    emit: ready
    path "nodes.dmp"
    path "merged.dmp"
    path "delnodes.dmp"
    path "division.dmp"
    path "gencode.dmp"
    path "citations.dmp"

    script:
    """
    set -euo pipefail
    wget --no-verbose https://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
    tar zxf taxdump.tar.gz
    rm taxdump.tar.gz
    """

    stub:
    """
    for f in names.dmp nodes.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp citations.dmp; do
        touch \$f
    done
    """
}
