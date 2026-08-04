// Download, checksum, and extract NCBI taxdump once; storeDir caches it at
// params.taxondb (a local persistent directory, e.g. lib/taxdump/) so
// subsequent runs skip this entirely. This is the classic taxdump (not
// new_taxdump) -- taxonkit's lineage/reformat/name2taxid only need
// names.dmp/nodes.dmp/merged.dmp/delnodes.dmp, all present here; new_taxdump
// only adds extras (host.dmp, images.dmp, etc.) taxonkit doesn't use.
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
    wget --no-verbose https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
    wget --no-verbose https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz.md5
    md5sum -c taxdump.tar.gz.md5
    tar zxf taxdump.tar.gz
    rm taxdump.tar.gz taxdump.tar.gz.md5
    """

    stub:
    """
    for f in names.dmp nodes.dmp merged.dmp delnodes.dmp division.dmp gencode.dmp citations.dmp; do
        touch \$f
    done
    """
}
