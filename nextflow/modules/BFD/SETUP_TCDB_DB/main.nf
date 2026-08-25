// Download the Transporter Classification Database (tcdb.org) reference
// FASTA and build a BLAST db from it. storeDir caches this at
// params.tcdb_dbdir so it downloads/builds at most once across all runs, same
// pattern as SETUP_ANTISMASH_DB in the funannotate profile.
process SETUP_TCDB_DB {
    label 'setup'

    storeDir { params.tcdb_dbdir }

    output:
        path("tcdb.fasta"),                emit: fasta
        path("tcdb.{phr,pin,psq}"),        emit: blastdb

    script:
    """
    set -euo pipefail
    curl -fsSL https://tcdb.org/public/tcdb -o tcdb.fasta
    module load apptainer
    export TMPDIR=\${SCRATCH:-/tmp}
    SING_BINDS="--bind \${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} ${params.blastp_sif} makeblastdb -in tcdb.fasta -dbtype prot -out tcdb
    """

    stub:
    """
    printf '>1.A.1.1.1|test\\nMKV\\n' > tcdb.fasta
    touch tcdb.phr tcdb.pin tcdb.psq
    """
}
