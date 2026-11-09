// Download the Transporter Classification Database (tcdb.org) reference
// FASTA and build a BLAST db from it. storeDir caches this at
// params.tcdb_dbdir so it downloads/builds at most once across all runs, same
// pattern as SETUP_ANTISMASH_DB in the funannotate profile.
process SETUP_TCDB_DB {
    label 'setup'

    storeDir { params.tcdb_dbdir }

    output:
        path("tcdb.fasta"),                emit: fasta
        // tcdb.p* (not tcdb.*): the pattern must not match tcdb.fasta, or the
        // same file is staged twice into RUN_TCDB. Covers every makeblastdb
        // sidecar (phr/pin/psq, plus pdb/pot/ptf/pto/pjs for a v5 db).
        path("tcdb.p*"),                   emit: blastdb

    script:
    """
    set -euo pipefail
    curl -fsSL https://tcdb.org/public/tcdb -o tcdb.fasta.raw
    # Some TCDB headers have a space before the last pipe, e.g.
    # '>gnl|TC-DB|CHK16680 |1.C.123.1.8 ...'. blastp reports only the text up
    # to the first space as the subject ID, so the TC number would be lost
    # (19 of 24,281 headers on 2026-11-06). Remove that space.
    sed -E 's/^(>gnl\\|TC-DB\\|[^| ]+) +\\|/\\1|/' tcdb.fasta.raw > tcdb.fasta
    rm -f tcdb.fasta.raw
    # Fail loudly if any header still lacks a TC number in its subject ID.
    BAD=\$(grep '>' tcdb.fasta | awk '{print \$1}' | awk -F'|' '\$4 !~ /^[0-9]+\\.[A-Z]\\.[0-9]+\\.[0-9]+\\.[0-9]+/' | wc -l)
    if [ "\$BAD" -ne 0 ]; then echo "ERROR: \$BAD TCDB headers without a TC number" >&2; exit 1; fi
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    SING_BINDS="--bind \${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    apptainer exec \${SING_BINDS} ${params.blastp_sif} makeblastdb -in tcdb.fasta -dbtype prot -blastdb_version 4 -out tcdb
    """

    stub:
    """
    printf '>1.A.1.1.1|test\\nMKV\\n' > tcdb.fasta
    touch tcdb.phr tcdb.pin tcdb.psq
    """
}
