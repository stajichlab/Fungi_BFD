include { tablesDir } from '../../common/utils.nf'

process MERGE_GENE_STATS {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "gene_info.parquet",        emit: gene_info
    path "gene_exons.parquet",       emit: gene_exons
    path "gene_CDS.parquet",         emit: gene_CDS
    path "gene_introns.parquet",     emit: gene_introns
    path "gene_transcripts.parquet", emit: gene_transcripts
    path "gene_trnas.parquet",       emit: gene_trnas
    path "gene_proteins.parquet",    emit: gene_proteins

    script:
    """
    set -euo pipefail
    TMPD=\${SCRATCH:-\$PWD}
    for type in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        first=1
        while IFS=\$'\\t' read -r f _mtime _size; do
            case "\$f" in *.\${type}.csv.gz) ;; *) continue ;; esac
            if [ "\$first" = "1" ]; then zcat "\$f"; first=0
            else zcat "\$f" | tail -n +2; fi
        done < ${manifest} > "\$TMPD/\${type}.csv"
        # Plain CSV, not .csv.gz: with sample_size=-1 DuckDB re-reads the file
        # and fails on gzip with "Trying to seek a compressed CSV File"
        # (2026-11-05, ~23k genomes). The plain copy lives in node-local
        # \$SCRATCH when the job has one, and is removed right after.
        module load duckdb 2>/dev/null || true
        duckdb -c "COPY (SELECT * FROM read_csv_auto('\$TMPD/\${type}.csv', sample_size=-1)) TO '\${type}.parquet' (FORMAT PARQUET);"
        rm -f "\$TMPD/\${type}.csv"
    done
    """

    stub:
    """
    for f in gene_info gene_exons gene_CDS gene_introns gene_transcripts gene_trnas gene_proteins; do
        printf 'id\\n' | gzip > \${f}.csv.gz
        module load duckdb 2>/dev/null || true
        duckdb -c "COPY (SELECT * FROM read_csv_auto('\${f}.csv.gz', sample_size=-1)) TO '\${f}.parquet' (FORMAT PARQUET);"
        rm -f \${f}.csv.gz
    done
    """
}
