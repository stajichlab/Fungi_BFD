include { tablesDir; coverageCheck } from '../../common/utils.nf'

process MERGE_IDP {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
        path 'idp/*'
        path 'sum/*'

    output:
        path("idp.parquet"),         emit: idp
        path("idp_summary.parquet"), emit: summary

    script:
    // DuckDB reads the staged .csv.gz files directly (bin/manifest_to_parquet.py);
    // the file lists are written with the printf builtin (no exec, no ARG_MAX
    // limit). On 1,018 genomes: 98 s -> 8.6 s, identical output (2026-11-06).
    def m2p = "python3 ${projectDir}/bin/manifest_to_parquet.py --threads ${task.cpus} " +
              "--memory-limit ${(task.memory.toGiga() * 0.75) as int}GB"
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    printf '%s\\n' idp/*.idp.csv.gz         > idp.list
    printf '%s\\n' sum/*.idp_summary.csv.gz > sum.list
    ${m2p} --manifest idp.list --out idp.parquet         --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('idp.parquet')}
    ${m2p} --manifest sum.list --out idp_summary.parquet --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('idp_summary.parquet')}
    rm -f idp.list sum.list
    """

    stub:
    """
    printf 'species_prefix,protein_id,IDP_start,IDP_end,IDP_length,mean_score\\n' | gzip > idp.csv.gz
    printf 'species_prefix,protein_id,IDP_residues,IDP_fraction,length\\n'        | gzip > idp_summary.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('idp.csv.gz', sample_size=-1)) TO 'idp.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('idp_summary.csv.gz', sample_size=-1)) TO 'idp_summary.parquet' (FORMAT PARQUET);"
    rm -f idp.csv.gz idp_summary.csv.gz
    """
}
