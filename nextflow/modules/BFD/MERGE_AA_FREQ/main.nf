include { tablesDir; coverageCheck } from '../../common/utils.nf'

process MERGE_AA_FREQ {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "aa_freq.parquet", emit: parquet

    script:
    // DuckDB reads the per-genome .csv.gz files listed in the manifest directly
    // (bin/manifest_to_parquet.py). The former zcat-concatenate loop +
    // read_csv_auto(sample_size=-1) was 4.5x slower on 1,020 genomes and failed
    // at ~23k genomes with "Trying to seek a compressed CSV File" (2026-11-05).
    """
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    python3 ${projectDir}/bin/manifest_to_parquet.py --manifest ${manifest} \\
        --out aa_freq.parquet --threads ${task.cpus} \\
        --memory-limit ${(task.memory.toGiga() * 0.75) as int}GB --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('aa_freq.parquet')}
    """

    stub:
    """
    printf 'species_prefix,amino_acid,frequency\\n' | gzip > aa_freq.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('aa_freq.csv.gz', sample_size=-1)) TO 'aa_freq.parquet' (FORMAT PARQUET);"
    """
}
