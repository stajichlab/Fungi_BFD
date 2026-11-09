include { tablesDir; coverageCheck; codeMd5 } from '../../common/utils.nf'

process MERGE_INTERGENIC {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy'

    input:
    path manifest

    output:
    path "gene_intergenic_distances.parquet", emit: parquet

    script:
    // DuckDB reads the per-genome .csv.gz files listed in the manifest directly
    // (bin/manifest_to_parquet.py). The former zcat-concatenate loop +
    // read_csv_auto(sample_size=-1) was 4.5x slower on 1,020 genomes and failed
    // at ~23k genomes with "Trying to seek a compressed CSV File" (2026-11-05).
    """
    # code md5 (task cache key): ${codeMd5(['bin/manifest_to_parquet.py'])}
    set -euo pipefail
    module load duckdb 2>/dev/null || true
    python3 ${projectDir}/bin/manifest_to_parquet.py --manifest ${manifest} \\
        --out gene_intergenic_distances.parquet --threads ${task.cpus} \\
        --memory-limit ${(task.memory.toGiga() * 0.75) as int}GB --temp-dir \${SCRATCH:-\$PWD}
    ${coverageCheck('gene_intergenic_distances.parquet')}
    """

    stub:
    """
    printf 'species_prefix,left_gene,right_gene,distance\\n' | gzip > gene_intergenic_distances.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('gene_intergenic_distances.csv.gz', sample_size=-1)) TO 'gene_intergenic_distances.parquet' (FORMAT PARQUET);"
    """
}
