include { tablesDir } from '../../common/utils.nf'

// Per T-014 §D.2: always builds the single, full/unscoped master samples+species
// table, regardless of --taxon. taxonRowFilter() already restricts which genomes
// are *processed* under --taxon (at the base genome channel); MERGE_SAMPLES no
// longer also restricts the table it writes, so a --taxon run no longer produces
// a separate tables/<Taxon>/ subset -- see also the retired `matched` input this
// process used to take.
//
// species.parquet is published atomically with a freshness manifest (see
// docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md); samples.parquet
// stays on the old publishDir copy path -- it isn't part of the view catalog yet.
process MERGE_SAMPLES {
    label      'merge'
    publishDir path: { tablesDir() }, mode: 'copy', pattern: 'samples.parquet'

    input:
    path(samples)

    output:
    path "samples.parquet", emit: samples
    path "species.parquet", emit: species

    script:
    def shrinkFlag = params.merge_all.toBoolean() ? '--enforce-no-shrink' : ''
    """
    python3 ${projectDir}/bin/subset_samples.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${projectDir}/bin/build_species_table.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('samples.csv.gz', sample_size=-1)) TO 'samples.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('species.csv.gz', sample_size=-1)) TO 'species.parquet' (FORMAT PARQUET);"
    rm -f samples.csv.gz species.csv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         species \\
        --local-parquet species.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_SAMPLES \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json \\
        ${shrinkFlag}
    """

    stub:
    """
    python3 ${projectDir}/bin/subset_samples.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        samples.csv.gz
    python3 ${projectDir}/bin/build_species_table.py \\
        --samples ${samples} \\
        --key     LOCUSTAG \\
        -o        species.csv.gz
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT * FROM read_csv_auto('samples.csv.gz', sample_size=-1)) TO 'samples.parquet' (FORMAT PARQUET);"
    duckdb -c "COPY (SELECT * FROM read_csv_auto('species.csv.gz', sample_size=-1)) TO 'species.parquet' (FORMAT PARQUET);"
    rm -f samples.csv.gz species.csv.gz
    python3 ${projectDir}/bin/publish_table.py \\
        --table         species \\
        --local-parquet species.parquet \\
        --tables-dir    ${tablesDir()} \\
        --built-by      MERGE_SAMPLES \\
        --merge-run-id  ${workflow.sessionId} \\
        --schema        ${projectDir}/../sql/table_schema.json
    """
}
