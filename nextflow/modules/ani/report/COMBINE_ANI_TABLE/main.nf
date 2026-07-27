// Merge every group's pair table + names lookup into one labeled, queryable
// table (CSV + DuckDB) so pairwise identities can be looked up across the
// whole run without re-parsing per-group files.
process COMBINE_ANI_TABLE {
    label 'report'

    publishDir "${params.outdir}/${params.ani_method}/${params.compare}", mode: 'copy'

    input:
        path ani_manifest
        path names_manifest

    output:
        path("all_pairs.csv")
        path("ani.duckdb")

    script:
    """
    # duckdb is not in the python:3.12 base image — install it first.
    # --target installs to a persistent bind mount (work/ANI/python_packages -> /tmp/python_packages)
    # so pip doesn't try to write to the read-only /rhome home dir inside the container.
    pip install --quiet --target /tmp/python_packages duckdb
    PYTHONPATH="/tmp/python_packages:\${PYTHONPATH}" \
    python3 ${projectDir}/bin/combine_ani_table.py \
        --ani-manifest   ${ani_manifest} \
        --names-manifest ${names_manifest} \
        --compare-level "${params.compare}" \
        --csv-output all_pairs.csv \
        --db-output  ani.duckdb
    """

    stub:
    """
    touch all_pairs.csv ani.duckdb
    """
}
