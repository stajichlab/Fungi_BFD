// Merge all per-species SRA query CSVs into a single named manifest.
// Output: {stem}.rnaseq_sra.csv written alongside the input samples file.
// Columns: species_tag, taxonid, sra_accession, spots, platform
process COLLECT_SRA_QUERY {
    publishDir { file(params.samples).parent.toAbsolutePath().toString() }, mode: 'copy'

    cpus   1
    memory '1 GB'
    time   '10m'

    input:
    path(query_csvs)
    val(stem)

    output:
    path("${stem}.rnaseq_sra.csv"), emit: manifest

    script:
    """
    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${stem}.rnaseq_sra.csv
    for f in ${query_csvs}; do
        tail -n +2 "\$f" >> ${stem}.rnaseq_sra.csv
    done
    NSPECIES=\$(awk -F',' 'NR>1{print \$1}' ${stem}.rnaseq_sra.csv | sort -u | wc -l)
    NACCESSIONS=\$(awk 'NR>1' ${stem}.rnaseq_sra.csv | wc -l)
    echo "[INFO] ${stem}.rnaseq_sra.csv: \$NACCESSIONS accessions across \$NSPECIES species with RNA-seq data"
    """

    stub:
    """
    printf 'species_tag,taxonid,sra_accession,spots,platform,layout\n' > ${stem}.rnaseq_sra.csv
    """
}
