// Write zero-byte paired FASTQ placeholder files for species with no SRA data.
// Called only for species whose SRA_QUERY CSV has no data rows, avoiding a
// SLURM job allocation for what would be an immediate empty-file write.
process WRITE_EMPTY_READS {
    tag "$species_tag"

    storeDir "${launchDir}/rnaseq_reads"

    cpus   1
    memory '1 GB'
    time   '5m'

    input:
    val(species_tag)

    output:
    tuple val(species_tag), path("${species_tag}_norm_R1.fastq.gz"), path("${species_tag}_norm_R2.fastq.gz"),
          path("${species_tag}_norm_SE.fastq.gz"), emit: reads

    script:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    echo "[INFO] No SRA data for ${species_tag}; created empty read placeholders"
    """

    stub:
    """
    : > ${species_tag}_norm_R1.fastq.gz
    : > ${species_tag}_norm_R2.fastq.gz
    : > ${species_tag}_norm_SE.fastq.gz
    """
}
