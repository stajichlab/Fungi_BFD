// Count transcripts in a species' shared Trinity-GG fasta so FUNANNOTATE_RNASEQ can
// branch low-yield genome-guided assemblies into TRINITY_STANDALONE. Deliberately NOT
// folded into RNASEQ_PREPARE's own storeDir outputs: adding a new required storeDir
// output there would make every already-cached species (not just the failing ones)
// look incomplete and force a full funannotate-train re-run on the next pipeline pass.
// Has its own separate storeDir instead, so a species already counted is skipped on
// the next pipeline pass without touching RNASEQ_PREPARE's cache group at all.
process COUNT_TRINITY_TRANSCRIPTS {
    tag "$species_tag"

    cpus   1
    memory '1 GB'
    time   '15m'

    storeDir "${launchDir}/rnaseq_data/counts"

    input:
    tuple val(species_tag), path(trinity_fa)

    output:
    tuple val(species_tag), path("${species_tag}.n_transcripts.txt"), emit: counted

    script:
    """
    grep -c '^>' ${trinity_fa} > ${species_tag}.n_transcripts.txt || echo 0 > ${species_tag}.n_transcripts.txt
    """
}
