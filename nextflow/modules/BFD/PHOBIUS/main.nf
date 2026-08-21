include { hashBucketForType } from '../../common/utils.nf'

// Phobius has no public container -- it's an academic-license binary
// (phobius.binf.ku.dk), not redistributable on bioconda/biocontainers, so
// unlike the other new tools in this batch it stays on the host Lmod module
// (phobius/1.01, confirmed present) rather than singularity. Same constraint
// funannotate's own phobius/1.01 dependency (see profile_funannotate.config's
// `Loading requirement: ... phobius/1.01 ...`) already lives with.
process RUN_PHOBIUS {
    tag        "${meta.locustag}"
    label      'phobius'
    storeDir   { "${params.outdir}/phobius/${hashBucketForType('phobius', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.phobius.short.txt.gz"), emit: short_txt

    script:
    """
    module load phobius
    phobius.pl -short ${proteins} > ${meta.locustag}.phobius.short.txt
    pigz ${meta.locustag}.phobius.short.txt
    """

    stub:
    """
    printf 'SEQENCE_ID  TM  SP  Prediction\\n' | gzip > ${meta.locustag}.phobius.short.txt.gz
    """
}
