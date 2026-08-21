include { hashBucketForType } from '../../common/utils.nf'

// ncRNA (non-coding RNA) gene calling against Rfam covariance models --
// complements RUN_TRNASCAN (tRNAs specifically) with rRNAs, snoRNAs,
// riboswitches, and other structured ncRNA families. Reuses the already
// cmpress'd shared Rfam 14.5 DB at /srv/projects/db/rfam/14.5 (Rfam.cm +
// Rfam.clanin both confirmed present 2026-08-20) -- no download needed.
// Follows the standard genome-wide Rfam annotation recipe: cmscan --cut_ga
// (trusted per-family bit-score thresholds, not a single global E-value) with
// --rfam (heuristic filters tuned for whole-genome scans) and --clanin
// (keeps only the best-scoring member of each overlapping clan, e.g. the many
// SSU/LSU rRNA subfamilies).
process RUN_INFERNAL {
    tag        "${meta.id}"
    label      'infernal'
    storeDir   { "${params.outdir}/infernal/${hashBucketForType('infernal', meta.id)}" }

    input:
        tuple val(meta), path(genome)

    output:
        path("${meta.id}.rfam.tblout.gz"), emit: tblout

    script:
    """
    module load singularity
    SING_BINDS="--bind ${params.rfam_dbdir}:${params.rfam_dbdir}"
    singularity exec \${SING_BINDS} ${params.infernal_sif} cmscan \\
        --cut_ga --rfam --nohmmonly \\
        --clanin ${params.rfam_dbdir}/Rfam.clanin \\
        --tblout ${meta.id}.rfam.tblout \\
        --cpu ${task.cpus} \\
        ${params.rfam_dbdir}/Rfam.cm ${genome} > /dev/null
    pigz ${meta.id}.rfam.tblout
    """

    stub:
    """
    printf '#target name\\n' | gzip > ${meta.id}.rfam.tblout.gz
    """
}
