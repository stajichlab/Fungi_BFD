include { hashBucketForType } from '../../common/utils.nf'

// DeepTMHMM (deep-learning TM-topology caller, supersedes classic TMHMM --
// see RUN_TMHMM alongside this). No Lmod module or bioconda package exists;
// this reuses the same prebuilt academic-license SIF already used by
// ~/projects/nf/nf_funannotate1's DEEPTMHMM_ANNOTATION
// (/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif), so no
// second image/license copy is needed. GPU-capable (Ada arch) but also runs
// CPU-only, just slower -- gated by params.deeptmhmm_gpu like the source repo.
process RUN_DEEPTMHMM {
    tag        "${meta.locustag}"
    label      'deeptmhmm'
    storeDir   { "${params.outdir}/deeptmhmm/${hashBucketForType('deeptmhmm', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.deeptmhmm.gff3.gz"), emit: gff3

    script:
    """
    module load singularity
    SING="singularity exec ${params.deeptmhmm_gpu ? '--nv' : ''} ${params.deeptmhmm_sif}"
    OUTD=\$(mktemp -d)
    \${SING} bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta \${PWD}/${proteins} --output-dir \${OUTD}" || \\
    \${SING} bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta ${proteins} --output-dir \${OUTD}"
    pigz -c \$OUTD/TMRs.gff3 > ${meta.locustag}.deeptmhmm.gff3.gz
    rm -rf \$OUTD
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${meta.locustag}.deeptmhmm.gff3.gz
    """
}
