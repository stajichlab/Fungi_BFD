include { hashBucketForType; plainName; gzInflate } from '../../common/utils.nf'

// DeepTMHMM (deep-learning TM-topology caller, supersedes classic TMHMM --
// see RUN_TMHMM alongside this). No Lmod module or bioconda package exists;
// this reuses the same prebuilt academic-license SIF already used by
// ~/projects/nf/nf_funannotate1's DEEPTMHMM_ANNOTATION
// (/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif), so no
// second image/license copy is needed. GPU-capable (Ada arch) but also runs
// CPU-only, just slower -- gated by params.deeptmhmm_gpu like the source repo.
// The faster driver predict_fast.py lives only in the private derived image
// DeepTMHMM-1.0-fast.sif (license: never publish that image or its recipe);
// select it with deeptmhmm_sif + deeptmhmm_script = 'predict_fast.py'.
process RUN_DEEPTMHMM {
    tag        "${meta.locustag}"
    label      'deeptmhmm'
    storeDir   { "${params.outdir}/deeptmhmm/${hashBucketForType('deeptmhmm', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.deeptmhmm.gff3.gz"), emit: gff3

    script:
    def proteins_plain = plainName(proteins)
    // predict.py (stock image) or predict_fast.py (private DeepTMHMM-1.0-fast.sif:
    // batched ESM-1b + batched topology models, ~4.5x faster on GPU at
    // --topo-batch 16, measured 2026-11-05). predict_fast.py takes its thread
    // count from NCPU; stock predict.py always uses every CPU on the node.
    def fast = params.deeptmhmm_script == 'predict_fast.py'
    def fastArgs = fast ? "--topo-batch ${params.deeptmhmm_topo_batch}" : ''
    """
    set -euo pipefail
    ${gzInflate(proteins)}
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    # ${proteins_plain} may be a symlink Nextflow staged from outside \$PWD (e.g. a
    # shared genome_annotation dir elsewhere under /bigdata); bind its real
    # parent dir too, or apptainer won't see the symlink target inside the
    # container.
    PROT_REAL_DIR=\$(dirname \$(readlink -f ${proteins_plain}))
    SING_BINDS="--bind \${PWD}:\${PWD},\${PROT_REAL_DIR}:\${PROT_REAL_DIR},\$TMPDIR:\$TMPDIR"
    SING="apptainer exec ${params.deeptmhmm_gpu ? '--nv' : ''} \${SING_BINDS} ${params.deeptmhmm_sif}"
    # predict.py exits 1 when --output-dir already exists, so pass a path that
    # does not exist yet. (Before 2026-11-05 this used OUTD=\$(mktemp -d): every
    # run failed and left an empty .deeptmhmm.gff3.gz in the storeDir.)
    # Embeddings (~5 KB per residue with predict.py) go under \$TMPDIR = \$SCRATCH.
    WORKD=\$(mktemp -d)
    OUTD=\$WORKD/out
    \${SING} env NCPU=${task.cpus} bash -c "cd /opt/deeptmhmm && python3 ${params.deeptmhmm_script} --fasta \${PWD}/${proteins_plain} --output-dir \${OUTD} ${fastArgs}"
    if [ ! -s \$OUTD/TMRs.gff3 ]; then
        echo "ERROR: DeepTMHMM wrote no TMRs.gff3 for ${meta.locustag}" >&2
        exit 1
    fi
    pigz -c \$OUTD/TMRs.gff3 > ${meta.locustag}.deeptmhmm.gff3.gz
    rm -rf \$WORKD
    """

    stub:
    """
    printf '##gff-version 3\\n' | gzip > ${meta.locustag}.deeptmhmm.gff3.gz
    """
}
