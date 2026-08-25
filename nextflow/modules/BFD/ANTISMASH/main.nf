include { hashBucketForType } from '../../common/utils.nf'

// Secondary-metabolite biosynthetic gene cluster (BGC) prediction. Runs
// directly off BFD's existing gff3 + genome-fasta inputs via
// --genefinding-gff3 (no GBK conversion step, unlike
// ~/projects/nf/nf_funannotate1's ANTISMASH_RUN, which starts from a
// funannotate-produced .gbk). Reuses the already-populated shared DB at
// params.antismash_dbdir (/bigdata/stajichlab/shared/lib/antismash_db,
// 9.6 GB, confirmed populated 2026-08-20) and the prebuilt SIF in the shared
// singularity cache -- no new download needed.
process RUN_ANTISMASH {
    tag        "${meta.id}"
    label      'antismash'
    storeDir   { "${params.outdir}/antismash/${hashBucketForType('antismash', meta.id)}" }

    input:
        tuple val(meta), path(gff3), path(genome)

    output:
        path("${meta.id}.antismash.json.gz"), emit: json

    script:
    """
    module load apptainer
    export TMPDIR=\${SCRATCH:-/tmp}
    SING_BINDS="--bind ${params.antismash_dbdir}:${params.antismash_dbdir},\${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    OUTD=\$(mktemp -d)
    apptainer exec \${SING_BINDS} ${params.antismash_sif} \\
        antismash --taxon ${params.antismash_taxon} \\
            --databases ${params.antismash_dbdir} \\
            --output-dir \${OUTD} \\
            --genefinding-gff3 ${gff3} \\
            --fullhmmer --clusterhmmer --cb-general --pfam2go \\
            -c ${task.cpus} \\
            ${genome}
    pigz -c \${OUTD}/*.json > ${meta.id}.antismash.json.gz
    rm -rf \${OUTD}
    """

    stub:
    """
    printf '{}' | gzip > ${meta.id}.antismash.json.gz
    """
}
