include { hashBucketForType } from '../../common/utils.nf'

// CAZyme Gene Cluster (CGC) / Polysaccharide Utilization Locus (PUL)
// prediction -- item #7 (dbCAN CGC/PUL) from the functional-annotation
// followups. dbcanlight (RUN_CAZY's tool) has no cgc mode, so this runs the
// full run_dbcan CLI instead, via dbCAN 5.2.9 (params.dbcan_cgc_sif -- see
// profile_BFD.config for why this isn't the originally-suggested
// haidyi/run_dbcan:latest image), against the separate full DB set
// SETUP_DBCAN_DB stages.
//
// CLI verified 2026-08-20 by pulling the real image and running
// `run_dbcan easy_CGC --help` -- flags below (--input_raw_data/--mode/
// --input_gff/--gff_type/--output_dir/--db_dir/--threads) all confirmed
// present, and match nf-core/funcscan's RUNDBCAN_EASYCGC module (same image
// version) exactly. --gff_type NCBI_prok is required when --mode=protein;
// funannotate's own gff3 output is NCBI-style (locus_tag'd gene/mRNA/CDS
// features), which is the closest of run_dbcan's dialect labels -- not
// independently verified against a real run_dbcan parse, so double-check
// cgc_standard_out.tsv on the first real genome.
process RUN_CAZY_CGC {
    tag        "${meta.id}"
    label      'cazy_cgc'
    storeDir   { "${params.outdir}/cazy_cgc/${hashBucketForType('cazy_cgc', meta.id)}" }

    input:
        tuple val(meta), path(gff3), path(proteins)
        path(dbcan_db)

    output:
        path("${meta.id}.cgc.tsv.gz"), emit: cgc

    script:
    """
    module load apptainer
    export TMPDIR=\${SCRATCH:-/tmp}
    SING_BINDS="--bind \${PWD}:\${PWD},${params.dbcan_dbdir}:${params.dbcan_dbdir},\$TMPDIR:\$TMPDIR"
    OUTD=\$(mktemp -d)
    apptainer exec \${SING_BINDS} ${params.dbcan_cgc_sif} run_dbcan easy_CGC \\
        --input_raw_data ${proteins} \\
        --mode protein \\
        --input_gff ${gff3} \\
        --gff_type NCBI_prok \\
        --output_dir \${OUTD} \\
        --db_dir ${params.dbcan_dbdir} \\
        --threads ${task.cpus}
    pigz -c \${OUTD}/cgc_standard_out.tsv > ${meta.id}.cgc.tsv.gz
    rm -rf \${OUTD}
    """

    stub:
    """
    printf 'CGC_id\\n' | gzip > ${meta.id}.cgc.tsv.gz
    """
}
