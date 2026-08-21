include { hashBucketForType } from '../../common/utils.nf'

// eggnog-mapper (diamond mode) against the EGGNOG_DATA_DIR eggnog.db +
// eggnog_proteins.dmnd already staged for funannotate's own emapper.py
// (params.eggnog_data_dir, see profile_funannotate.config) -- reused as-is so
// there is exactly one eggNOG copy on disk, not a second BFD-local one.
// Containerized (not `module load`) per the "less HPCC-dependent" direction:
// the biocontainers image needs only EGGNOG_DATA_DIR bind-mounted, no host
// module/license dependency, unlike phobius below.
process RUN_EGGNOG {
    tag        "${meta.locustag}"
    label      'eggnog'
    storeDir   { "${params.outdir}/eggnog/${hashBucketForType('eggnog', meta.locustag)}" }

    input:
        tuple val(meta), path(proteins)

    output:
        path("${meta.locustag}.emapper.annotations.gz"), emit: annotations

    script:
    """
    module load singularity
    SING_BINDS="--bind \${EGGNOG_DATA_DIR:-${params.eggnog_data_dir}}:\${EGGNOG_DATA_DIR:-${params.eggnog_data_dir}}"
    SING="singularity exec \${SING_BINDS} ${params.eggnog_sif}"
    \${SING} emapper.py -i ${proteins} \\
        --itype proteins -m diamond \\
        --data_dir ${params.eggnog_data_dir} \\
        --cpu ${task.cpus} \\
        --output ${meta.locustag} \\
        --output_dir . \\
        --no_file_comments \\
        --temp_dir \${TMPDIR:-.}
    pigz -c ${meta.locustag}.emapper.annotations > ${meta.locustag}.emapper.annotations.gz
    """

    stub:
    """
    printf '#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\tKEGG_Module\\tKEGG_Reaction\\tKEGG_rclass\\tBRITE\\tKEGG_TC\\tCAZy\\tBiGG_Reaction\\tPFAMs\\n' | gzip > ${meta.locustag}.emapper.annotations.gz
    """
}
