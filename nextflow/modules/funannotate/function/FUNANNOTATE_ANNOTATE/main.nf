process FUNANNOTATE_ANNOTATE {
    tag "$out"

    cpus   16
    memory '32 GB'
    time   '48h'

    input:
    tuple val(out), val(asmid), val(species), val(strain), val(locustag),
          val(busco_lineage), val(header_length), val(transl_table)

    output:
    tuple val(out), path("${out}.annotate.done"), emit: marker

    script:
    def antiSm    = file("${params.target}/${out}/antismash_local/${out}.gbk")
    def antiSmArg = antiSm.exists() ? "--antismash ${antiSm}" : ""
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load singularity
    # ── Containerized funannotate ──────────────────────────────────────────────
    # Same swap/rationale as FUNANNOTATE_TRAIN/main.nf: funannotate annotate runs
    # via \$SING (singularity exec \${params.funannotate_sif}) instead of
    # `module load funannotate`. No GeneMark/mysql involvement at this stage.
    # Binds cover target (annotate reads+writes there, incl. the optional
    # antismash_local GBK), augustus_config, funannotate_db, sbt_template,
    # eggnog_data_dir, and \$TMPDIR.
    unset -f which 2>/dev/null || true
    unset which_declare 2>/dev/null || true

    export AUGUSTUS_CONFIG_PATH=${params.augustus_config}
    export FUNANNOTATE_DB=${params.funannotate_db}
    # EggNOG-mapper's DB (used by funannotate annotate's auto-run of emapper.py --
    # see get_emapper_version()/annotate.py) is not baked into the funannotate-live
    # image (too large -- same reason the upstream Docker image omits it). Without
    # EGGNOG_DATA_DIR set, emapper.py --version prints a "couldn't find eggnog.db"
    # warning to stderr ahead of its actual version line, which funannotate's naive
    # stderr-prefix parse misreads as no version at all (LooseVersion(False) has no
    # .version attribute), crashing annotate outright. Matches the same
    # EGGNOG_DATA_DIR value the host funannotate/dev-1.9 module already sets.
    # Confirmed against a real beta.6 run, 2026-08-18.
    export EGGNOG_DATA_DIR=${params.eggnog_data_dir}
    TMPDIR=\${SCRATCH:-/tmp}
    SING_BINDS="--bind ${params.target}:${params.target},${params.augustus_config}:${params.augustus_config},${params.funannotate_db}:${params.funannotate_db},${params.sbt_template}:${params.sbt_template},${params.eggnog_data_dir}:${params.eggnog_data_dir},\$TMPDIR:\$TMPDIR"
    SING="singularity exec \${SING_BINDS} ${params.funannotate_sif}"

    \$SING funannotate annotate -i ${params.target}/${out} -o ${params.target}/${out} \\
        --species "${species}" --strain "${strain}" \\
        --busco_db ${busco_lineage} --rename ${locustag} \\
        --sbt ${params.sbt_template} \\
        --header_length ${header_length} \\
        ${antiSmArg} \\
        --cpu ${task.cpus} --tmpdir \$TMPDIR

    EXPECTED_GBK="${params.target}/${out}/annotate_results/${out}.gbk"
    if [ ! -f "\$EXPECTED_GBK" ]; then
        echo "ERROR: funannotate annotate did not produce expected GBK: \$EXPECTED_GBK" >&2
        exit 1
    fi
    touch ${out}.annotate.done
    """

    stub:
    """
    echo "[STUB] Would run funannotate annotate for ${out}"
    mkdir -p ${params.target}/${out}/annotate_results ${params.target}/${out}/annotate_misc
    touch ${params.target}/${out}/annotate_results/${out}.gbk
    touch ${out}.annotate.done
    """
}
