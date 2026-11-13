process REPEATMASK_STRAIN {
    tag    "${asmid} (${species})"
    label  'earlgrey_repeatmask'

    storeDir { "${params.outdir}/${species.replaceAll(/[^A-Za-z0-9._-]+/, '_')}/strains" }

    input:
        tuple val(asmid), val(species), path(genome), path(library)

    output:
        tuple val(asmid), path("${asmid}.masked.fasta.gz")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load RepeatMasker/${params.repeatmasker_version}

    GENOME_IN="${genome}"
    case "${genome}" in
        *.gz) echo "[INFO] Inflating compressed genome ${genome}"; gzip -dc "${genome}" > genome_input.fa; GENOME_IN=genome_input.fa ;;
    esac

    mkdir -p rmask_out
    # -engine rmblast: since 2026-11-11 the RepeatMasker/4.1.8 module defaults to the HMMER
    # engine, which fails in hmmpress on a FASTA library (exit 2, no output). A custom FASTA
    # library needs rmblast.
    RepeatMasker \
        -engine rmblast \
        -lib ${library} \
        -xsmall \
        -pa ${task.cpus} \
        -dir rmask_out \
        "\$GENOME_IN" \
        || { echo "ERROR: RepeatMasker failed for ${asmid} (exit \$?)" >&2; exit 1; }

    # RepeatMasker writes <genome>.masked even when no repeat matches (checked with 4.1.8), so a
    # missing file means a failure. Do not deliver an unmasked genome as a masked one.
    MASKED="rmask_out/\$(basename "\$GENOME_IN").masked"
    if [ ! -s "\$MASKED" ]; then
        echo "ERROR: RepeatMasker made no masked file for ${asmid}" >&2
        exit 1
    fi
    gzip -c "\$MASKED" > ${asmid}.masked.fasta.gz
    """

    stub:
    """
    printf '>stub_${asmid}\\nacgtACGTacgt\\n' | gzip -c > ${asmid}.masked.fasta.gz
    """
}
