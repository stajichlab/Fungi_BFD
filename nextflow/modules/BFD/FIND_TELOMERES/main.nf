include { hashBucketForType } from '../../common/utils.nf'

process FIND_TELOMERES {
    label    'telomeres'
    tag      "${meta.asmid}"
    storeDir { "${params.genome_stats_outdir}/telomeres/${hashBucketForType('telomeres', meta.asmid)}" }

    input:
    tuple val(meta), path(genome)

    output:
    path "${meta.asmid}.telomeres.tsv.gz", emit: tsv

    script:
    def patterns_arg = params.telomere_patterns ? "-p ${params.telomere_patterns.join(' ')}" : ''
    def fuzzy_arg    = params.telomere_fuzzy.toBoolean() ? '--fuzzy' : ''
    def both_ends    = params.telomere_both_ends.toBoolean() ? '--both-ends' : ''
    def allow_internal = params.telomere_allow_internal.toBoolean() ? '--allow-internal' : ''
    """
    module load biopython
    python3 ${projectDir}/bin/find_telomeres.py \
        ${patterns_arg} \
        -n ${params.telomere_min_repeats} \
        -l ${params.telomere_min_length} \
        -w ${params.telomere_flank_window} \
        -s ${params.telomere_search_window} \
        ${fuzzy_arg} \
        --max-mismatch ${params.telomere_max_mismatch} \
        --max-indel ${params.telomere_max_indel} \
        --terminal-tolerance ${params.telomere_terminal_tolerance} \
        ${both_ends} \
        ${allow_internal} \
        ${genome} \
        -o ${meta.asmid}.telomeres.tsv
    gzip -n ${meta.asmid}.telomeres.tsv
    """

    stub:
    """
    cat > ${meta.asmid}.telomeres.tsv <<'EOF'
scaffold	end	strand	monomer	repeat_count	tract_length	start	end_coord	terminal	distance_to_end	tract_seq	flank_seq
scaffold_1	5prime	+	TTAGGG	12	72	0	72	True	0	ttagggttagggttagggttagggttagggttagggttagggttagggttagggttagggttagggttaggg	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
scaffold_1	3prime	-	TTAGGG	10	60	394	454	True	0	ccctaaccctaaccctaaccctaaccctaaccctaa	ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
EOF
    gzip -n ${meta.asmid}.telomeres.tsv
    """
}
