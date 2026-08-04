process SELECT_REPS {
    tag    'select'
    label  'earlgrey_select'

    storeDir "${launchDir}/misc"

    input:
        path samples
        path asm_stats

    output:
        path 'repeat_representatives.csv'

    script:
    def suppress_arg = params.suppress_list ? "--suppress-list ${params.suppress_list}" : ""
    """
    python ${projectDir}/bin/select_repeat_representatives.py \
        --samples ${samples} \
        --asm-stats ${asm_stats} \
        --genome-dir ${params.genome_dir} \
        --genome-suffix ${params.genome_suffix} \
        --cutoff-mb ${params.cutoff_mb} \
        ${suppress_arg} \
        --output repeat_representatives.csv
    """

    stub:
    """
    printf 'SPECIES,REP_ASMID,REP_SIZE_MB,N_MEMBERS,MEMBER_ASMIDS\\n' > repeat_representatives.csv
    printf 'Stub species one,STUBASM_REP1,250.0,1,STUBASM_MEM1\\n'   >> repeat_representatives.csv
    printf 'Stub species two,STUBASM_REP2,300.0,0,\\n'              >> repeat_representatives.csv
    """
}
