process SELECT_REPS {
    tag    'select'
    label  'earlgrey_select'

    storeDir "${launchDir}/misc"

    input:
        path samples
        path asm_stats_parquet
        path busco_genome_parquet

    output:
        path 'repeat_representatives.csv'

    script:
    def suppress_arg = params.suppress_list ? "--suppress-list ${params.suppress_list}" : ""
    // Join the merged BFD database tables (same asm_stats.parquet +
    // busco_genome.parquet MERGE_ASM_STATS / MERGE_BUSCO_GENOME publish for
    // pick_representative_strain.py) into the flat TSV
    // select_repeat_representatives.py ranks candidates from -- BUSCO
    // completeness first, N50 second, fewest contigs as final tiebreak,
    // mirroring PICK_REPRESENTATIVE_STRAIN's funannotate-side selection.
    """
    module load duckdb 2>/dev/null || true
    duckdb -c "COPY (SELECT a.ASMID, a.total_length_bp, a.N50_bp, a.contig_count, COALESCE(CAST(b.complete_pct AS DOUBLE), -1) AS complete_pct FROM read_parquet('${asm_stats_parquet}') a LEFT JOIN read_parquet('${busco_genome_parquet}') b USING (ASMID)) TO 'asm_stats_busco.tsv' (FORMAT CSV, DELIMITER '\\t', HEADER TRUE);"
    python ${projectDir}/bin/select_repeat_representatives.py \
        --samples ${samples} \
        --asm-stats asm_stats_busco.tsv \
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
