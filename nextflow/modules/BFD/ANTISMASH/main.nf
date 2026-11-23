include { hashBucketForType; plainName; gzInflate } from '../../common/utils.nf'

// Secondary-metabolite biosynthetic gene cluster (BGC) prediction. Runs
// directly off BFD's existing gff3 + genome-fasta inputs via
// --genefinding-gff3 (no GBK conversion step, unlike
// ~/projects/nf/nf_funannotate1's ANTISMASH_RUN, which starts from a
// funannotate-produced .gbk). Reuses the already-populated shared DB at
// params.antismash_dbdir (/bigdata/stajichlab/shared/lib/antismash_db,
// 9.6 GB, confirmed populated 2026-08-20) and the prebuilt SIF in the shared
// singularity cache -- no new download needed.
//
// Two modes (params.antismash_mode, default 'lite'):
//   lite -- cluster detection + MIBiG known-cluster comparison only. Drops
//           --fullhmmer, --clusterhmmer (Pfam annotation of genes; they do not
//           detect clusters) and --pfam2go (GO, unused). Measured on 6 fungal
//           genomes 2026-11-21: 13-19 min/genome at 8 CPUs, against 109 min for
//           the full flag set on Af293, with identical regions on Af293 (44 = 44,
//           same coordinates and classes). The raw antiSMASH output (~240 MB per
//           genome) is NOT kept. Only three small gzip TSV tables are:
//             <LOCUSTAG>.antismash.regions.tsv.gz  location, class, contig edge, core
//             <LOCUSTAG>.antismash.cds.tsv.gz      genes inside regions (locus tag,
//                                            coordinates, strand, kind, function)
//             <LOCUSTAG>.antismash.known.tsv.gz    top-5 MIBiG hits per region (with entry size, synteny)
//             <LOCUSTAG>.antismash.known_pairs.tsv.gz  our gene -> MIBiG gene pairs (identity, coverage)
//             <LOCUSTAG>.antismash.run.json        antiSMASH and MIBiG versions, mode
//           ~2 MB per genome for all of the subset tables; the three minimal
//           tables are smaller (UNVERIFIED size for the minimal set alone).
//   full -- the original flag set (--fullhmmer --clusterhmmer --cb-general
//           --pfam2go) and the raw JSON kept as <LOCUSTAG>.antismash.json.gz, plus
//           the same three tables.
// Outputs are hash-bucketed by LOCUSTAG under params.outdir/antismash/<bucket>/ (same
// scheme as cazy and pfam_hmmscan; width 2).
// The reduction is bin/extract_antismash_subset.py (--minimal). A copy with
// more tables (domains, modules, A-domains) lives in
// DeltaGain_Fungi/analysis/bgc_subset/; keep the two in step.

process RUN_ANTISMASH {
    tag        "${meta.locustag}"
    label      'antismash'
    storeDir   { "${params.outdir}/antismash/${hashBucketForType('antismash', meta.locustag as String)}" }

    input:
        tuple val(meta), path(gff3), path(genome)

    output:
        path("${meta.locustag}.antismash.regions.tsv.gz"), emit: regions
        path("${meta.locustag}.antismash.cds.tsv.gz"),     emit: cds
        path("${meta.locustag}.antismash.known.tsv.gz"),   emit: known
        path("${meta.locustag}.antismash.known_pairs.tsv.gz"), emit: known_pairs   // our gene -> MIBiG gene hits
        path("${meta.locustag}.antismash.run.json"),       emit: run_info          // version, MIBiG data version, mode
        path("${meta.locustag}.antismash.json.gz"), optional: true, emit: json   // full mode only

    script:
    def gff3_plain = plainName(gff3)
    // The genome may be <id>.scaffolds.fa.gz (bgzip, FUNANNOTATE_PREDICT since 2026-11-07).
    // Inflate it into the task dir: this tool may not read gzip.
    def genome_plain = plainName(genome)
    def genome_rm = genome.name.endsWith('.gz') ? "rm -f '${genome_plain}'" : ''
    def full = (params.antismash_mode ?: 'lite') == 'full'
    // antiSMASH 8.0.4 has no genetic-code option and reads the GFF3 as code 1 (CUG = Leu).
    // For TRANSL_TABLE != 1 (samples.csv; 1,164 genomes use 12), tag every CDS line with
    // transl_table so translations match BFD. Measured on D. hansenii CBS767 (table 12):
    // 3,923 of 5,951 CDS differ without the tag, 0 with it. Table 26 is untested.
    def tt = (meta.transl_table ?: '1') as String
    def tt_cmd = tt == '1' ? '' :
        "awk -F'\\t' -v OFS='\\t' -v t=${tt} '\$0 !~ /^#/ && \$3==\"CDS\" {\$9=\$9\";transl_table=\"t} {print}' ${gff3_plain} > tt.gff3 && mv tt.gff3 ${gff3_plain}"
    // lite: no --cb-general. Measured 2026-11-22 on 3 genomes: 3.8x-6.9x faster, regions,
    // known and known_pairs identical. Only the unstored blast_top table needs it.
    // --no-*: smaller raw output in $SCRATCH (206 -> 102 MB on A1163); no time saving.
    def flags = full ? '--fullhmmer --clusterhmmer --cb-general --cb-knownclusters --pfam2go'
                     : '--cb-knownclusters --no-zip-output --no-summary-gbk --no-region-gbks'
    """
    ${gzInflate(gff3)}
    ${gzInflate(genome)}
    ${tt_cmd}
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    export MPLCONFIGDIR=\$TMPDIR/mpl; mkdir -p \$MPLCONFIGDIR
    # ${gff3_plain}/${genome} may be symlinks Nextflow staged from outside \$PWD (e.g.
    # a shared genome_annotation dir elsewhere under /bigdata); bind their real
    # parent dirs too, or apptainer won't see the symlink targets inside the
    # container.
    GFF3_REAL_DIR=\$(dirname \$(readlink -f ${gff3_plain}))
    GENOME_REAL_DIR=\$(dirname \$(readlink -f ${genome_plain}))
    SING_BINDS="--bind ${params.antismash_dbdir}:${params.antismash_dbdir},\${GFF3_REAL_DIR}:\${GFF3_REAL_DIR},\${GENOME_REAL_DIR}:\${GENOME_REAL_DIR},\${PWD}:\${PWD},\$TMPDIR:\$TMPDIR"
    # Work in node-local scratch; only the reduced tables are written to the task dir.
    OUTD=\$(mktemp -d -p \$TMPDIR)
    apptainer exec \${SING_BINDS} ${params.antismash_standalone_sif} \
        antismash --taxon ${params.antismash_taxon} \
            --databases ${params.antismash_dbdir} \
            --output-dir \${OUTD} \
            --genefinding-gff3 ${gff3_plain} \
            ${flags} \
            -c ${task.cpus} \
            ${genome_plain}
    J=\$(ls \${OUTD}/*.json | head -1)
    # Fail loudly if antiSMASH wrote no JSON: an empty result must not look like "no clusters".
    test -s "\$J"
    python3 ${projectDir}/bin/extract_antismash_subset.py "\$J" \${OUTD}/sub ${meta.locustag} --minimal --mode ${full ? "full" : "lite"}
    for t in regions cds known known_pairs; do mv \${OUTD}/sub/${meta.locustag}.\$t.tsv.gz ${meta.locustag}.antismash.\$t.tsv.gz; done
    mv \${OUTD}/sub/${meta.locustag}.run.json ${meta.locustag}.antismash.run.json
    ${full ? "pigz -c \$J > ${meta.locustag}.antismash.json.gz" : ''}
    rm -rf \${OUTD}
    ${genome_rm}
    """

    stub:
    """
    for t in regions cds known known_pairs; do printf '' | gzip > ${meta.locustag}.antismash.\$t.tsv.gz; done
    echo '{}' > ${meta.locustag}.antismash.run.json
    """
}
