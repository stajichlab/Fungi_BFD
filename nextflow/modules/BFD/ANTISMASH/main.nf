include { hashBucketForType; globAnyOf } from '../../common/utils.nf'

// Secondary-metabolite biosynthetic gene cluster (BGC) prediction. Runs
// directly off BFD's existing gff3 + genome-fasta inputs via
// --genefinding-gff3 (no GBK conversion step, unlike
// ~/projects/nf/nf_funannotate1's ANTISMASH_RUN, which starts from a
// funannotate-produced .gbk). Reuses the already-populated shared DB at
// params.antismash_dbdir (/bigdata/stajichlab/shared/lib/antismash_db,
// 9.6 GB, confirmed populated 2026-08-20) and the prebuilt SIF in the shared
// singularity cache -- no new download needed.
//
// BATCHED: one task handles params.antismash_batch_size genomes (default 10) and
// runs params.antismash_batch_parallel of them at a time (default 2), each with
// task.cpus / antismash_batch_parallel CPUs. A genome takes 13-19 min in lite mode at
// 8 CPUs (6 genomes measured 2026-11-21), so 10 genomes at 2 in parallel is
// 5 rounds, about 65-95 min, which fits the 2 h limit of the 'short' queue.
// Per-genome time at other CPU counts and the slowest genomes are NOT measured.
// One failed genome fails the task and the whole batch is lost (storeDir stores
// outputs only after a task succeeds); the other genomes in the batch still run
// to completion before the task exits non-zero.
// The caller (BFD_FUNCTIONAL) passes only genomes with no stored regions table.
// storeDir treats a task as done when its output pattern matches anything in the
// store, so every output pattern names this batch's LOCUSTAGs.
//
// Two modes (params.antismash_mode, default 'lite'):
//   lite -- cluster detection + MIBiG known-cluster comparison only. Drops
//           --fullhmmer, --clusterhmmer (Pfam annotation of genes; they do not
//           detect clusters) and --pfam2go (GO, unused). Measured on 6 fungal
//           genomes 2026-11-21: 13-19 min/genome at 8 CPUs, against 109 min for
//           the full flag set on Af293, with identical regions on Af293 (44 = 44,
//           same coordinates and classes). The raw antiSMASH output (~240 MB per
//           genome) is NOT kept. Only small gzip TSV tables are:
//             <LOCUSTAG>.antismash.regions.tsv.gz  location, class, contig edge, core
//             <LOCUSTAG>.antismash.cds.tsv.gz      genes inside regions (locus tag,
//                                            coordinates, strand, kind, function)
//             <LOCUSTAG>.antismash.known.tsv.gz    top-5 MIBiG hits per region (with entry size, synteny)
//             <LOCUSTAG>.antismash.known_pairs.tsv.gz  our gene -> MIBiG gene pairs (identity, coverage)
//             <LOCUSTAG>.antismash.run.json        antiSMASH and MIBiG versions, mode
//   full -- the original flag set (--fullhmmer --clusterhmmer --cb-general
//           --pfam2go) and the raw JSON kept as <LOCUSTAG>.antismash.json.gz, plus
//           the same tables.
// Outputs are hash-bucketed by LOCUSTAG under params.outdir/antismash/<bucket>/ (same
// scheme as cazy and pfam_hmmscan; width 2). Each batch creates its genomes' bucket
// directories itself (same approach as BATCH_AA_FREQ).
// The reduction is bin/extract_antismash_subset.py (--minimal). A copy with
// more tables (domains, modules, A-domains) lives in
// DeltaGain_Fungi/analysis/bgc_subset/; keep the two in step.

process RUN_ANTISMASH {
    tag        "${metas instanceof List ? metas[0].locustag + '+' + (metas.size() - 1) : metas.locustag}"
    label      'antismash'
    storeDir   "${params.outdir}/antismash"

    input:
        tuple val(metas), path(gff3s), path(genomes)

    output:
        path("*/${globAnyOf((metas instanceof List ? metas : [metas]).collect { m -> m.locustag })}.antismash.regions.tsv.gz"),     emit: regions
        path("*/${globAnyOf((metas instanceof List ? metas : [metas]).collect { m -> m.locustag })}.antismash.cds.tsv.gz"),         emit: cds
        path("*/${globAnyOf((metas instanceof List ? metas : [metas]).collect { m -> m.locustag })}.antismash.known.tsv.gz"),       emit: known
        path("*/${globAnyOf((metas instanceof List ? metas : [metas]).collect { m -> m.locustag })}.antismash.known_pairs.tsv.gz"), emit: known_pairs   // our gene -> MIBiG gene hits
        path("*/${globAnyOf((metas instanceof List ? metas : [metas]).collect { m -> m.locustag })}.antismash.run.json"),           emit: run_info          // version, MIBiG data version, mode
        path("*/${globAnyOf((metas instanceof List ? metas : [metas]).collect { m -> m.locustag })}.antismash.json.gz"), optional: true, emit: json   // full mode only

    script:
    def metaList = metas instanceof List ? metas : [metas]
    def gffList  = gff3s instanceof List ? gff3s : [gff3s]
    def genList  = genomes instanceof List ? genomes : [genomes]
    def full = (params.antismash_mode ?: 'lite') == 'full'
    def par  = Math.max(1, (params.antismash_batch_parallel ?: 2) as int)
    def cpus_each = Math.max(1, task.cpus.intdiv(par))
    // One line per genome: LOCUSTAG GFF3 GENOME TRANSL_TABLE BUCKET (xargs -n 5).
    def rows = [metaList, gffList, genList].transpose().collect { m, g, d ->
        "${m.locustag} ${g.name} ${d.name} ${(m.transl_table ?: '1') as String} ${hashBucketForType('antismash', m.locustag as String)}"
    }.join('\n')
    // antiSMASH 8.0.4 has no genetic-code option and reads the GFF3 as code 1 (CUG = Leu).
    // For TRANSL_TABLE != 1 (samples.csv; 1,164 genomes use 12), tag every CDS line with
    // transl_table so translations match BFD. Measured on D. hansenii CBS767 (table 12):
    // 3,923 of 5,951 CDS differ without the tag, 0 with it. Table 26 is untested.
    // lite: no --cb-general. Measured 2026-11-22 on 3 genomes: 3.8x-6.9x faster, regions,
    // known and known_pairs identical. Only the unstored blast_top table needs it.
    // --no-*: smaller raw output in $SCRATCH (206 -> 102 MB on A1163); no time saving.
    def flags = full ? '--fullhmmer --clusterhmmer --cb-general --cb-knownclusters --pfam2go'
                     : '--cb-knownclusters --no-zip-output --no-summary-gbk --no-region-gbks'
    """
    module load apptainer
    export TMPDIR=\${SCRATCH:?SCRATCH is not set (expected the SLURM job scratch dir)}
    export MPLCONFIGDIR=\$TMPDIR/mpl; mkdir -p \$MPLCONFIGDIR
    export AS_DB='${params.antismash_dbdir}'
    export AS_SIF='${params.antismash_standalone_sif}'
    export AS_TAXON='${params.antismash_taxon}'
    export AS_CPUS=${cpus_each}
    export AS_MODE=${full ? 'full' : 'lite'}
    export AS_EXTRACT='${projectDir}/bin/extract_antismash_subset.py'
    export AS_ROOT=\$(mktemp -d -p \$TMPDIR)

    # Process one genome. Work in node-local scratch; only the reduced tables are
    # written to the task dir, under <bucket>/.
    one() {
        set -euo pipefail
        local lt=\$1 gff=\$2 gen=\$3 tt=\$4 bucket=\$5
        local w=\$AS_ROOT/\$lt
        mkdir -p \$w/in \$w/out
        local gffp=\${gff%.gz} genp=\${gen%.gz}
        # This tool may not read gzip, so inflate into scratch (FASTA/GFF3 may be .gz).
        if [[ \$gff == *.gz ]]; then pigz -dc "\$gff" > \$w/in/\$gffp; else cp -L "\$gff" \$w/in/\$gffp; fi
        if [[ \$gen == *.gz ]]; then pigz -dc "\$gen" > \$w/in/\$genp; else cp -L "\$gen" \$w/in/\$genp; fi
        if [ "\$tt" != "1" ]; then
            awk -F'\\t' -v OFS='\\t' -v t="\$tt" '\$0 !~ /^#/ && \$3=="CDS" {\$9=\$9";transl_table="t} {print}' \$w/in/\$gffp > \$w/in/tt.gff3
            mv \$w/in/tt.gff3 \$w/in/\$gffp
        fi
        apptainer exec --bind "\$AS_DB:\$AS_DB,\$TMPDIR:\$TMPDIR" "\$AS_SIF" \\
            antismash --taxon "\$AS_TAXON" \\
                --databases "\$AS_DB" \\
                --output-dir \$w/out/r \\
                --genefinding-gff3 \$w/in/\$gffp \\
                ${flags} \\
                -c \$AS_CPUS \\
                \$w/in/\$genp
        local J
        J=\$(ls \$w/out/r/*.json | head -1)
        # Fail loudly if antiSMASH wrote no JSON: an empty result must not look like "no clusters".
        test -s "\$J"
        python3 "\$AS_EXTRACT" "\$J" \$w/sub \$lt --minimal --mode \$AS_MODE
        mkdir -p "\$bucket"
        for t in regions cds known known_pairs; do mv \$w/sub/\$lt.\$t.tsv.gz "\$bucket/\$lt.antismash.\$t.tsv.gz"; done
        mv \$w/sub/\$lt.run.json "\$bucket/\$lt.antismash.run.json"
        ${full ? '[ "\$AS_MODE" = full ] && pigz -c "\$J" > "\$bucket/\$lt.antismash.json.gz"' : 'true'}
        rm -rf \$w
    }
    export -f one

    cat > genomes.list <<'LIST_EOF'
${rows}
LIST_EOF

    # xargs keeps running the other genomes when one fails and exits 123 at the end.
    xargs -a genomes.list -P ${par} -n 5 bash -c 'one "\$@"' _
    rm -rf \$AS_ROOT
    """

    stub:
    def metaList = metas instanceof List ? metas : [metas]
    """
    ${ metaList.collect { m ->
        def b = hashBucketForType('antismash', m.locustag as String)
        "mkdir -p ${b} && for t in regions cds known known_pairs; do printf '' | gzip > ${b}/${m.locustag}.antismash.\$t.tsv.gz; done && echo '{}' > ${b}/${m.locustag}.antismash.run.json"
    }.join('\n') }
    """
}
