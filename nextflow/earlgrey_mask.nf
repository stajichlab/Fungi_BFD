#!/usr/bin/env nextflow

/*
 * earlgrey_mask — curated repeat masking for large fungal genomes
 *
 * Runs EarlGrey once per species on the single best representative genome
 * (>cutoff_mb), then applies the resulting curated TE library with RepeatMasker to
 * every conspecific strain. Output is a soft-masked input_clean_genomes/<asmid>.masked.fasta
 * per strain — exactly the file funannotate.nf consumes (its tantan MASKREPEAT_TANTAN_RUN
 * uses storeDir on input_clean_genomes, so it is skipped wherever this file already exists).
 *
 * Runs independently of funannotate.nf. Submit from the project root.
 *
 * Usage:
 *   nextflow run nextflow/earlgrey_mask.nf \
 *       -c nextflow/nextflow.config -profile earlgrey -resume
 *
 *   # limit to N species for testing
 *   nextflow run nextflow/earlgrey_mask.nf \
 *       -c nextflow/nextflow.config -profile earlgrey --n_test 1 -resume
 *
 * Stub/dry-run:
 *   nextflow run nextflow/earlgrey_mask.nf \
 *       -c nextflow/nextflow.config -profile earlgrey -stub-run --n_test 2
 */

// ════════════════════════════════════════════════════════════════════════════
// PARAMETER DEFAULTS  (override via --param value)
// ════════════════════════════════════════════════════════════════════════════

params.genome_dir          = "${launchDir}/input_clean_genomes"
params.genome_suffix       = '.fa'                       // clean (unmasked) genome suffix
params.asm_stats           = "${launchDir}/tables/asm_stats.tsv.gz"
params.cutoff_mb           = 200                         // species qualifies if rep > this
params.repeat_taxon        = 'fungi'                     // EarlGrey -r RepeatMasker search term
params.earlgrey_version    = '7.2.6'
params.repeatmasker_version = '4.1.8'
params.outdir              = "${launchDir}/results/repeatlibrary"
params.masked_dir          = "${launchDir}/input_clean_genomes"   // where <asmid>.masked.fasta land

// ════════════════════════════════════════════════════════════════════════════
// PROCESSES
// ════════════════════════════════════════════════════════════════════════════

// ── SELECT_REPS ──────────────────────────────────────────────────────────────
// Joins samples.csv + asm_stats to pick the best representative per species and
// keep only species whose representative exceeds the size cutoff. storeDir caches
// the CSV in misc/ so re-runs skip this step (delete the file to refresh).
process SELECT_REPS {
    tag    'select'
    label  'select'

    storeDir "${launchDir}/misc"

    input:
        path samples
        path asm_stats

    output:
        path 'repeat_representatives.csv'

    script:
    """
    python ${launchDir}/scripts/select_repeat_representatives.py \\
        --samples ${samples} \\
        --asm-stats ${asm_stats} \\
        --genome-dir ${params.genome_dir} \\
        --genome-suffix ${params.genome_suffix} \\
        --cutoff-mb ${params.cutoff_mb} \\
        --output repeat_representatives.csv
    """

    stub:
    """
    printf 'SPECIES,REP_ASMID,REP_SIZE_MB,N_MEMBERS,MEMBER_ASMIDS\\n' > repeat_representatives.csv
    printf 'Stub species one,STUBASM_REP1,250.0,1,STUBASM_MEM1\\n'   >> repeat_representatives.csv
    printf 'Stub species two,STUBASM_REP2,300.0,0,\\n'              >> repeat_representatives.csv
    """
}

// ── EARLGREY_BUILD_LIB ───────────────────────────────────────────────────────
// Heavy: de-novo TE discovery + curation on the representative genome.
// -d yes produces the soft-masked representative genome; the curated consensus
// library (<species>-families.fa) is reused to mask conspecific strains.
//
// storeDir is our own results/repeatlibrary/<species> namespace — it gives
// cross-run persistence (the multi-day EarlGrey run is skipped if its products
// already exist) WITHOUT being confused by any pre-existing tantan mask sitting
// in input_clean_genomes. Delivery into input_clean_genomes is a separate
// overwrite step (DELIVER_MASK), so EarlGrey masks replace older tantan masks.
process EARLGREY_BUILD_LIB {
    tag    "${species} (${rep_asmid})"
    label  'earlgrey'

    storeDir { "${params.outdir}/${species.replaceAll(/[^A-Za-z0-9._-]+/, '_')}" }

    input:
        tuple val(species), val(rep_asmid), path(genome)

    output:
        tuple val(species), val(rep_asmid), path("${rep_asmid}.families.fa"), emit: lib
        tuple val(rep_asmid), path("${rep_asmid}.masked.fasta"), emit: masked
        path("*_RepeatLandscape"), emit: landscape, optional: true
        path("*_summaryFiles"),    emit: summary,   optional: true

    script:
    def sp_safe = species.replaceAll(/[^A-Za-z0-9._-]+/, '_')
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load earlgrey/${params.earlgrey_version}

    # rm -rf earlgrey_out
    mkdir -p earlgrey_out

    earlGrey \\
        -g ${genome} \\
        -s ${sp_safe} \\
        -o earlgrey_out \\
        -r ${params.repeat_taxon} \\
        -d yes \\
        -c yes \\
        -v yes \\
        -q yes \\
	-M 3200 \\
        -t ${task.cpus}

    # Locate EarlGrey products (paths are version-specific; glob defensively).
    lib=\$(find earlgrey_out -name '*-families.fa' | head -n1)
    soft=\$(find earlgrey_out -name '*.softmasked.fasta' | head -n1)

    if [ -z "\$lib" ] || [ -z "\$soft" ]; then
        echo "ERROR: EarlGrey outputs not found (lib=\$lib soft=\$soft)" >&2
        find earlgrey_out -maxdepth 3 -type f >&2
        exit 1
    fi

    cp "\$lib"  ${rep_asmid}.families.fa
    cp "\$soft" ${rep_asmid}.masked.fasta

    # Preserve the RepeatLandscape and summaryFiles report directories so storeDir
    # keeps them under results/repeatlibrary/<species>/ (named with the EarlGrey
    # species prefix). Copy to the task root; emitted via the optional path globs.
    for d in \$(find earlgrey_out -maxdepth 3 -type d \\( -name '*_RepeatLandscape' -o -name '*_summaryFiles' \\)); do
        cp -r "\$d" ./
    done
    """

    stub:
    """
    printf '>stub_family_1#Unknown\\nACGTACGTACGT\\n' > ${rep_asmid}.families.fa
    printf '>stub_${rep_asmid}\\nacgtACGTacgt\\n'      > ${rep_asmid}.masked.fasta

    sp_safe=\$(echo '${species}' | sed 's/[^A-Za-z0-9._-]\\+/_/g')
    mkdir -p "\${sp_safe}_RepeatLandscape" "\${sp_safe}_summaryFiles"
    printf 'stub landscape\\n' > "\${sp_safe}_RepeatLandscape/stub.txt"
    printf 'stub summary\\n'   > "\${sp_safe}_summaryFiles/stub.txt"
    """
}

// ── REPEATMASK_STRAIN ────────────────────────────────────────────────────────
// Applies the species curated library to a conspecific strain with RepeatMasker.
// -xsmall → soft-masked (lowercase) output, matching funannotate's expectation.
// storeDir caches in our results/<species>/strains namespace (persistent across
// runs); DELIVER_MASK handles the overwrite into input_clean_genomes.
process REPEATMASK_STRAIN {
    tag    "${asmid} (${species})"
    label  'repeatmask'

    storeDir { "${params.outdir}/${species.replaceAll(/[^A-Za-z0-9._-]+/, '_')}/strains" }

    input:
        tuple val(asmid), val(species), path(genome), path(library)

    output:
        tuple val(asmid), path("${asmid}.masked.fasta")

    script:
    """
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load RepeatMasker/${params.repeatmasker_version}

    mkdir -p rmask_out
    RepeatMasker \\
        -lib ${library} \\
        -xsmall \\
        -pa ${task.cpus} \\
        -dir rmask_out \\
        ${genome}

    # RepeatMasker writes <genome>.masked; if nothing was masked it may be absent,
    # in which case the (unmasked) input is the correct soft-masked result.
    if [ -f rmask_out/${genome}.masked ]; then
        cp rmask_out/${genome}.masked ${asmid}.masked.fasta
    else
        echo "WARN: no repeats masked for ${asmid}; using unmasked genome" >&2
        cp ${genome} ${asmid}.masked.fasta
    fi
    """

    stub:
    """
    printf '>stub_${asmid}\\nacgtACGTacgt\\n' > ${asmid}.masked.fasta
    """
}

// ── DELIVER_MASK ─────────────────────────────────────────────────────────────
// Copies a finished soft-masked genome into input_clean_genomes, OVERWRITING any
// pre-existing tantan mask (overwrite: true). Kept separate from the masking
// processes so their storeDir cache is never confused by foreign masks, and so a
// stub run can be redirected to a sandbox dir instead of the real genome dir.
process DELIVER_MASK {
    tag    "${asmid}"
    label  'deliver'

    publishDir path: { workflow.stubRun ? "${launchDir}/work/stub_masked" : params.masked_dir },
               mode: 'copy', overwrite: true

    input:
        tuple val(asmid), path(masked)

    output:
        path("${asmid}.masked.fasta")

    script:
    'true'

    stub:
    'true'
}

// ════════════════════════════════════════════════════════════════════════════
// MAIN WORKFLOW
// ════════════════════════════════════════════════════════════════════════════

workflow {

    // ── Select representatives ────────────────────────────────────────────────
    def reps = SELECT_REPS(
        file(params.samples,   glob: false),
        file(params.asm_stats, glob: false),
    )

    // ── Per-species records (n_test limits *species*) ─────────────────────────
    def records = reps
        .splitCsv(header: true)
        .take(params.n_test > 0 ? params.n_test as int : -1)
        .map { row ->
            tuple(row.SPECIES?.trim(), row.REP_ASMID?.trim(), (row.MEMBER_ASMIDS ?: '').trim())
        }

    // ── EarlGrey library build on the representative ──────────────────────────
    def build_in = records.map { species, rep_asmid, _members ->
        def g = file("${params.genome_dir}/${rep_asmid}${params.genome_suffix}", glob: false)
        if (!g.exists() && !workflow.stubRun) {
            log.warn "Skipping ${species}: representative genome not found at ${g}"
            return null
        }
        tuple(species, rep_asmid, g)
    }.filter { it != null }

    EARLGREY_BUILD_LIB(build_in)

    // ── Expand conspecific members and pair each with its species library ─────
    def members_ch = records.flatMap { species, _rep, members ->
        (members ? members.split(';') : [])
            .findAll { it?.trim() }
            .collect { asm -> tuple(species, asm.trim()) }
    }

    def lib_by_species = EARLGREY_BUILD_LIB.out.lib
        .map { species, _rep, lib -> tuple(species, lib) }

    def mask_in = members_ch
        .combine(lib_by_species, by: 0)     // tuple(species, asmid, library)
        .map { species, asmid, library ->
            def g = file("${params.genome_dir}/${asmid}${params.genome_suffix}", glob: false)
            if (!g.exists() && !workflow.stubRun) {
                log.warn "Skipping member ${asmid} (${species}): genome not found at ${g}"
                return null
            }
            tuple(asmid, species, g, library)
        }
        .filter { it != null }

    REPEATMASK_STRAIN(mask_in)

    // ── Deliver all soft-masked genomes into input_clean_genomes ──────────────
    // Both the EarlGrey-masked representative and the RepeatMasked members are
    // copied in (overwriting any older tantan mask). Each is tuple(asmid, masked).
    def masked_ch = EARLGREY_BUILD_LIB.out.masked
        .mix(REPEATMASK_STRAIN.out)

    DELIVER_MASK(masked_ch)
}
