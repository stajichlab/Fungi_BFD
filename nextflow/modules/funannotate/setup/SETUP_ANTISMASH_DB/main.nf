// Populate the antiSMASH reference database once; storeDir caches it at
// params.antismash_databases so subsequent runs skip this entirely, exactly
// like SETUP_TAXONDB.
//
// ── DEFAULT BEHAVIOR: reuses the HPCC module's DB, downloads nothing ──────
// params.antismash_databases defaults (see profile_funannotate.config) to
// the path `module load antismash/8.0.4` resolves to on this cluster -- an
// already-populated (~9.5 GB) database dir owned by pkgadmin. storeDir only
// checks whether the declared output (README) already exists there; since it
// already does, THIS TASK NEVER RUNS by default -- ANTISMASH_RUN just binds
// that directory straight into the container. No download, no write access
// to that dir needed (world-readable is enough).
//
// This task only actually executes (and only then downloads anything) when
// params.antismash_databases is overridden to a not-yet-populated directory
// -- e.g. for a portable run off this cluster, per profile_funannotate.config's
// antismash_databases comment.
//
// ── IMPORTANT ASSUMPTION when reusing a pre-existing directory ────────────
// storeDir does NOT verify database version or schema, only that README
// exists. Reusing the module's DB is only safe if the antiSMASH build inside
// params.antismash_sif is the exact same version+build as whatever populated
// that directory -- antiSMASH's database schema is not guaranteed compatible
// across versions/builds. params.antismash_sif is pinned to
// quay.io/biocontainers/antismash:8.0.4--pyhdfd78af_1 specifically because it
// is the exact bioconda build backing the HPCC `antismash/8.0.4` module
// (confirmed via `conda list` in that module env: antismash 8.0.4
// pyhdfd78af_1). If either pin ever moves, re-verify the match (or point
// antismash_databases at a fresh, not-yet-populated dir so this task
// downloads a version-matched copy instead).
process SETUP_ANTISMASH_DB {
    storeDir params.antismash_databases

    cpus   2
    memory '8 GB'
    time   '8h'

    output:
    path "README", emit: ready

    script:
    """
    set -euo pipefail
    source /etc/profile.d/modules.sh 2>/dev/null || true
    module load singularity
    singularity exec ${params.antismash_sif} download-antismash-databases --database-dir \$(pwd)
    """

    stub:
    """
    touch README
    """
}
