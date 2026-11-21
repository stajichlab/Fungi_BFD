#!/usr/bin/bash -l
#SBATCH -p epyc -N 1 -n 1 -c 16 --mem 64G --time 2:00:00
#SBATCH -J rebuild_pipeline_duckdb
#SBATCH -o logs/rebuild_pipeline_duckdb.%j.log
#
# Rebuild a run directory's db/BFD.duckdb from its tables/*.parquet, using the
# pipeline builder (nextflow/bin/build_BFD_duckDB.sh). It does what the
# Nextflow BUILD_DUCKDB process does, without Nextflow, so no MERGE_* step
# runs. Use it when the builder or the SQL changed but the data did not.
#
# Run it from the run directory, naming this script by its path in the
# checkout you want to build with:
#
#     cd /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/BFD_function
#     sbatch /path/to/Fungi_BFD/scripts/rebuild_pipeline_duckdb.sh
#
# Builds into a staging directory and moves the result into place only at the
# end, so a failed build leaves the previous database untouched.
#
# Do NOT run while a MERGE_* step is rewriting tables/: publishDir copies in
# the background, so a parquet file can exist at its final path while still
# being written. The Nextflow guard below covers a run in this directory; a
# merge launched from elsewhere is not detectable here.
#
# Lived untracked in Fungi_BFD_runs/BFD_function/ until 2026-11-19. It is
# generic code -- the run directory it acts on is $PWD and the checkout it
# builds with is $PIPELINE -- so it belongs under review like anything else
# that defines how a build happens.
set -euo pipefail

RUNDIR=$PWD
case "$RUNDIR" in
    */BFD_function) ;;
    *) echo "ERROR: run from a run directory ending in /BFD_function (got $RUNDIR)" >&2; exit 1 ;;
esac

# Prefer the checkout this script was invoked from, so the build matches the
# code under review. sbatch runs a spooled copy, where BASH_SOURCE points at
# the spool rather than the repo, so fall back to the shared checkout and
# accept PIPELINE as an explicit override.
if [ -z "${PIPELINE:-}" ]; then
    _self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)"
    _candidate="${_self_dir%/scripts}"
    if [ -n "$_candidate" ] && [ -f "$_candidate/nextflow/bin/build_BFD_duckDB.sh" ]; then
        PIPELINE="$_candidate"
    else
        PIPELINE=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
    fi
fi
echo "pipeline checkout: $PIPELINE"

BUILDER=$PIPELINE/nextflow/bin/build_BFD_duckDB.sh
[ -f "$BUILDER" ] || { echo "ERROR: missing $BUILDER" >&2; exit 1; }

# Staleness guards: files the builder cannot run without. These are structural
# -- either the checkout has them or the build is impossible -- which is why
# they are checked here and claims about the CONTENT of the SQL are not. The
# GPI-anchor filter used to be grepped out of $BUILDER here; it is now asserted
# against the database that was actually built, in verify_BFD_duckDB.py, so the
# check survives reformatting and also covers a database built elsewhere.
# PR #81 moved the compat and helper SQL into the pipeline. Without them the
# builder produces far fewer views and silently lacks every v_genome_* helper.
# The labels and traits layer lives in BFD_dashboard (2026-11-20).
for f in sql/bfd_compat_views.sql sql/bfd_extra_views.sql \
         nextflow/bin/filter_sql_for_db.py nextflow/bin/verify_BFD_duckDB.py; do
    [ -s "$PIPELINE/$f" ] \
        || { echo "ERROR: $PIPELINE/$f missing or empty (Fungi_BFD PR #81); update the checkout" >&2; exit 1; }
done

if pgrep -u "$USER" -f "nextflow.*main.nf" | xargs -r -I{} readlink /proc/{}/cwd | grep -qx "$RUNDIR"; then
    echo "ERROR: a Nextflow run is live in $RUNDIR; it may be rewriting tables/" >&2
    exit 1
fi

[ -d "$RUNDIR/tables" ] || { echo "ERROR: no tables/ in $RUNDIR" >&2; exit 1; }

module load duckdb
mkdir -p "$RUNDIR/db" logs
STAGE=$RUNDIR/db/.stage.${SLURM_JOB_ID:-manual}
rm -rf "$STAGE"; mkdir -p "$STAGE"
trap 'rm -rf "$STAGE"' EXIT

# The database being replaced, so the builder's verifier can refuse a build
# that has lost views rather than leaving the regression for a consumer.
PREV_DB=""
[ -f "$RUNDIR/db/BFD.duckdb" ] && PREV_DB="$RUNDIR/db/BFD.duckdb"

echo "Build start $(date -Is) rundir=$RUNDIR stage=$STAGE"
TABLES=$RUNDIR/tables DBDIR=$STAGE \
BIN_DIR=$PIPELINE/nextflow/bin SCHEMA=$PIPELINE/sql/table_schema.json \
PREV_DB="$PREV_DB" \
    bash "$BUILDER"

[ -s "$STAGE/BFD.duckdb" ] || { echo "ERROR: no database built" >&2; exit 1; }

# build_BFD_duckDB.sh runs verify_BFD_duckDB.py itself and exits non-zero on a
# problem, so reaching here means the build was checked. Print the shape for
# the log.
duckdb -readonly "$STAGE/BFD.duckdb" -csv -noheader \
    -c "SELECT count(*) FROM duckdb_views() WHERE NOT internal;" | sed 's/^/views in new db: /'

mv -f "$STAGE/BFD.duckdb" "$RUNDIR/db/BFD.duckdb.new"
mv -f "$RUNDIR/db/BFD.duckdb.new" "$RUNDIR/db/BFD.duckdb"
echo "Published $RUNDIR/db/BFD.duckdb ($(du -h "$RUNDIR/db/BFD.duckdb" | cut -f1)) at $(date -Is)"
