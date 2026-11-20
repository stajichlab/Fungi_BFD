#!/usr/bin/bash
#
# Reload ref_class_subphylum and ref_phylum_subphylum in an existing BFD
# DuckDB from lib/refs/*.tsv, without rebuilding the database.
#
#     scripts/refresh_subphylum_refs.sh /path/to/BFD.duckdb
#
# Why this is safe and why it is fast: v_genome_traits is a VIEW that resolves
# SUBPHYLUM at query time as
#
#     coalesce(species.SUBPHYLUM, ref_class_subphylum.SUBPHYLUM,
#                                 ref_phylum_subphylum.SUBPHYLUM)
#
# so replacing the two lookup TABLES updates every view above them -- including
# v_genome_labels -- on the next query. The two tables hold on the order of a
# hundred rows between them. Measured 2026-11-19 against the 445 MB dashboard
# database: 3.7 s, and the per-genome result was identical to the 9 m 51 s full
# rebuild it replaces.
#
# The statements are SLICED OUT of sql/bfd_traits.sql rather than duplicated
# here, so this cannot drift from the canonical SQL. The same REF_DIR= header
# substitution that nextflow/bin/build_BFD_duckDB.sh performs is applied first.
#
# WHAT THIS DOES NOT DO
#   * It does not change view definitions. If the body of v_genome_traits (or
#     anything else) changed in sql/bfd_traits.sql, run the full build instead
#     -- scripts/rebuild_pipeline_duckdb.sh. The guard below refuses to run if
#     the database's view no longer matches the shape this script assumes.
#   * It does not touch species.SUBPHYLUM, the FIRST and dominant arm of that
#     coalesce (24,096 of 24,363 genomes as of 2026-11-19). That column comes
#     from tables/species.parquet via the MERGE_SAMPLES pipeline step. This
#     script only moves genomes that fall through to the lookup tables.
#   * It does not reload ref_fungaltraits or ref_growth_form_override.
#
# DuckDB allows a single writer. Do not run this against a database a build or
# an interactive session is holding; the lock will refuse and nothing changes.
set -euo pipefail

DB="${1:-${BFD_DB:-}}"
[ -n "$DB" ] || { echo "usage: $0 /path/to/BFD.duckdb   (or set BFD_DB)" >&2; exit 2; }
[ -f "$DB" ] || { echo "ERROR: no database at $DB" >&2; exit 1; }
[ -w "$DB" ] || { echo "ERROR: $DB is not writable" >&2; exit 1; }

# Prefer the checkout this script was invoked from, so the refs and the SQL
# match the code under review; fall back to the shared checkout. Same rule as
# scripts/rebuild_pipeline_duckdb.sh.
if [ -z "${PIPELINE:-}" ]; then
    _self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd || true)"
    _candidate="${_self_dir%/scripts}"
    if [ -n "$_candidate" ] && [ -f "$_candidate/sql/bfd_traits.sql" ]; then
        PIPELINE="$_candidate"
    else
        PIPELINE=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD
    fi
fi
TRAITS_SQL="$PIPELINE/sql/bfd_traits.sql"
REF_DIR="${REF_DIR:-$PIPELINE/lib/refs}"
[ -f "$TRAITS_SQL" ] || { echo "ERROR: missing $TRAITS_SQL (set PIPELINE)" >&2; exit 1; }
for f in class_subphylum.tsv phylum_subphylum.tsv; do
    [ -s "$REF_DIR/$f" ] || { echo "ERROR: missing or empty $REF_DIR/$f (set REF_DIR)" >&2; exit 1; }
done
REF_DIR="$(cd "$REF_DIR" && pwd)"
echo "database:  $DB"
echo "checkout:  $PIPELINE"
echo "ref dir:   $REF_DIR"

command -v duckdb >/dev/null 2>&1 || module load duckdb
# DuckDB >= 1.4 turns on an interactive "agent mode" when AI_AGENT is set and
# stdout is not a terminal, which would corrupt the values parsed below.
duck() { env -u AI_AGENT duckdb "$@"; }

# The database must already have the view this script is a shortcut for, and
# that view must still read BOTH lookup tables. If it does not, the SQL has
# moved on from what this script assumes and a full rebuild is the honest
# answer. Checked against the stored (normalised) view SQL, so reformatting
# the source file does not trip it.
# Existence is asked as a count rather than as "is the SQL text empty":
# duckdb's CSV writer renders an empty string as a quoted "", which is not
# empty to the shell, so an absent view would otherwise fall through to the
# shape check below and be reported as the wrong problem.
HAVE_VIEW="$(duck -readonly -csv -noheader "$DB" \
    -c "SELECT count(*) FROM duckdb_views() WHERE view_name = 'v_genome_traits';")"
[ "$HAVE_VIEW" = "1" ] || {
    echo "ERROR: $DB has no v_genome_traits view; build it first" >&2; exit 1; }
VIEW_SQL="$(duck -readonly -csv -noheader "$DB" \
    -c "SELECT sql FROM duckdb_views() WHERE view_name = 'v_genome_traits';")"
for t in ref_class_subphylum ref_phylum_subphylum; do
    case "$VIEW_SQL" in
        *"$t"*) ;;
        *) echo "ERROR: v_genome_traits no longer reads $t -- the SQL has changed shape." >&2
           echo "       Run a full build (scripts/rebuild_pipeline_duckdb.sh) instead." >&2
           exit 1 ;;
    esac
done

counts() {
    duck -readonly -csv -noheader "$DB" -c \
        "SELECT (SELECT count(*) FROM ref_class_subphylum)
             || ' ' || (SELECT count(*) FROM ref_phylum_subphylum)
             || ' ' || (SELECT count(*) FROM v_genome_traits WHERE SUBPHYLUM IS NULL);"
}
read -r OLD_CLASS OLD_PHYLUM OLD_NULL <<<"$(counts)"
echo "before:    ref_class=$OLD_CLASS ref_phylum=$OLD_PHYLUM subphylum_null=$OLD_NULL"

# Substitute the REF_DIR= header exactly as build_BFD_duckDB.sh does, then keep
# only the two subphylum statements (and their COMMENTs).
DEFAULT_REF="$(sed -n '1s/^-- REF_DIR=//p' "$TRAITS_SQL")"
[ -n "$DEFAULT_REF" ] || { echo "ERROR: $TRAITS_SQL has no '-- REF_DIR=' header on line 1" >&2; exit 1; }
case "$REF_DIR$DEFAULT_REF" in
    *[\#\&\\]*) echo "ERROR: ref paths must not contain '#', '&' or a backslash" >&2; exit 1 ;;
esac
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT
sed "s#${DEFAULT_REF}#${REF_DIR}#g" "$TRAITS_SQL" \
  | sed -n '/^CREATE OR REPLACE TABLE ref_class_subphylum/,/^COMMENT ON TABLE ref_phylum_subphylum/p' \
  > "$STAGE/subphylum_only.sql"

# The slice must be exactly the two tables. Too few means bfd_traits.sql was
# reordered or renamed and the slice silently grabbed nothing; too many means
# it swallowed an unrelated statement.
N_CREATE=$(grep -c '^CREATE OR REPLACE TABLE ' "$STAGE/subphylum_only.sql" || true)
[ "$N_CREATE" -eq 2 ] || {
    echo "ERROR: expected 2 CREATE statements in the slice, found $N_CREATE." >&2
    echo "       sql/bfd_traits.sql has been reordered; fix the slice in this script." >&2
    exit 1; }
grep -qE "read_csv\('${REF_DIR}/class_subphylum\.tsv'"  "$STAGE/subphylum_only.sql" \
  && grep -qE "read_csv\('${REF_DIR}/phylum_subphylum\.tsv'" "$STAGE/subphylum_only.sql" \
  || { echo "ERROR: the slice does not read both TSVs from $REF_DIR" >&2; exit 1; }

# One transaction, so a failure part-way leaves both tables as they were.
{ echo "BEGIN TRANSACTION;"; cat "$STAGE/subphylum_only.sql"; echo "COMMIT;"; } \
    > "$STAGE/apply.sql"
duck "$DB" < "$STAGE/apply.sql"

read -r NEW_CLASS NEW_PHYLUM NEW_NULL <<<"$(counts)"
echo "after:     ref_class=$NEW_CLASS ref_phylum=$NEW_PHYLUM subphylum_null=$NEW_NULL"

# Losing lookup rows is almost always an edit mistake rather than an intent,
# so say so loudly. ALLOW_SHRINK=1 to proceed anyway, matching the
# --allow-shrink idiom in nextflow/bin/publish_table.py.
if [ "${ALLOW_SHRINK:-0}" != "1" ] \
   && { [ "$NEW_CLASS" -lt "$OLD_CLASS" ] || [ "$NEW_PHYLUM" -lt "$OLD_PHYLUM" ]; }; then
    echo "ERROR: lookup rows dropped (class $OLD_CLASS->$NEW_CLASS, phylum $OLD_PHYLUM->$NEW_PHYLUM)." >&2
    echo "       The tables have ALREADY been replaced; re-run with the intended" >&2
    echo "       $REF_DIR TSVs, or set ALLOW_SHRINK=1 if the removal is deliberate." >&2
    exit 1
fi
echo "subphylum refs refreshed: $((OLD_NULL - NEW_NULL)) genome(s) gained a SUBPHYLUM"
