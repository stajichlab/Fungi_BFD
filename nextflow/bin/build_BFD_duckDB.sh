#!/usr/bin/bash -l
#SBATCH -p short --mem 64gb -c 16 -N 1 -n 1 --out logs/build_duckdb.%j.log

# Build the master BFD DuckDB from tables/*.parquet.
#
# Usage:
#   sbatch scripts/build_BFD_duckDB.sh
#   bash   scripts/build_BFD_duckDB.sh
#
# Per T-014 §D.1/§D.2, MERGE_SAMPLES (and every other MERGE_* process) now
# always builds one full, unscoped master table set in tables/ (Parquet, one
# unpartitioned file per type -- see #27) -- there is no more per-taxon
# tables/<Taxon>/ subset or db/BFD.<SUBSET>.duckdb. A taxon-scoped view is a
# post-hoc extraction from this master DB (EXTRACT_BFD_TAXONOMIC_SUBSET,
# tracked separately), not a build-time option here.
#
# Full rebuild every run (rm -f "$DB"), no incremental upsert, per §D.1.
# Since 2026-11-07 the DB holds VIEWs over tables/*.parquet (absolute paths),
# not copies; only swissprot_subject (575k rows) is materialized.

set -euo pipefail

module load duckdb

TABLES="${TABLES:-tables}"
DBDIR="${DBDIR:-db}"
DBNAME=BFD

SRC="$TABLES"

if [ ! -d "$SRC" ]; then
    echo "ERROR: tables folder '$SRC' does not exist." >&2
    exit 1
fi

mkdir -p "$DBDIR"
mkdir -p logs

DB="$DBDIR/$DBNAME.duckdb"
# Remove existing DB to allow clean rebuild
rm -f "$DB"

echo "Building $DB from $SRC"

# ─── species / asm_stats / busco_genome: schema-contract-driven VIEWs ───────
# These three are VIEWs over tables/*.parquet, not materialized copies -- see
# sql/table_schema.json and
# docs/superpowers/specs/2026-09-09-bfd-duckdb-datalake-design.md. The other
# types below use the same VIEW pattern (mkview).
BIN_DIR="${BIN_DIR:?BIN_DIR must be set to an absolute path to nextflow/bin}"
SCHEMA="${SCHEMA:?SCHEMA must be set to an absolute path to sql/table_schema.json}"
# A VIEW's read_parquet('...') path is evaluated at *query* time by whatever
# process opens db/BFD.duckdb (e.g. the MCP server, possibly from a
# different cwd), so it must be absolute -- resolve $SRC defensively even
# though the production caller (BUILD_DUCKDB/main.nf) already passes
# params.tables, which is itself already absolute.
SRC_ABS="$(cd "$SRC" && pwd)"
CATALOG_SQL="$(mktemp)"
python3 "$BIN_DIR/generate_bfd_catalog_sql.py" --schema "$SCHEMA" --tables-dir "$SRC_ABS" > "$CATALOG_SQL"
duckdb "$DB" < "$CATALOG_SQL"
rm -f "$CATALOG_SQL"

# ─── per-type VIEWs over tables/*.parquet ─────────────────────────────────────
# Changed 2026-11-07: every table below is a VIEW, not a CREATE TABLE copy.
# The copies made db/BFD.duckdb 394 GB for ~127 GB of Parquet (2026-11-06
# rebuild, job 31810018). A VIEW reads the Parquet at query time, so:
#   - the paths are absolute ($SRC_ABS), and db/BFD.duckdb works only where
#     tables/ is at the same path;
#   - there are no ART indexes. Measured 2026-11-07 (16 threads, same
#     queries): one genome in gene_proteins 8.5 s (views) vs 1.9 s (copies);
#     one protein in pfam 15.9 vs 13.3 s; v_protein_annotation for one genome
#     68.8 vs 54.5 s; v_species_summary 34.9 vs 30.0 s. Row-group pruning
#     needs rows sorted by id: pfam.sql now writes ORDER BY protein_id (the
#     2026-11-06 pfam.parquet was not sorted: 2,987 of 2,988 row groups overlap);
#   - the UNIQUE indexes were also the duplicate-key check. That check is now
#     an explicit query (check_unique), on by default; CHECK_UNIQUE=0 skips it.
CHECK_UNIQUE="${CHECK_UNIQUE:-1}"

# mkview NAME FILE REQUIRED [EXTRA_SELECT]
#   REQUIRED=1: a missing FILE stops the build; 0: the view is skipped.
mkview() {
    local name=$1 file=$2 required=$3 extra=${4:-}
    if [ ! -f "$SRC/$file" ]; then
        if [ "$required" = 1 ]; then
            echo "ERROR: $SRC/$file not found (required for $name)." >&2; exit 1
        fi
        echo "SKIP: $SRC/$file not found; $name view omitted."
        return 0
    fi
    duckdb -c "CREATE VIEW $name AS SELECT *${extra:+, $extra} FROM read_parquet('$SRC_ABS/$file');" "$DB"
}

# check_unique NAME KEYCOLS -- fail the build when KEYCOLS repeat in NAME.
check_unique() {
    local name=$1 keys=$2
    [ "$CHECK_UNIQUE" = 1 ] || return 0
    duckdb -c "SELECT CASE WHEN count(*) > 0 THEN error('duplicate ($keys) in $name: ' || count(*) || ' keys') END
               FROM (SELECT 1 FROM $name GROUP BY $keys HAVING count(*) > 1);" "$DB" >/dev/null \
        || { echo "ERROR: duplicate keys ($keys) in $name" >&2; exit 1; }
}
has_view() {
    [ "$(duckdb -csv -noheader -c "SELECT count(*) FROM duckdb_views() WHERE view_name = '$1'" "$DB")" = 1 ]
}

# ─── telomeres ────────────────────────────────────────────────────────────────
# Produced by the dedicated telomere-finder pipeline (MERGE_TELOMERES), separate
# from AAFTF's simple motif-count heuristic already embedded in asm_stats as
# t2t_scaffolds/telomere_fwd/telomere_rev -- the two can and will disagree; this
# is the detailed, per-tract-capable source. Both tables carry ASMID and
# LOCUSTAG directly (looked up from samples.csv at merge time).
mkview telomere_summary telomere_summary.parquet 0
has_view telomere_summary && check_unique telomere_summary ASMID
mkview telomere_tracts telomere_tracts.parquet 0
if has_view telomere_tracts; then
    check_unique telomere_tracts tract_id
    # Scaffolds carrying a *terminal* tract on both ends -- the chromosome-
    # completeness signal. terminal=true excludes interstitial telomeric
    # sequences (ITS). Column is end_type, not end (END is reserved SQL).
    duckdb -c "
CREATE VIEW v_telomere_both_ends AS
SELECT ASMID, scaffold
FROM telomere_tracts
WHERE terminal
GROUP BY ASMID, scaffold
HAVING count(DISTINCT end_type) = 2;
" "$DB"
fi

# ─── gene tables ──────────────────────────────────────────────────────────────
# Sequence text is not here: proteins and CDS are in seqdb/ (bgzip + faidx).
mkview gene_info        gene_info.parquet        1
mkview gene_proteins    gene_proteins.parquet    1 "string_split(protein_id,'_')[1] AS species_prefix"
mkview gene_transcripts gene_transcripts.parquet 1 "string_split(transcript_id,'_')[1] AS species_prefix"
mkview gene_exons       gene_exons.parquet       1 "string_split(exon_id,'_')[1] AS species_prefix"
mkview gene_CDS         gene_CDS.parquet         1 "string_split(cds_id,'_')[1] AS species_prefix"
mkview gene_introns     gene_introns.parquet     1 "string_split(transcript_id,'_')[1] AS species_prefix"
mkview gene_trna        gene_trnas.parquet       1 "string_split(gene_id,'_')[1] AS species_prefix"
mkview gene_intergenic_distances gene_intergenic_distances.parquet 1
check_unique gene_info        gene_id
check_unique gene_proteins    protein_id
check_unique gene_transcripts transcript_id
check_unique gene_exons       exon_id
check_unique gene_CDS         cds_id
check_unique gene_introns     intron_id
check_unique gene_trna        gene_id

# ─── functional annotation ────────────────────────────────────────────────────
mkview pfam          pfam.parquet              1 "string_split(protein_id,'_')[1] AS species_prefix"
mkview cazy          cazy.cazymes_hmm.parquet  1
mkview cazy_overview cazy.overview.parquet     1
mkview merops        merops.parquet            1
mkview tf_inventory  tf_inventory.parquet      0
check_unique pfam          "protein_id, pfam_id, domain_num"
check_unique cazy          "protein_id, HMM_id, s_start"
check_unique cazy_overview protein_id

# ─── subcellular localisation / secretome ─────────────────────────────────────
mkview signalp   signalp.signal_peptide.parquet 1
mkview targetp   targetP.parquet                1
mkview wolfpsort wolfpsort.parquet              1
mkview predgpi   predgpi.parquet                1
mkview tmhmm     tmhmm.parquet                  1

# ─── disorder ─────────────────────────────────────────────────────────────────
mkview idp_summary idp_summary.parquet 1
mkview idp         idp.parquet         1
check_unique idp_summary protein_id
check_unique idp         "protein_id, IDP_start"

# ─── composition ──────────────────────────────────────────────────────────────
mkview aa_frequency    aa_freq.parquet    1
mkview codon_frequency codon_freq.parquet 1

# ─── paralogoscope / WGD ───────────────────────────────────────────────────────
# Optional: only exists after the paralogoscope workflow's terminal merge
# (MERGE_WGD_KSD -> tables/wgd.ks.parquet). genome is the sampletag; join on
# species_prefix (the LOCUSTAG) for BFD cross-tables. Rows with NULL dS are
# non-dup paralog pairs without a valid Ka/Ks estimate, not errors.
mkview wgd_ks         wgd.ks.parquet         0
mkview wgd_ks_summary wgd.ks.summary.parquet 0
has_view wgd_ks_summary && check_unique wgd_ks_summary genome

# ─── SwissProt homology ───────────────────────────────────────────────────────
# Design: Fungi_BFD_runs/analysis/swissprot_design/REVIEW.md (2026-11-06).
# swissprot_subject (575k rows) is materialized; the hit tables are VIEWs over
# Parquet (2.2 B rows in swissprot_hits; a copy inside the DB would duplicate
# ~40 GB). Paths are absolute ($SRC_ABS) because VIEWs resolve at query time.
if [ -f "$SRC/swissprot_subject.parquet" ]; then
duckdb -c "
CREATE TABLE swissprot_subject AS SELECT * FROM read_parquet('$SRC_ABS/swissprot_subject.parquet');
CREATE UNIQUE INDEX idx_sp_subject_id  ON swissprot_subject(sp_id);
CREATE UNIQUE INDEX idx_sp_subject_acc ON swissprot_subject(accession);
" "$DB"
else
    echo "SKIP: $SRC/swissprot_subject.parquet not found; SwissProt tables omitted."
fi
if [ -f "$SRC/swissprot_subject.parquet" ] && ls "$SRC"/swissprot_hits/swissprot_hits.*.parquet >/dev/null 2>&1; then
duckdb -c "
CREATE VIEW swissprot_hits AS
    SELECT * FROM read_parquet('$SRC_ABS/swissprot_hits/swissprot_hits.*.parquet');
CREATE VIEW v_swissprot_80_80 AS
    SELECT * FROM swissprot_hits
    WHERE pident_x10 >= 800 AND qcov_x10 >= 800 AND scov_x10 >= 800;
-- For one protein or genome use swissprot_model_fungi_for(id) (below): this
-- VIEW opens all 256 bucket files (about 10 s per protein).
CREATE VIEW v_swissprot_model_fungi AS
    SELECT h.*, s.org_code, s.accession, s.entry_name, s.protein_name
    FROM swissprot_hits h JOIN swissprot_subject s USING (sp_id)
    WHERE s.org_code IN ('YEAST','SCHPO','NEUCR','EMENI')
    QUALIFY row_number() OVER (PARTITION BY h.protein_id, s.org_code ORDER BY h.rk) = 1;
" "$DB"
fi
# Fast lookups by protein or genome (table macros). swissprot_hits is 256
# bucket files; a filter on protein_id through the VIEW still opens all 256
# footers (7-10 s per protein, 2026-11-08). The macros compute the one bucket
# file from the id instead: bucket = first 2 hex of sha1(LOCUSTAG), the same
# as hashBucketForType('swissprot', locustag) in modules/common/utils.nf
# (bucketWidthFor('swissprot') = 2). Measured 0.04 s for one protein and
# 0.16 s for its model-fungi best hits.
#   SELECT * FROM swissprot_hits_for('F77FFF1D_000002-T1');   -- one protein
#   SELECT * FROM swissprot_hits_for('F77FFF1D');             -- one genome
#   SELECT * FROM swissprot_model_fungi_for('F77FFF1D_000002-T1');
if [ -f "$SRC/swissprot_subject.parquet" ] && ls "$SRC"/swissprot_hits/swissprot_hits.*.parquet >/dev/null 2>&1; then
duckdb -c "
CREATE MACRO swissprot_hits_file(id) AS
    '$SRC_ABS/swissprot_hits/swissprot_hits.' || left(sha1(split_part(id, '_', 1)), 2) || '.parquet';
CREATE MACRO swissprot_hits_for(id) AS TABLE
    SELECT * FROM read_parquet(swissprot_hits_file(id))
    WHERE protein_id = id OR species_prefix = id;
CREATE MACRO swissprot_model_fungi_for(id) AS TABLE
    SELECT h.*, s.org_code, s.accession, s.entry_name, s.protein_name
    FROM swissprot_hits_for(id) h JOIN swissprot_subject s USING (sp_id)
    WHERE s.org_code IN ('YEAST','SCHPO','NEUCR','EMENI')
    QUALIFY row_number() OVER (PARTITION BY h.protein_id, s.org_code ORDER BY h.rk) = 1;
" "$DB"
fi
for t in swissprot_best swissprot_subject_summary; do
    if [ -f "$SRC/$t.parquet" ]; then
        duckdb -c "CREATE VIEW $t AS SELECT * FROM read_parquet('$SRC_ABS/$t.parquet');" "$DB"
    else
        echo "SKIP: $SRC/$t.parquet not found; $t view omitted."
    fi
done

# ─── analytical views ─────────────────────────────────────────────────────────
# Per-species summary: genome stats + functional domain counts
duckdb -c "
CREATE VIEW v_species_summary AS
SELECT
    sp.LOCUSTAG,
    sp.SPECIES,
    sp.GENUS,
    sp.FAMILY,
    sp.\"ORDER\",
    sp.CLASS,
    sp.SUBPHYLUM,
    sp.PHYLUM,
    sp.BUSCO_LINEAGE,
    asm.TOTAL_LENGTH,
    asm.GC_PERCENT,
    asm.contig_count,
    asm.N50_bp,
    asm.masked_pct,
    prot.gene_count,
    prot.mean_protein_length,
    COALESCE(pfam_c.pfam_proteins,    0) AS pfam_proteins,
    COALESCE(cazy_c.cazy_proteins,    0) AS cazy_proteins,
    COALESCE(merops_c.merops_proteins,0) AS merops_proteins,
    COALESCE(signalp_c.signalp_count, 0) AS signalp_count,
    COALESCE(tmhmm_c.tmhmm_count,     0) AS tmhmm_count,
    COALESCE(idp_c.idp_high_disorder,  0) AS idp_high_disorder
FROM species sp
JOIN asm_stats asm USING(LOCUSTAG)
JOIN (SELECT species_prefix AS LOCUSTAG,
             count(*)       AS gene_count,
             avg(length)    AS mean_protein_length
      FROM gene_proteins GROUP BY species_prefix) prot USING(LOCUSTAG)
LEFT JOIN (SELECT species_prefix AS LOCUSTAG, count(DISTINCT protein_id) AS pfam_proteins
           FROM pfam GROUP BY species_prefix) pfam_c USING(LOCUSTAG)
LEFT JOIN (SELECT species_prefix AS LOCUSTAG, count(DISTINCT protein_id) AS cazy_proteins
           FROM cazy GROUP BY species_prefix) cazy_c USING(LOCUSTAG)
LEFT JOIN (SELECT species_prefix AS LOCUSTAG, count(DISTINCT protein_id) AS merops_proteins
           FROM merops GROUP BY species_prefix) merops_c USING(LOCUSTAG)
LEFT JOIN (SELECT species_prefix AS LOCUSTAG, count(*) AS signalp_count
           FROM signalp GROUP BY species_prefix) signalp_c USING(LOCUSTAG)
LEFT JOIN (SELECT species_prefix AS LOCUSTAG, count(*) AS tmhmm_count
           FROM tmhmm WHERE PredHel > 0 GROUP BY species_prefix) tmhmm_c USING(LOCUSTAG)
LEFT JOIN (SELECT species_prefix AS LOCUSTAG, count(*) AS idp_high_disorder
           FROM idp_summary WHERE IDP_fraction > 0.8 GROUP BY species_prefix) idp_c USING(LOCUSTAG);
" "$DB"

# Per-protein secretome/annotation flags
duckdb -c "
CREATE VIEW v_protein_annotation AS
SELECT
    gp.protein_id,
    gp.species_prefix  AS LOCUSTAG,
    gp.length,
    COALESCE(pfam_c.pfam_count, 0)    AS pfam_domain_count,
    COALESCE(cazy_p.cazyme, FALSE)     AS has_cazyme,
    COALESCE(merops_p.peptidase, FALSE)AS has_merops,
    COALESCE(sp.has_signal, FALSE)     AS has_signal_peptide,
    COALESCE(tp.prediction, 'OTHER')   AS targetp_class,
    wp.localization                    AS wolfpsort_loc,
    COALESCE(tm.PredHel, 0)            AS tmhmm_helices,
    COALESCE(gpi.is_gpi, FALSE)        AS is_gpi_anchored,
    idps.IDP_fraction
FROM gene_proteins gp
LEFT JOIN (SELECT protein_id, count(*) AS pfam_count FROM pfam GROUP BY protein_id) pfam_c USING(protein_id)
LEFT JOIN (SELECT DISTINCT protein_id, TRUE AS cazyme FROM cazy_overview) cazy_p USING(protein_id)
LEFT JOIN (SELECT DISTINCT protein_id, TRUE AS peptidase FROM merops) merops_p USING(protein_id)
LEFT JOIN (SELECT DISTINCT protein_id, TRUE AS has_signal FROM signalp) sp USING(protein_id)
LEFT JOIN (SELECT protein_id, prediction FROM targetp) tp USING(protein_id)
LEFT JOIN (SELECT protein_id, localization FROM wolfpsort) wp USING(protein_id)
LEFT JOIN (SELECT protein_id, PredHel FROM tmhmm) tm USING(protein_id)
LEFT JOIN (SELECT DISTINCT protein_id, TRUE AS is_gpi FROM predgpi) gpi USING(protein_id)
LEFT JOIN (SELECT protein_id, IDP_fraction FROM idp_summary) idps USING(protein_id);
" "$DB"

echo "Build complete: $DB ($(du -h "$DB" | cut -f1))"
duckdb -c "SELECT table_name, estimated_size FROM duckdb_tables() ORDER BY table_name;" "$DB"
duckdb -c "SELECT view_name FROM duckdb_views() WHERE NOT internal ORDER BY view_name;" "$DB"
