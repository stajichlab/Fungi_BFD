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
# Full rebuild every run (rm -f "$DB" + CREATE TABLE AS SELECT), same as
# before Parquet -- no incremental upsert, per §D.1's decision.

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

# ─── species ──────────────────────────────────────────────────────────────────
duckdb -c "
CREATE TABLE species AS SELECT * FROM read_parquet('$SRC/species.parquet');
CREATE UNIQUE INDEX idx_species_locustag ON species(LOCUSTAG);
CREATE UNIQUE INDEX idx_species_asm      ON species(ASMID);
-- Non-unique: a genus-level taxonomy tuple is shared by multiple species.
CREATE INDEX idx_species_taxonomy ON species(PHYLUM, CLASS, \"ORDER\", FAMILY, GENUS);
" "$DB"

# ─── asm_stats ────────────────────────────────────────────────────────────────
# Source has ASMID/gc_pct/total_length_bp; we join to species for LOCUSTAG
# and alias columns to match existing R-script conventions (GC_PERCENT, TOTAL_LENGTH).
duckdb -c "
CREATE TABLE asm_stats AS
SELECT
    sp.LOCUSTAG,
    a.ASMID,
    a.contig_count,
    a.total_length_bp   AS TOTAL_LENGTH,
    a.min_contig_bp,
    a.max_contig_bp,
    a.median_contig_bp,
    a.mean_contig_bp,
    a.L50, a.N50_bp, a.L90, a.N90_bp,
    a.gc_pct            AS GC_PERCENT,
    a.n_gap_count,
    a.total_n_bases,
    a.masked_bases,
    a.masked_pct,
    a.t2t_scaffolds,
    a.telomere_fwd,
    a.telomere_rev
FROM read_parquet('$SRC/asm_stats.parquet') AS a
JOIN species AS sp USING(ASMID);
CREATE UNIQUE INDEX idx_asm_locustag ON asm_stats(LOCUSTAG);
CREATE UNIQUE INDEX idx_asm_asmid    ON asm_stats(ASMID);
" "$DB"

# ─── telomeres ────────────────────────────────────────────────────────────────
# Produced by the dedicated telomere-finder pipeline (MERGE_TELOMERES), separate
# from AAFTF's simple motif-count heuristic already embedded in asm_stats as
# t2t_scaffolds/telomere_fwd/telomere_rev -- the two can and will disagree; this
# is the detailed, per-tract-capable source. Both tables carry ASMID and
# LOCUSTAG directly (looked up from samples.csv at merge time), rather than the
# species_prefix-derived-from-id pattern used by gene_*/functional tables,
# since there is no encoded id string to split here.
if [ -f "$SRC/telomere_summary.parquet" ]; then
    duckdb -c "
CREATE TABLE telomere_summary AS
SELECT * FROM read_parquet('$SRC/telomere_summary.parquet');
CREATE UNIQUE INDEX idx_telosum_asmid    ON telomere_summary(ASMID);
CREATE        INDEX idx_telosum_locustag ON telomere_summary(LOCUSTAG);
" "$DB"
else
    echo "SKIP: $SRC/telomere_summary.parquet not found; telomere_summary table omitted."
fi

if [ -f "$SRC/telomere_tracts.parquet" ]; then
    duckdb -c "
CREATE TABLE telomere_tracts AS
SELECT * FROM read_parquet('$SRC/telomere_tracts.parquet');
CREATE UNIQUE INDEX idx_telotracts_id       ON telomere_tracts(tract_id);
CREATE        INDEX idx_telotracts_asmid    ON telomere_tracts(ASMID);
CREATE        INDEX idx_telotracts_locustag ON telomere_tracts(LOCUSTAG);
CREATE        INDEX idx_telotracts_scaffold ON telomere_tracts(ASMID, scaffold);

-- Scaffolds carrying a *terminal* tract on both ends -- the chromosome-
-- completeness signal. terminal=true excludes interstitial telomeric
-- sequences (ITS), which would otherwise falsely register a scaffold as
-- capped on both ends. Column is end_type, not end -- END is a reserved
-- SQL keyword and 'end' would force every ad hoc query to quote it.
CREATE VIEW v_telomere_both_ends AS
SELECT ASMID, scaffold
FROM telomere_tracts
WHERE terminal
GROUP BY ASMID, scaffold
HAVING count(DISTINCT end_type) = 2;
" "$DB"
else
    echo "SKIP: $SRC/telomere_tracts.parquet not found; telomere_tracts table omitted."
fi

# ─── busco_genome ───────────────────────────────────────────────────────────────
# Assembly-level BUSCO completeness (BUSCO_GENOME), keyed by ASMID -- joined to
# species for LOCUSTAG, same pattern as asm_stats. Not to be confused with
# BUSCO_PEP (annotation-level, protein-set completeness), which is a separate
# downstream QC metric and is not merged here.
# run_busco_genome defaults false (cost at 5K+ species scale), so this table is
# skipped -- not created empty -- when the parquet hasn't been produced yet.
if [ -f "$SRC/busco_genome.parquet" ]; then
    duckdb -c "
CREATE TABLE busco_genome AS
SELECT
    sp.LOCUSTAG,
    b.ASMID,
    b.complete_pct,
    b.single_pct,
    b.duplicated_pct,
    b.fragmented_pct,
    b.missing_pct,
    b.n_markers,
    b.lineage
FROM read_parquet('$SRC/busco_genome.parquet') AS b
JOIN species AS sp USING(ASMID);
CREATE UNIQUE INDEX idx_busco_genome_locustag ON busco_genome(LOCUSTAG);
CREATE UNIQUE INDEX idx_busco_genome_asmid    ON busco_genome(ASMID);
" "$DB"
else
    echo "SKIP: $SRC/busco_genome.parquet not found (run_busco_genome likely false); busco_genome table omitted."
fi

# ─── gene tables ──────────────────────────────────────────────────────────────
duckdb -c "
CREATE TABLE gene_info AS
SELECT * FROM read_parquet('$SRC/gene_info.parquet');
CREATE UNIQUE INDEX idx_gene_info_gene_id  ON gene_info(gene_id);
CREATE        INDEX idx_gene_info_locustag ON gene_info(LOCUSTAG);
" "$DB"

duckdb -c "
CREATE TABLE gene_proteins AS
SELECT *, string_split(protein_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/gene_proteins.parquet');
CREATE UNIQUE INDEX idx_gene_prot_protein_id  ON gene_proteins(protein_id);
CREATE        INDEX idx_gene_prot_locustag    ON gene_proteins(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE gene_transcripts AS
SELECT *, string_split(transcript_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/gene_transcripts.parquet');
CREATE UNIQUE INDEX idx_tx_transcript_id ON gene_transcripts(transcript_id);
CREATE        INDEX idx_tx_locustag      ON gene_transcripts(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE gene_exons AS
SELECT *, string_split(exon_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/gene_exons.parquet');
CREATE UNIQUE INDEX idx_exon_exon_id  ON gene_exons(exon_id);
CREATE        INDEX idx_exon_locustag ON gene_exons(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE gene_CDS AS
SELECT *, string_split(cds_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/gene_CDS.parquet');
CREATE UNIQUE INDEX idx_cds_cds_id   ON gene_CDS(cds_id);
CREATE        INDEX idx_cds_locustag ON gene_CDS(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE gene_introns AS
SELECT *, string_split(transcript_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/gene_introns.parquet');
CREATE UNIQUE INDEX idx_intron_intron_id  ON gene_introns(intron_id);
CREATE        INDEX idx_intron_locustag   ON gene_introns(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE gene_trna AS
SELECT *, string_split(gene_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/gene_trnas.parquet');
CREATE UNIQUE INDEX idx_trna_gene_id  ON gene_trna(gene_id);
CREATE        INDEX idx_trna_locustag ON gene_trna(species_prefix);
" "$DB"

# intergenic distances (large — only species_prefix index)
duckdb -c "
CREATE TABLE gene_intergenic_distances AS
SELECT * FROM read_parquet('$SRC/gene_intergenic_distances.parquet');
CREATE INDEX idx_intergenic_locustag ON gene_intergenic_distances(species_prefix);
" "$DB"

# ─── functional annotation ────────────────────────────────────────────────────
duckdb -c "
CREATE TABLE pfam AS
SELECT *, string_split(protein_id,'_')[1] AS species_prefix
FROM read_parquet('$SRC/pfam.parquet');
CREATE UNIQUE INDEX idx_pfam_domain_hit    ON pfam(protein_id, pfam_id, domain_num);
CREATE        INDEX idx_pfam_locustag      ON pfam(species_prefix);
CREATE        INDEX idx_pfam_pfam_id       ON pfam(pfam_id);
" "$DB"

duckdb -c "
CREATE TABLE cazy AS
SELECT * FROM read_parquet('$SRC/cazy.cazymes_hmm.parquet');
CREATE UNIQUE INDEX idx_cazy_hit        ON cazy(protein_id, HMM_id, s_start);
CREATE        INDEX idx_cazy_locustag   ON cazy(species_prefix);
CREATE        INDEX idx_cazy_hmm_id     ON cazy(HMM_id);
" "$DB"

duckdb -c "
CREATE TABLE cazy_overview AS
SELECT * FROM read_parquet('$SRC/cazy.overview.parquet');
CREATE UNIQUE INDEX idx_cazyov_protein_id ON cazy_overview(protein_id);
CREATE        INDEX idx_cazyov_locustag   ON cazy_overview(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE merops AS
SELECT * FROM read_parquet('$SRC/merops.parquet');
CREATE        INDEX idx_merops_locustag  ON merops(species_prefix);
CREATE        INDEX idx_merops_merops_id ON merops(merops_id);
" "$DB"

# ─── subcellular localisation / secretome ─────────────────────────────────────
duckdb -c "
CREATE TABLE signalp AS
SELECT * FROM read_parquet('$SRC/signalp.signal_peptide.parquet');
CREATE INDEX idx_signalp_locustag  ON signalp(species_prefix);
CREATE INDEX idx_signalp_protein   ON signalp(protein_id);
" "$DB"

duckdb -c "
CREATE TABLE targetp AS
SELECT * FROM read_parquet('$SRC/targetP.parquet');
CREATE INDEX idx_targetp_locustag ON targetp(species_prefix);
CREATE INDEX idx_targetp_protein  ON targetp(protein_id);
" "$DB"

duckdb -c "
CREATE TABLE wolfpsort AS
SELECT * FROM read_parquet('$SRC/wolfpsort.parquet');
CREATE INDEX idx_wolfpsort_locustag ON wolfpsort(species_prefix);
CREATE INDEX idx_wolfpsort_protein  ON wolfpsort(protein_id);
" "$DB"

duckdb -c "
CREATE TABLE predgpi AS
SELECT * FROM read_parquet('$SRC/predgpi.parquet');
CREATE INDEX idx_predgpi_locustag ON predgpi(species_prefix);
CREATE INDEX idx_predgpi_protein  ON predgpi(protein_id);
" "$DB"

duckdb -c "
CREATE TABLE tmhmm AS
SELECT * FROM read_parquet('$SRC/tmhmm.parquet');
CREATE INDEX idx_tmhmm_locustag ON tmhmm(species_prefix);
CREATE INDEX idx_tmhmm_protein  ON tmhmm(protein_id);
" "$DB"

# ─── disorder ─────────────────────────────────────────────────────────────────
duckdb -c "
CREATE TABLE idp_summary AS
SELECT * FROM read_parquet('$SRC/idp_summary.parquet');
CREATE UNIQUE INDEX idx_idpsum_protein   ON idp_summary(protein_id);
CREATE        INDEX idx_idpsum_locustag  ON idp_summary(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE idp AS
SELECT * FROM read_parquet('$SRC/idp.parquet');
CREATE UNIQUE INDEX idx_idp_region    ON idp(protein_id, IDP_start);
CREATE        INDEX idx_idp_locustag  ON idp(species_prefix);
" "$DB"

# ─── composition ──────────────────────────────────────────────────────────────
duckdb -c "
CREATE TABLE aa_frequency AS
SELECT * FROM read_parquet('$SRC/aa_freq.parquet');
CREATE INDEX idx_aa_locustag ON aa_frequency(species_prefix);
" "$DB"

duckdb -c "
CREATE TABLE codon_frequency AS
SELECT * FROM read_parquet('$SRC/codon_freq.parquet');
CREATE INDEX idx_codon_locustag ON codon_frequency(species_prefix);
" "$DB"

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

echo "Build complete: $DB"
duckdb -c "SELECT table_name, estimated_size FROM duckdb_tables() ORDER BY table_name;" "$DB"
duckdb -c "SELECT view_name FROM duckdb_views() ORDER BY view_name;" "$DB"
