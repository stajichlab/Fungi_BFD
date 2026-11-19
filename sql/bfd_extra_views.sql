-- TABLES_DIR=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
-- extra_views.sql: helper views for the local BFD.duckdb. Each view removes a join or a trap that
-- repeats in the BFD_dashboard analyses (evidence: REPORT.md section 4).
-- Run after build_BFD_duckDB.sh AND after compat_views.sql (v_genome_labels reads the `samples` view):
--     duckdb BFD.duckdb < compat_views.sql
--     duckdb BFD.duckdb < bfd_traits.sql     (v_genome_labels reads v_genome_traits)
--     duckdb BFD.duckdb < extra_views.sql
-- No statement in this file reads a parquet path directly, so the TABLES_DIR line is informative only.
-- All statements are idempotent. Nothing is materialized. No view reads swissprot_hits.
-- Tested 2026-11-11 with DuckDB CLI 1.1.3 and python duckdb 1.5.5 (timings in REPORT.md section 6).
--
-- Genome key: per-genome views (v_genome_*) use LOCUSTAG, as species and v_species_summary do.
-- Views with one row per protein or per family keep species_prefix, as their source tables do.

-- ---------------------------------------------------------------------------
-- in_genome(id, locustag): genome filter that prunes row groups
-- ---------------------------------------------------------------------------
-- Pattern: left(id,8) IN (...), split_part(id,'_',1) = ..., LIKE 'X%' and species_prefix = 'X' on pfam and
-- gene_* (docs/analyses/kickxellomycota/README.md:447-457, kingdom_domains/README.md:48). On pfam,
-- gene_proteins, gene_transcripts, gene_exons, gene_CDS, gene_introns and gene_trna the species_prefix column
-- is computed by the view, so a filter on it reads every row group. This macro expands to a range on the
-- sorted id column, which DuckDB pushes into the parquet scan.
--   SELECT * FROM pfam WHERE in_genome(protein_id, 'F07D2B4E');
--   SELECT * FROM gene_introns WHERE in_genome(transcript_id, 'F07D2B4E');
-- '`' is the ASCII character after '_', so the range holds exactly the ids that start with '<LOCUSTAG>_'.
CREATE OR REPLACE MACRO in_genome(id, locustag) AS (id >= locustag || '_' AND id < locustag || '`');

-- ---------------------------------------------------------------------------
-- v_genome_labels: one row per genome with taxonomy, assembly, BUSCO, quality flags, representatives
-- ---------------------------------------------------------------------------
-- Pattern: every analysis first joins species + asm_stats + busco_genome and then picks one genome per
-- species (kingdom_domains/README.md:37-42 kw_rep_genomes; kickxellomycota/README.md:447-448 kx_genomes;
-- bin/atlas/sql/J3_protein_core.sql:35-36; lib/subset.py:106-113). kickxellomycota/README.md:55 states that
-- representatives and quality flags are not in the database.
-- Rules are those of BFD_dashboard/bin/atlas/build_labels.py:
--   is_representative = best busco_complete per SPECIES, then N50_bp, then LOCUSTAG; annotated genomes only
--   is_genus_rep      = best species representative per GENUS (same order)
--   q_ann = annotated; q_dup = busco_duplicated <= 10; q_n50 = N50_bp >= 50000
--   q_busco_rel = busco_complete >= 25th percentile of its CLASS (classes with >= 20 annotated genomes), else of its PHYLUM
-- Column names follow atlas_data/labels/genome_labels.parquet (ORDER_, total_length_bp, gc_pct, busco_*).
-- FungalTraits columns (growth_form, primary_lifestyle, ...) and SUBPHYLUM come from v_genome_traits, which
-- sql/bfd_traits.sql defines; apply that file BEFORE this one.
CREATE OR REPLACE VIEW v_genome_labels AS
WITH base AS (
  SELECT s.LOCUSTAG, s.ASMID, s.PHYLUM, s.CLASS, s."ORDER" AS ORDER_, s.FAMILY, s.GENUS, s.SPECIES, s.STRAIN,
         s.NCBI_TAXONID, sa.TRANSL_TABLE AS transl_table,
         a.annotation_status, a.TOTAL_LENGTH AS total_length_bp, a.contig_count, a.N50_bp,
         a.GC_PERCENT AS gc_pct, a.masked_pct,
         b.complete_pct AS busco_complete, b.single_pct AS busco_single, b.duplicated_pct AS busco_duplicated,
         b.fragmented_pct AS busco_fragmented, b.missing_pct AS busco_missing,
         t.* EXCLUDE (LOCUSTAG)
  FROM species s
  LEFT JOIN v_genome_traits t ON s.LOCUSTAG = t.LOCUSTAG
  LEFT JOIN samples sa ON s.LOCUSTAG = sa.LOCUSTAG
  LEFT JOIN asm_stats a ON s.ASMID = a.ASMID
  LEFT JOIN busco_genome b ON s.ASMID = b.ASMID
), cls_thr AS (
  SELECT CLASS, quantile_cont(busco_complete, 0.25) AS p25, quantile_cont(busco_complete, 0.5) AS p50
  FROM base WHERE annotation_status = 'annotated' AND CLASS IS NOT NULL
  GROUP BY CLASS HAVING count(*) >= 20
), phy_thr AS (
  SELECT PHYLUM, quantile_cont(busco_complete, 0.25) AS p25, quantile_cont(busco_complete, 0.5) AS p50
  FROM base WHERE annotation_status = 'annotated' AND PHYLUM IS NOT NULL
  GROUP BY PHYLUM
), flagged AS (
  SELECT l.*,
         coalesce(c.p25, p.p25) AS busco_thr_p25,
         coalesce(c.p50, p.p50) AS busco_thr_p50,
         CASE WHEN c.CLASS IS NOT NULL THEN 'class' WHEN p.PHYLUM IS NOT NULL THEN 'phylum' END AS busco_thr_level,
         (l.annotation_status = 'annotated') AS q_ann,
         (l.busco_complete >= coalesce(c.p25, p.p25)) AS q_busco_rel,
         (l.busco_complete >= coalesce(c.p50, p.p50)) AS q_busco_rel50,
         (l.busco_duplicated <= 10) AS q_dup,
         (l.N50_bp >= 50000) AS q_n50,
         (row_number() OVER (PARTITION BY l.SPECIES
                             ORDER BY l.busco_complete DESC NULLS LAST, l.N50_bp DESC NULLS LAST, l.LOCUSTAG) = 1
          AND l.annotation_status = 'annotated' AND l.SPECIES IS NOT NULL) AS is_representative
  FROM base l
  LEFT JOIN cls_thr c ON l.CLASS = c.CLASS
  LEFT JOIN phy_thr p ON l.PHYLUM = p.PHYLUM
)
SELECT f.*,
       (f.is_representative AND f.GENUS IS NOT NULL AND
        row_number() OVER (PARTITION BY f.GENUS, f.is_representative
                           ORDER BY f.busco_complete DESC NULLS LAST, f.N50_bp DESC NULLS LAST, f.LOCUSTAG) = 1) AS is_genus_rep
FROM flagged f;
COMMENT ON VIEW v_genome_labels IS 'One row per genome: taxonomy (ORDER_ has a trailing underscore), assembly and BUSCO values, quality flags q_ann, q_busco_rel, q_busco_rel50, q_dup, q_n50, and is_representative (one genome per SPECIES) and is_genus_rep (one per GENUS). Same rules and column names as atlas genome_labels.parquet (bin/atlas/build_labels.py); FungalTraits columns (growth_form, growth_form_source, primary_lifestyle, secondary_lifestyle, plant_pathogenic_capacity, decay_substrate, animal_biotrophic_capacity, in_fungaltraits, yeast_origin) are a GENUS-level label from FungalTraits 1.2 (in_fungaltraits false = genus not listed); SUBPHYLUM is NULL when the class has no entry in ref_class_subphylum (subphylum_source says where it came from). Use WHERE is_representative before counting genomes per clade, or strains are counted as species. Key LOCUSTAG.';

-- ---------------------------------------------------------------------------
-- v_genome_coverage: which result tables have rows for each genome
-- ---------------------------------------------------------------------------
-- Pattern/trap: tmhmm, signalp, targetp, cazy_overview, merops and tf_inventory list only positive proteins,
-- and some tools ran on a subset (targetp: 125 genomes). "No rows" is then read as "zero"
-- (kickxellomycota/README.md:59; docs/FUNCTIONAL_ANALYSES.md:22-32; lib/bfd_data.py:61-75 counts genomes
-- per table for the same reason). has_* = the table has at least one row for the genome.
-- pfam, gene_* and idp are not included: their genome key needs a scan of the whole table.
-- Cost: one pass over the species_prefix column of each table (13-26 s with DuckDB 1.1.3 and 4 threads).
CREATE OR REPLACE VIEW v_genome_coverage AS
SELECT s.LOCUSTAG, s.ASMID, a.annotation_status,
       (b.ASMID IS NOT NULL)              AS has_busco,
       (pr.species_prefix IS NOT NULL)    AS has_proteins,
       (sg.species_prefix IS NOT NULL)    AS has_signalp,
       (tm.species_prefix IS NOT NULL)    AS has_tmhmm,
       (gp.species_prefix IS NOT NULL)    AS has_predgpi,
       (wp.species_prefix IS NOT NULL)    AS has_wolfpsort,
       (tp.species_prefix IS NOT NULL)    AS has_targetp,
       (cz.species_prefix IS NOT NULL)    AS has_cazy,
       (me.species_prefix IS NOT NULL)    AS has_merops,
       (tf.species_prefix IS NOT NULL)    AS has_tf,
       (ids.species_prefix IS NOT NULL)   AS has_idp_summary,
       (sb.species_prefix IS NOT NULL)    AS has_swissprot_best,
       (wk.species_prefix IS NOT NULL)    AS has_wgd_ks,
       (te.LOCUSTAG IS NOT NULL)          AS has_telomere_summary
FROM species s
LEFT JOIN asm_stats a ON a.ASMID = s.ASMID
LEFT JOIN busco_genome b ON b.ASMID = s.ASMID
LEFT JOIN (SELECT DISTINCT species_prefix FROM aa_frequency)   pr  ON pr.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM signalp)        sg  ON sg.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM tmhmm)          tm  ON tm.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM predgpi)        gp  ON gp.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM wolfpsort)      wp  ON wp.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM targetp)        tp  ON tp.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM cazy_overview)  cz  ON cz.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM merops)         me  ON me.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM tf_inventory)   tf  ON tf.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM idp_summary)    ids ON ids.species_prefix = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM swissprot_best) sb  ON sb.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT species_prefix FROM wgd_ks_summary) wk  ON wk.species_prefix  = s.LOCUSTAG
LEFT JOIN (SELECT DISTINCT LOCUSTAG FROM telomere_summary)     te  ON te.LOCUSTAG        = s.LOCUSTAG;
COMMENT ON VIEW v_genome_coverage IS 'One row per genome: has_* is true when the named table has at least one row for the genome. Check it before reading a missing row as zero: tmhmm, signalp, targetp, cazy_overview, merops and tf_inventory list only positive proteins, and targetp covers 125 genomes. has_proteins comes from aa_frequency. pfam, gene_* and idp are not covered. Key LOCUSTAG.';

-- ---------------------------------------------------------------------------
-- v_secretome_proteins / v_genome_secretome: SignalP + TMHMM + PredGPI secretome
-- ---------------------------------------------------------------------------
-- Pattern: the same three-table join is written twice in bin/atlas/sql/J4_function_counts.sql (lines 5-34
-- and 39-46) and again as separate counts in lib/bfd_data.py:445-458. Trap: predgpi has a 'Chain' row for
-- nearly every protein; only feature = 'GPI-anchor' marks an anchor. Definitions are those of J4:
--   secreted_noTM   = SignalP probability >= 0.9 and at most 1 TMHMM helix (the signal peptide can count as one)
--   secreted_strict = secreted_noTM and no GPI anchor
-- One row per protein with a SignalP signal peptide. Proteins without a signal peptide are not listed.
-- Filter by species_prefix (all three sources are sorted by it), not by a LIKE on protein_id.
CREATE OR REPLACE VIEW v_secretome_proteins AS
SELECT s.species_prefix, s.protein_id,
       s.probability AS signalp_prob, s.peptide_end AS signalp_end,
       coalesce(t.PredHel, 0) AS tm_helices,
       (g.protein_id IS NOT NULL) AS is_gpi,
       (s.probability >= 0.9 AND coalesce(t.PredHel, 0) <= 1) AS secreted_noTM,
       (s.probability >= 0.9 AND coalesce(t.PredHel, 0) <= 1 AND g.protein_id IS NULL) AS secreted_strict
FROM signalp s
LEFT JOIN tmhmm t ON t.species_prefix = s.species_prefix AND t.protein_id = s.protein_id
LEFT JOIN (SELECT DISTINCT species_prefix, protein_id FROM predgpi WHERE feature = 'GPI-anchor') g
       ON g.species_prefix = s.species_prefix AND g.protein_id = s.protein_id;
COMMENT ON VIEW v_secretome_proteins IS 'One row per protein with a SignalP signal peptide: signalp_prob, signalp_end (last residue of the peptide), tm_helices (TMHMM PredHel, 0 when absent), is_gpi (PredGPI feature GPI-anchor), secreted_noTM (prob >= 0.9 and <= 1 TM helix) and secreted_strict (also no GPI anchor). Same rule as atlas job J4. Proteins without a signal peptide are not listed. Filter with species_prefix = <LOCUSTAG>.';

CREATE OR REPLACE VIEW v_genome_secretome AS
SELECT species_prefix AS LOCUSTAG,
       count(*) AS n_signalp,
       count(*) FILTER (WHERE signalp_prob >= 0.9) AS n_signalp_09,
       count(*) FILTER (WHERE secreted_noTM) AS n_secreted_noTM,
       count(*) FILTER (WHERE secreted_strict) AS n_secreted_strict,
       count(*) FILTER (WHERE is_gpi) AS n_signalp_gpi
FROM v_secretome_proteins
GROUP BY species_prefix;
COMMENT ON VIEW v_genome_secretome IS 'Per-genome secretome counts from v_secretome_proteins: n_signalp, n_signalp_09, n_secreted_noTM, n_secreted_strict, n_signalp_gpi. Equal to atlas A5/secretome.parquet for all 23,410 genomes (checked 2026-11-11). Genomes without SignalP rows are absent. Key LOCUSTAG.';

-- ---------------------------------------------------------------------------
-- Per-genome family counts for CAZy, TF and MEROPS
-- ---------------------------------------------------------------------------
-- Pattern: bin/atlas/sql/J4_function_counts.sql:55-100 and lib/bfd_data.py:338-359 build genome x family
-- counts; the analyses then read those atlas files (antarctolithica/scripts/04_build_features.py:46-50,
-- kickxellomycota/scripts/k6_pf.py:6). Gates are those of J4: CAZy toolcount >= 2, TF domain_i_evalue <= 1e-5,
-- MEROPS evalue <= 1e-10 and the best hit per protein by bitscore.
-- A Pfam count per genome is NOT here: it needs a scan of all of pfam (368 M rows). See extra_views_atlas.sql.

-- cazyme_fam looks like 'GH27(12-262)+CBM35(...)': one row per family token.
CREATE OR REPLACE VIEW v_cazy_family_counts AS
WITH u AS (
  SELECT species_prefix, protein_id, unnest(string_split(cazyme_fam, '+')) AS fam_tok
  FROM cazy_overview
  WHERE toolcount >= 2 AND cazyme_fam IS NOT NULL AND cazyme_fam <> '-'
)
SELECT species_prefix,
       regexp_extract(fam_tok, '^([A-Za-z]+)', 1) AS cazy_class,
       regexp_extract(fam_tok, '^([A-Za-z]+[0-9]+)', 1) AS cazy_family,
       count(DISTINCT protein_id) AS n_proteins
FROM u
GROUP BY 1, 2, 3;
COMMENT ON VIEW v_cazy_family_counts IS 'Proteins per genome and CAZy family (dbCAN calls with toolcount >= 2; cazyme_fam split on +). cazy_class is the letter prefix (GH, GT, AA, CBM, CE, PL). A family that is absent has no row (true zero only if v_genome_coverage.has_cazy). Equal to atlas A5/cazy_class_counts.parquet (checked 2026-11-11). Key species_prefix.';

CREATE OR REPLACE VIEW v_tf_family_counts AS
SELECT species_prefix, tf_family, count(DISTINCT protein_id) AS n_proteins
FROM tf_inventory
WHERE domain_i_evalue <= 1e-5
GROUP BY 1, 2;
COMMENT ON VIEW v_tf_family_counts IS 'Proteins per genome and transcription-factor family (tf_inventory, domain_i_evalue <= 1e-5). Equal to atlas A5/tf_counts.parquet (checked 2026-11-11). Key species_prefix.';

-- merops holds several BLAST hits per protein. One row per protein: the hit with the highest bitscore.
CREATE OR REPLACE VIEW v_merops_best AS
SELECT species_prefix, protein_id,
       arg_max(merops_id, bitscore) AS merops_id,
       max(bitscore) AS bitscore,
       arg_max(percent_identity, bitscore) AS percent_identity,
       arg_max(evalue, bitscore) AS evalue
FROM merops
WHERE evalue <= 1e-10
GROUP BY species_prefix, protein_id;
COMMENT ON VIEW v_merops_best IS 'Best MEROPS hit per protein (highest bitscore, evalue <= 1e-10). merops itself has several hits per protein, so count(*) on merops counts hits, not peptidases. Filter with species_prefix = <LOCUSTAG>. Key species_prefix, protein_id.';

-- One row per genome with the counts that the analyses join to the label table
-- (kickxellomycota/scripts/k1.py:4-11 joins labels + secretome + cazy_genome + merops_genome + tf_genome).
-- A count is NULL when the source table has no row for the genome (not run, or no positive protein).
CREATE OR REPLACE VIEW v_genome_function_counts AS
WITH sp AS MATERIALIZED (
  SELECT species_prefix, protein_id, signalp_prob, is_gpi, secreted_noTM, secreted_strict FROM v_secretome_proteins
), gs AS (
  SELECT species_prefix,
         count(*) AS n_signalp,
         count(*) FILTER (WHERE signalp_prob >= 0.9) AS n_signalp_09,
         count(*) FILTER (WHERE secreted_noTM) AS n_secreted_noTM,
         count(*) FILTER (WHERE secreted_strict) AS n_secreted_strict,
         count(*) FILTER (WHERE is_gpi) AS n_signalp_gpi
  FROM sp GROUP BY species_prefix
), cz AS (
  SELECT c.species_prefix, count(*) AS n_cazy, count(s.protein_id) AS n_cazy_secreted
  FROM (SELECT DISTINCT species_prefix, protein_id FROM cazy_overview WHERE toolcount >= 2) c
  LEFT JOIN sp s ON s.species_prefix = c.species_prefix AND s.protein_id = c.protein_id AND s.secreted_strict
  GROUP BY c.species_prefix
), me AS (
  SELECT m.species_prefix, count(*) AS n_merops, count(s.protein_id) AS n_merops_secreted
  FROM v_merops_best m
  LEFT JOIN sp s ON s.species_prefix = m.species_prefix AND s.protein_id = m.protein_id AND s.secreted_strict
  GROUP BY m.species_prefix
), tf AS (
  SELECT species_prefix, count(DISTINCT protein_id) AS n_tf_proteins
  FROM tf_inventory WHERE domain_i_evalue <= 1e-5 GROUP BY species_prefix
)
SELECT g.LOCUSTAG,
       gs.n_signalp, gs.n_signalp_09, gs.n_secreted_noTM, gs.n_secreted_strict, gs.n_signalp_gpi,
       cz.n_cazy, cz.n_cazy_secreted,
       me.n_merops, me.n_merops_secreted,
       tf.n_tf_proteins
FROM species g
LEFT JOIN gs ON gs.species_prefix = g.LOCUSTAG
LEFT JOIN cz ON cz.species_prefix = g.LOCUSTAG
LEFT JOIN me ON me.species_prefix = g.LOCUSTAG
LEFT JOIN tf ON tf.species_prefix = g.LOCUSTAG;
COMMENT ON VIEW v_genome_function_counts IS 'One row per genome: secretome counts (v_genome_secretome), n_cazy and n_cazy_secreted (dbCAN toolcount >= 2), n_merops and n_merops_secreted (best hit, evalue <= 1e-10), n_tf_proteins (domain_i_evalue <= 1e-5). Secreted = secreted_strict. NULL means the source table has no row for the genome. Same gates as atlas A5 (J4). Reads signalp, tmhmm, predgpi, cazy_overview, merops and tf_inventory in full (13-20 s with DuckDB 1.1.3 and 4 threads). Key LOCUSTAG.';

-- ---------------------------------------------------------------------------
-- v_gene_transcripts: gene_transcripts with boolean start/stop flags and the complete-model flag
-- ---------------------------------------------------------------------------
-- Pattern/trap: has_start_codon and has_stop_codon are VARCHAR with the values 'TRUE', 'FALSE' and 'NULL'
-- (the text NULL). Every analysis rewrites the same test:
--   (NOT is_partial AND has_start_codon = 'TRUE' AND has_stop_codon = 'TRUE') AS q_model
-- (bin/atlas/sql/J1_gene_stats.sql:21-26, J3_protein_core.sql:10-12, kickxellomycota/README.md:451;
-- kickxellomycota/README.md:58 lists the strings as a key-form trap).
-- has_start and has_stop are NULL where the table holds the text 'NULL'. q_model is false for those rows.
-- Filter with transcript_id = '...' or in_genome(transcript_id, '<LOCUSTAG>'), as for gene_transcripts.
CREATE OR REPLACE VIEW v_gene_transcripts AS
SELECT t.*,
       CASE has_start_codon WHEN 'TRUE' THEN true WHEN 'FALSE' THEN false END AS has_start,
       CASE has_stop_codon  WHEN 'TRUE' THEN true WHEN 'FALSE' THEN false END AS has_stop,
       (NOT is_partial AND has_start_codon = 'TRUE' AND has_stop_codon = 'TRUE') IS TRUE AS q_model
FROM gene_transcripts t;
COMMENT ON VIEW v_gene_transcripts IS 'gene_transcripts plus has_start and has_stop (BOOLEAN; NULL where the table holds the text NULL) and q_model (complete model: not partial, start and stop codon present; the atlas Q-model flag). In gene_transcripts itself has_start_codon and has_stop_codon are VARCHAR TRUE/FALSE/NULL. Filter with transcript_id = or in_genome(transcript_id, LOCUSTAG).';

-- ---------------------------------------------------------------------------
-- v_genome_truncation: per-genome partial-gene rate next to assembly contiguity
-- ---------------------------------------------------------------------------
-- Why (2026-11-16): after gene_stats was rebuilt with the upper-case fix (Fungi_BFD fb1d7cd), 0.14% of
-- 249,366,496 transcripts are partial. The genomes that still have a high partial rate are fragmented
-- assemblies, so the partial calls are real truncation at contig ends, not a soft-masking artifact:
--   partial rate   genomes   median contigs   median N50 (bp)   median masked_pct
--   <= 1%           22,976        208            329,909             10.4
--   1-5%               348      3,458             13,717             12.8
--   5-10%               58      4,761              9,264              9.6
--   > 10%               30      3,596              4,185              8.8
-- (Same 250 genomes before the fix: 14.23% partial; after: 0.19%.) This is a correlation across genomes.
-- Individual genes were not checked against contig ends.
-- truncation_class thresholds on pct_partial are a descriptive choice: high > 10, moderate > 5, low > 1, else none.
-- q_truncated = pct_partial > 5 (moderate or high): flag genomes whose gene counts and protein sets are
-- affected by truncation, e.g. when documenting problems in gene-level analyses.
-- Scans gene_transcripts (249 M rows, 3 columns): about a minute. One row per genome with transcripts
-- (23,412 on 2026-11-16; suppressed genomes are not in the tables).
CREATE OR REPLACE VIEW v_genome_truncation AS
WITH t AS (
  SELECT split_part(transcript_id, '_', 1) AS LOCUSTAG,
         count(*) AS n_transcripts,
         count(*) FILTER (WHERE is_partial) AS n_partial,
         count(*) FILTER (WHERE has_start_codon = 'FALSE') AS n_no_start,
         count(*) FILTER (WHERE has_stop_codon = 'FALSE') AS n_no_stop
  FROM gene_transcripts GROUP BY 1
), p AS (
  SELECT t.*, round(100.0 * n_partial / n_transcripts, 2) AS pct_partial FROM t
)
SELECT p.LOCUSTAG, s.ASMID, s.SPECIES, s.STRAIN, p.n_transcripts, p.n_partial, p.pct_partial,
       p.n_no_start, p.n_no_stop, a.contig_count, a.N50_bp, a.masked_pct,
       CASE WHEN p.pct_partial > 10 THEN 'high' WHEN p.pct_partial > 5 THEN 'moderate'
            WHEN p.pct_partial > 1 THEN 'low' ELSE 'none' END AS truncation_class,
       p.pct_partial > 5 AS q_truncated
FROM p LEFT JOIN species s ON s.LOCUSTAG = p.LOCUSTAG
       LEFT JOIN asm_stats a ON a.LOCUSTAG = p.LOCUSTAG;
COMMENT ON VIEW v_genome_truncation IS 'One row per genome: n_transcripts, n_partial, pct_partial, n_no_start, n_no_stop from gene_transcripts, plus contig_count, N50_bp, masked_pct. truncation_class = high (>10% partial), moderate (>5), low (>1), none. q_truncated = pct_partial > 5. High-partial genomes are fragmented assemblies (median N50 4,185 bp for >10%, 329,909 bp for <=1%); masked_pct is similar across classes, so the partials are real truncation, not a masking artifact. Scans gene_transcripts (about a minute).';

-- v_gene_introns_phase was removed on 2026-11-17: gene_introns.parquet was rebuilt on 2026-11-11 with Fungi_BFD
-- commit 34c8d14, so codon_position and codon_frame are correct in the table. Query gene_introns directly
-- (docs/BFD_VIEW_REVIEW.md section 4.3).
