-- MERGE_ANTISMASH: per-genome antiSMASH-lite tables -> four Parquet tables (schema v2, 2026-11-22).
-- This file is the primary definition of the bgc_* tables. Run by bin/merge_antismash.py, which
-- fills the placeholders below (plain text replacement):
--   {regions} {cds} {known} {pairs}   SQL list literals of the input files (<LOCUSTAG>.antismash.<kind>.tsv.gz)
--   {topn}                            keep MIBiG hits with rank <= this (BFD default 5)
--   {out_bgc_cluster} {out_bgc_cluster_gene} {out_bgc_mibig_hit} {out_bgc_mibig_pair}   output paths (quoted)
-- The `genome` column of every input is the LOCUSTAG (the module passes meta.locustag), used here as `locustag`.
-- Join key between the tables: (locustag, record, region). `region` is unique only inside a record and `record`
-- only inside a genome, so always carry locustag. gene_id joins to gene_info.gene_id.
-- No sequences are stored. Genomes with no cluster give header-only files and no rows.
-- Sort orders (ORDER BY of each COPY) are part of the schema: they keep Parquet row-group statistics useful.
-- "end" is a reserved word and needs quotes.
CREATE TEMP TABLE r AS SELECT * FROM read_csv({regions}, delim='\t', header=true, auto_detect=true, types={'contig_edge':'BOOLEAN'});
CREATE TEMP TABLE c AS SELECT * FROM read_csv({cds}, delim='\t', header=true, auto_detect=true, types={'locus_tag':'VARCHAR','protein_id':'VARCHAR'});
CREATE TEMP TABLE k AS SELECT * FROM read_csv({known}, delim='\t', header=true, auto_detect=true);
CREATE TEMP TABLE p AS SELECT * FROM read_csv({pairs}, delim='\t', header=true, auto_detect=true, types={'query_locus_tag':'VARCHAR','mibig_protein':'VARCHAR'});
-- Best hit per region: highest similarity among ranks 1..topn; ties go to the lower rank.
CREATE TEMP TABLE best AS
  SELECT genome, record, region, arg_max(mibig_accession, (similarity, -rank)) AS best_mibig, max(similarity) AS best_similarity,
         arg_max(rank, (similarity, -rank)) AS best_rank
  FROM k WHERE rank <= {topn} GROUP BY genome, record, region;

-- 1. bgc_cluster: one row per antiSMASH region. Sorted by locustag, record, region.
--    best_* are NULL for regions without a MIBiG hit. Only ONE hit per region is summarised here; all hits are in bgc_mibig_hit.
COPY (
  SELECT r.genome AS locustag, r.record, r.region::USMALLINT AS region, r.start::UINTEGER AS start, r."end"::UINTEGER AS "end", r.length::UINTEGER AS length,
         r.products, r.categories, r.contig_edge, r.core_start::UINTEGER AS core_start, r.core_end::UINTEGER AS core_end,
         r.candidate_kinds, r.n_cds::USMALLINT AS n_cds, r.n_biosynthetic::USMALLINT AS n_biosynthetic, r.n_additional::USMALLINT AS n_additional,
         b.best_mibig, b.best_similarity::UTINYINT AS best_similarity, b.best_rank::UTINYINT AS best_rank
  FROM r LEFT JOIN best b ON b.genome=r.genome AND b.record=r.record AND b.region=r.region
  ORDER BY r.genome, r.record, r.region) TO {out_bgc_cluster} (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);

-- 2. bgc_cluster_gene: one row per gene inside a region. Sorted by locustag, record, region, start.
COPY (
  SELECT c.genome AS locustag, c.record, c.region::USMALLINT AS region, c.locus_tag AS gene_id, c.start::UINTEGER AS start, c."end"::UINTEGER AS "end",
         (c.strand='-') AS is_minus, c.gene_kind, c.functions
  FROM c ORDER BY c.genome, c.record, c.region, c.start) TO {out_bgc_cluster_gene} (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);

-- 3. bgc_mibig_hit: one row per kept MIBiG hit (rank 1..topn; rank is antiSMASH order, not similarity order).
--    Sorted by locustag, record, region, rank.
COPY (
  SELECT k.genome AS locustag, k.record, k.region::USMALLINT AS region, k.rank::UTINYINT AS rank, k.mibig_accession, k.hit_description AS compounds, k.hit_cluster_type,
         k.similarity::UTINYINT AS similarity, k.core_gene_hits::USMALLINT AS core_gene_hits, k.gene_hits::USMALLINT AS gene_hits, k.n_entry_proteins::USMALLINT AS n_entry_proteins,
         k.synteny_score::SMALLINT AS synteny_score, k.blast_score::FLOAT AS blast_score
  FROM k WHERE k.rank <= {topn} ORDER BY k.genome, k.record, k.region, k.rank) TO {out_bgc_mibig_hit} (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);

-- 4. bgc_mibig_pair: our gene -> MIBiG gene, one row per pair. Sorted by locustag, record, region, rank, gene_id.
--    Query coordinates and strand are not stored (join bgc_cluster_gene on gene_id); the MIBiG annotation text is not stored.
COPY (
  SELECT p.genome AS locustag, p.record, p.region::USMALLINT AS region, p.rank::UTINYINT AS rank, p.mibig_accession, p.query_locus_tag AS gene_id,
         p.mibig_gene, p.mibig_protein, p.perc_ident::UTINYINT AS perc_ident, p.perc_coverage::FLOAT AS perc_coverage, p.blast_score::FLOAT AS blast_score, p.evalue::DOUBLE AS evalue
  FROM p WHERE p.rank <= {topn} ORDER BY p.genome, p.record, p.region, p.rank, p.query_locus_tag) TO {out_bgc_mibig_pair} (FORMAT parquet, COMPRESSION zstd, ROW_GROUP_SIZE 100000);
