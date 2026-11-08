-- MERGE_SWISSPROT_GATHER: one row per SwissProt subject that BFD proteins hit,
-- for subject-first questions (how many proteins / genomes hit this entry, how
-- often it is the best hit). {files} = swissprot_hits.<bucket>.parquet files.
COPY (
  SELECT sp_id,
         count(*) AS n_hits,
         count(DISTINCT species_prefix) AS n_genomes,
         sum((rk = 1)::INTEGER) AS n_best_hits,
         max(bitscore) AS max_bitscore,
         max(pident_x10) AS max_pident_x10
  FROM read_parquet({files})
  GROUP BY sp_id
  ORDER BY sp_id
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
