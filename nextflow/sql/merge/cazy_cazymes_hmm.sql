-- MERGE_CAZY (HMM hits): dbCAN <out>.cazymes.tsv.gz files -> cazy.cazymes_hmm.parquet
-- Same rows, columns and types as bin/merge_cazy.py (merge_hmm) followed by
-- read_csv_auto(sample_size=-1). Each file has one header line (skipped).
-- species_prefix is the part of Gene_ID (column 3) before the first '_'.
COPY (
  WITH l AS (
    SELECT string_split(line, E'\t') AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = true, delim = E'\x01',
                  quote = '', escape = '', auto_detect = false)
    WHERE line IS NOT NULL AND line <> ''
  )
  SELECT NULLIF(split_part(f[3], '_', 1), '') AS species_prefix,
         NULLIF(f[1], '') AS HMM_id,
         NULLIF(f[2], '')::BIGINT AS profile_length,
         NULLIF(f[3], '') AS protein_id,
         NULLIF(f[4], '')::BIGINT AS protein_length,
         NULLIF(f[5], '')::DOUBLE AS evalue,
         NULLIF(f[6], '')::BIGINT AS q_start,
         NULLIF(f[7], '')::BIGINT AS q_end,
         NULLIF(f[8], '')::BIGINT AS s_start,
         NULLIF(f[9], '')::BIGINT AS s_end,
         NULLIF(f[10], '')::DOUBLE AS coverage
  FROM l
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
