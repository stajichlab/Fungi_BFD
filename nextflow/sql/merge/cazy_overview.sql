-- MERGE_CAZY (overview): dbCAN <out>.overview.tsv.gz files -> cazy.overview.parquet
-- Same rows, columns and types as bin/merge_cazy.py (merge_overview) followed by
-- read_csv_auto(sample_size=-1). Each file has one header line (skipped).
-- Lines are read whole and split on TAB. Empty fields become NULL, as they do
-- when DuckDB re-reads the Python CSV.
COPY (
  WITH l AS (
    SELECT string_split(line, E'\t') AS f
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = true, delim = E'\x01',
                  quote = '', escape = '', auto_detect = false)
    WHERE line IS NOT NULL AND line <> ''
  )
  SELECT NULLIF(split_part(f[1], '_', 1), '') AS species_prefix,
         NULLIF(f[1], '') AS protein_id,
         NULLIF(f[2], '') AS EC,
         NULLIF(f[3], '') AS cazyme_fam,
         NULLIF(f[4], '') AS sub_fam,
         NULLIF(f[5], '') AS diamond_fam,
         NULLIF(f[6], '') AS substrate,
         NULLIF(f[7], '')::BIGINT AS toolcount
  FROM l
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
