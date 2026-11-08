-- MERGE_INFERNAL: <out>.rfam.tblout.gz (cmscan --cut_ga --fmt 2 --tblout, see
-- modules/BFD/INFERNAL/main.nf) -> infernal_rfam.parquet, one row per Rfam hit.
--
-- Fixed 2026-11-06: the former merge_infernal.py used the column positions of
-- --fmt 1, but the files are --fmt 2 (leading idx column, plus a clan column),
-- so every value after the first column was shifted (e.g. score held the GC
-- fraction, evalue the bias). It also cut the genome name at its first '.'
-- (Absidia_glauca_CBS_101.48... -> Absidia_glauca_CBS_101), so two genomes
-- could share a prefix. The genome column is now `genome_out`: the full file
-- name minus '.rfam.tblout.gz' (the `out` key, as in genome_annotation/).
--
-- --fmt 2 fields (1-based): 1 idx, 2 target name, 3 target accession, 4 query
-- name (contig), 5 query accession, 6 clan, 7 mdl, 8 mdl from, 9 mdl to,
-- 10 seq from, 11 seq to, 12 strand, 13 trunc, 14 pass, 15 gc, 16 bias,
-- 17 score, 18 E-value, 19 inc, 20 olp, 21 anyidx, 22 afrct1, 23 afrct2,
-- 24 winidx, 25 wfrct1, 26 wfrct2, 27 mdl len, 28 seq len, 29.. description.
COPY (
  WITH l AS (
    SELECT regexp_split_to_array(trim(line), '\s+') AS f,
           regexp_replace(parse_filename(filename), '\.rfam\.tblout(\.gz)?$', '') AS genome_out
    FROM read_csv({files}, columns = {'line': 'VARCHAR'}, header = false, delim = chr(1),
                  quote = '', escape = '', auto_detect = false, filename = true)
    WHERE line IS NOT NULL AND line NOT LIKE '#%' AND regexp_matches(line, '\S')
  )
  SELECT genome_out,
         f[2]  AS target_name,
         f[3]  AS target_acc,
         f[4]  AS contig,
         f[6]  AS clan,
         f[8]::INTEGER  AS mdl_from,
         f[9]::INTEGER  AS mdl_to,
         f[10]::BIGINT  AS seq_from,
         f[11]::BIGINT  AS seq_to,
         f[12] AS strand,
         f[13] AS trunc,
         f[15]::DOUBLE AS gc,
         f[16]::DOUBLE AS bias,
         f[17]::DOUBLE AS score,
         f[18]::DOUBLE AS evalue,
         f[19] AS inc,
         f[20] AS olp,
         f[27]::INTEGER AS mdl_len,
         f[28]::BIGINT  AS seq_len,
         array_to_string(f[29:], ' ') AS description
  FROM l
  WHERE len(f) >= 28
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD);
