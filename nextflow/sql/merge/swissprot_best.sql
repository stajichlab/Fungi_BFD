-- MERGE_SWISSPROT_GATHER: best SwissProt hit per BFD protein (rk = 1) with the
-- subject strings joined in, so most queries need no join. {files} = the
-- swissprot_hits.<bucket>.parquet files; {subject} = swissprot_subject.parquet.
COPY (
  SELECT h.protein_id, h.species_prefix, h.sp_id, s.accession, s.entry_name, s.org_code,
         s.protein_name, h.pident_x10, h.qcov_x10, h.scov_x10, h.aln_length,
         h.q_start, h.q_end, h.s_start, h.s_end, h.q_len, s.s_len,
         h.bitscore, h.neg_log10_evalue
  FROM read_parquet({files}) h JOIN read_parquet({subject}) s USING (sp_id)
  WHERE h.rk = 1
  ORDER BY h.protein_id
) TO {out} (FORMAT PARQUET, COMPRESSION ZSTD, ROW_GROUP_SIZE 1000000);
