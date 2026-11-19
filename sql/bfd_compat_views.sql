-- TABLES_DIR=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables
-- compat_views.sql: make SQL written for the MotherDuck database `bfd` run on the local BFD.duckdb.
-- Run after build_BFD_duckDB.sh:   duckdb BFD.duckdb < compat_views.sql
-- Every absolute parquet path below starts with the TABLES_DIR value on line 1. To move the
-- tables, replace that path in the whole file (for example with sed) before running it.
-- All statements are idempotent (CREATE OR REPLACE VIEW, COMMENT ON). The pipeline names
-- (aa_frequency, codon_frequency, gene_trna, signalp, cazy, targetp) are not changed.
-- Tested 2026-11-11 with DuckDB CLI 1.1.3 on a replica of the 2026-11-08 catalog, and with
-- python duckdb 1.5.5 (in-memory catalog, published BFD.duckdb attached READ_ONLY).

-- ---------------------------------------------------------------------------
-- 1. MotherDuck names that are missing in the pipeline database
-- ---------------------------------------------------------------------------
-- Aliases of existing views. The column lists are the same as the parquet files.
CREATE OR REPLACE VIEW aa_freq AS SELECT * FROM aa_frequency;
CREATE OR REPLACE VIEW codon_freq AS SELECT * FROM codon_frequency;
CREATE OR REPLACE VIEW signalp_signal_peptide AS SELECT * FROM signalp;
CREATE OR REPLACE VIEW cazy_cazymes_hmm AS SELECT * FROM cazy;
-- gene_trna adds species_prefix; gene_trnas.parquet (the MotherDuck table) does not have it.
CREATE OR REPLACE VIEW gene_trnas AS SELECT * EXCLUDE (species_prefix) FROM gene_trna;

-- Parquet files that have no view in the pipeline database.
CREATE OR REPLACE VIEW samples AS SELECT * FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/samples.parquet');
CREATE OR REPLACE VIEW swissprot_annot AS SELECT * FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/swissprot_annot.parquet');
CREATE OR REPLACE VIEW wgd_ks_nodes AS SELECT * FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/wgd.ks.nodes.parquet');
CREATE OR REPLACE VIEW infernal_rfam AS SELECT * FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/infernal_rfam.parquet');

-- ---------------------------------------------------------------------------
-- 2. Table descriptions (MotherDuck COMMENT ON TABLE texts from md_load_bfd.py DEFAULT_COMMENTS)
-- ---------------------------------------------------------------------------
-- 2a. Views created above.
COMMENT ON VIEW aa_freq IS 'Proteome amino-acid composition: one row per genome x amino acid, frequency = fraction of all residues. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW codon_freq IS 'Codon usage: one row per genome x codon, frequency = fraction of all codons in CDS. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW signalp_signal_peptide IS 'SignalP predicted signal peptides (secretion): peptide start/end, probability. Only proteins with a predicted signal peptide are listed. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW cazy_cazymes_hmm IS 'dbCAN HMM hits to CAZy families: HMM_id, profile/protein lengths, e-value, coordinates, coverage. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW gene_trnas IS 'tRNA genes: gene_id, amino acid, anticodon. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW samples IS 'Per-genome sample metadata: input species name, strain, BioProject, NCBI taxid, BUSCO lineage, translation table, taxonomy. Keys ASMID, LOCUSTAG.';
COMMENT ON VIEW swissprot_annot IS 'SwissProt annotation by accession: protein and gene names, organism, GO terms, EC numbers, InterPro, Pfam, length.';
COMMENT ON VIEW wgd_ks_nodes IS 'wgd gene-family tree nodes: pair counts, pass/sample flags, mean/min/max dS, mean identity. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW infernal_rfam IS 'PARTIAL COVERAGE: only 14 genomes as of 2026-11. Infernal cmscan hits to Rfam non-coding RNA models: model name/accession, contig, coordinates, strand, e-value, score. Genome key species_prefix = species.LOCUSTAG. DATA DEFECT seen 2026-11-11: species_prefix holds the sample name (for example ''Aaosphaeria_arxii_CBS_175''), not the LOCUSTAG, and the other columns are shifted by one (target_name is a row number, target_acc is the model name, query_name is the Rfam accession).';

-- 2b. Pipeline views that the build always creates (mkview ... 1, and the schema-contract views).
COMMENT ON VIEW species IS 'One row per genome. Taxonomy PHYLUM..SPECIES (SUBPHYLUM is mostly empty; use CLASS). Keys: ASMID; LOCUSTAG = protein_id prefix = species_prefix.';
COMMENT ON VIEW asm_stats IS 'Assembly statistics per genome (ASMID): contig count, total length, N50/L50, N90/L90, GC, gaps, repeat-masked bases, T2T scaffolds, telomere counts, annotation status. Local view columns: TOTAL_LENGTH (= total_length_bp in the parquet and in MotherDuck), GC_PERCENT (= gc_pct), plus LOCUSTAG; SPECIES and STRAIN are not in the view.';
COMMENT ON VIEW busco_genome IS 'BUSCO completeness per genome (ASMID): complete/single/duplicated/fragmented/missing %. All scored with fungi_odb12; budding-yeast classes (Saccharomycetes, Pichiomycetes) score ~75-80% from real gene loss, not poor assembly. This view adds LOCUSTAG.';
COMMENT ON VIEW aa_frequency IS 'Proteome amino-acid composition: one row per genome x amino acid, frequency = fraction of all residues. Genome key species_prefix = species.LOCUSTAG. MotherDuck name: aa_freq.';
COMMENT ON VIEW codon_frequency IS 'Codon usage: one row per genome x codon, frequency = fraction of all codons in CDS. Genome key species_prefix = species.LOCUSTAG. MotherDuck name: codon_freq.';
COMMENT ON VIEW gene_info IS 'One row per gene: coordinates (chrom, start, end, strand) and gene_type. LOCUSTAG joins species.';
COMMENT ON VIEW gene_transcripts IS 'One row per transcript: gene_id, coordinates, partial flag, start/stop codon presence, CDS length, md5. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW gene_proteins IS 'One row per protein: gene_id, protein_id, transcript_id, length (aa), md5 of sequence. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW gene_exons IS 'One row per exon: transcript_id, exon order, coordinates, GC content. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW gene_CDS IS 'One row per CDS segment: transcript_id, order, coordinates, GC content. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW gene_introns IS 'One row per intron: transcript_id, intron number, coordinates, 5''/3'' splice dinucleotides, GC. KNOWN BUG: codon_position and codon_frame describe the NEXT CDS junction (off by one exon; last intron gets protein length, frame 0) - recompute phase/position from gene_CDS cumulative lengths. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW gene_intergenic_distances IS 'Distance (bp) between each pair of adjacent genes (left_gene, right_gene). Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW gene_trna IS 'tRNA genes: gene_id, amino acid, anticodon. Genome = split_part(<id>,''_'',1) = species.LOCUSTAG. MotherDuck name: gene_trnas. This view adds species_prefix.';
COMMENT ON VIEW pfam IS 'Pfam hmmsearch domain hits, one row per domain: pfam_id (name), pfam_acc, e-values, scores, HMM/alignment/envelope coordinates. Genome = split_part(protein_id,''_'',1) = species.LOCUSTAG.';
COMMENT ON VIEW cazy_overview IS 'dbCAN CAZyme calls, one row per protein: CAZy family, subfamily, EC, predicted substrate, number of dbCAN tools agreeing (toolcount). Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW cazy IS 'dbCAN HMM hits to CAZy families: HMM_id, profile/protein lengths, e-value, coordinates, coverage. Genome key species_prefix = species.LOCUSTAG. MotherDuck name: cazy_cazymes_hmm.';
COMMENT ON VIEW merops IS 'MEROPS peptidase BLAST hits: merops_id, % identity, alignment stats, e-value, bitscore. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW signalp IS 'SignalP predicted signal peptides (secretion): peptide start/end, probability. Only proteins with a predicted signal peptide are listed. Genome key species_prefix = species.LOCUSTAG. MotherDuck name: signalp_signal_peptide.';
COMMENT ON VIEW targetp IS 'PARTIAL COVERAGE: only 125 genomes as of 2026-11 (check species_prefix before comparing taxa; use signalp_signal_peptide for secretion across ~23k genomes). TargetP targeting predictions: only proteins predicted SP (secretory signal peptide) or mTP (mitochondrial transit peptide) are listed; absent = no targeting peptide. Probability, cleavage site, motif. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW wolfpsort IS 'WoLF PSORT localization predictions: localization class and score per protein. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW predgpi IS 'PredGPI GPI-anchor predictions (GFF-style rows: feature, coordinates, score, attributes). Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW tmhmm IS 'TMHMM transmembrane helix predictions: only proteins with >=1 predicted helix are listed (about 20% of proteins); absent = no TM helix. Do not use for protein length of all proteins (use gene_proteins). Columns: len, ExpAA, First60, PredHel, Topology. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW idp IS 'Intrinsically disordered regions: one row per region with start, end, length, mean disorder score. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW idp_summary IS 'Per-protein disorder summary: disordered residues, IDP_fraction, protein length. Genome key species_prefix = species.LOCUSTAG.';

-- 2c. Pipeline objects that the build creates only when the parquet file exists (mkview ... 0 and the
--     SwissProt block). All of them exist in the 2026-11-08 database. If one is missing, its COMMENT
--     statement fails and the others still run; delete the line in that case.
COMMENT ON VIEW telomere_summary IS 'Telomere repeat summary per genome (ASMID, LOCUSTAG): scaffolds with telomeres at one/both ends, tract counts and lengths, 5''/3'' and terminal/internal counts, dominant repeat monomer.';
COMMENT ON VIEW telomere_tracts IS 'One row per telomeric repeat tract: scaffold, end (5''/3''), strand, monomer, repeat count, length, coordinates, whether terminal, distance to scaffold end, sequence. Keys ASMID, LOCUSTAG.';
COMMENT ON VIEW tf_inventory IS 'Transcription-factor inventory: proteins with TF-defining Pfam domains, assigned tf_family, domain e-value and score. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW wgd_ks IS 'Paralog pair Ks estimates (wgd): family, genes, dN, dS, dN/dS, alignment stats. dS is the synonymous divergence used for WGD peaks. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW wgd_ks_summary IS 'Paralog Ks summary per genome: pair/family/node counts, hybrid flag, identical-copy groups. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW swissprot_best IS 'Best SwissProt DIAMOND hit per protein, with subject annotation: accession, protein name, % identity/coverage x10 (divide by 10), bitscore, -log10 e-value. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW swissprot_hits IS 'All ranked SwissProt DIAMOND hits per protein (rk = rank). pident/qcov/scov stored x10. Join sp_id to swissprot_subject. Genome key species_prefix = species.LOCUSTAG.';
COMMENT ON VIEW swissprot_subject_summary IS 'Per SwissProt subject (sp_id): number of hits, number of genomes hit, best hits, max bitscore and identity (x10).';
-- swissprot_subject is a base table in the pipeline database, not a view.
COMMENT ON TABLE swissprot_subject IS 'SwissProt subject proteins referenced by hits (sp_id): accession, name, gene, organism, GO terms, EC, InterPro, Pfam.';

-- ---------------------------------------------------------------------------
-- 3. DECISION NEEDED, not active: asm_stats column names
-- ---------------------------------------------------------------------------
-- MotherDuck asm_stats was a copy of asm_stats.parquet (columns SPECIES, STRAIN, total_length_bp, gc_pct).
-- The pipeline view renames total_length_bp -> TOTAL_LENGTH and gc_pct -> GC_PERCENT, adds LOCUSTAG and
-- drops SPECIES and STRAIN (sql/table_schema.json). SQL that names total_length_bp or gc_pct on asm_stats
-- fails locally. The statement below keeps every current column and adds the two parquet names. It passed
-- the same tests as the statements above. It is commented out because asm_stats is a schema-contract view:
-- the durable fix is two more columns in sql/table_schema.json, not an override here.
-- SPECIES and STRAIN are left out on purpose: they would make `SPECIES` ambiguous in existing queries that
-- join species and asm_stats.
--
-- CREATE OR REPLACE VIEW asm_stats AS
-- SELECT CAST(species.LOCUSTAG AS VARCHAR) AS LOCUSTAG, CAST(t.ASMID AS VARCHAR) AS ASMID,
--        CAST(t.contig_count AS BIGINT) AS contig_count, CAST(t.total_length_bp AS BIGINT) AS TOTAL_LENGTH,
--        CAST(t.min_contig_bp AS BIGINT) AS min_contig_bp, CAST(t.max_contig_bp AS BIGINT) AS max_contig_bp,
--        CAST(t.median_contig_bp AS DOUBLE) AS median_contig_bp, CAST(t.mean_contig_bp AS DOUBLE) AS mean_contig_bp,
--        CAST(t.L50 AS BIGINT) AS L50, CAST(t.N50_bp AS BIGINT) AS N50_bp, CAST(t.L90 AS BIGINT) AS L90,
--        CAST(t.N90_bp AS BIGINT) AS N90_bp, CAST(t.gc_pct AS DOUBLE) AS GC_PERCENT,
--        CAST(t.n_gap_count AS BIGINT) AS n_gap_count, CAST(t.total_n_bases AS BIGINT) AS total_n_bases,
--        CAST(t.masked_bases AS BIGINT) AS masked_bases, CAST(t.masked_pct AS DOUBLE) AS masked_pct,
--        CAST(t.t2t_scaffolds AS BIGINT) AS t2t_scaffolds, CAST(t.telomere_fwd AS BIGINT) AS telomere_fwd,
--        CAST(t.telomere_rev AS BIGINT) AS telomere_rev, CAST(t.annotation_status AS VARCHAR) AS annotation_status,
--        CAST(t.total_length_bp AS BIGINT) AS total_length_bp, CAST(t.gc_pct AS DOUBLE) AS gc_pct
-- FROM read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/asm_stats.parquet') AS t
-- LEFT JOIN read_parquet('/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/tables/species.parquet') AS species ON t.ASMID = species.ASMID;
