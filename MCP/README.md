# BFD Protein/Genome Functional Database + MCP Server

Agentic investigation of protein composition, domain architecture, and genome characteristics across 22,412 fungal genomes.

---

## Overview

### Data pipeline
```
tables/*.csv.gz  →  scripts/build_duckdb.sh  →  functionalDB/function.duckdb
                                                          ↓
                                                    mcp_server.py  (stdio MCP)
                                                          ↓
                                               Claude / claude-code agent
```

### Scale
| Category | Count |
|---|---|
| Genomes (species) | 22,412 |
| Proteins | ~1.2 M |
| Pfam domain hits | ~1.2 M |
| CAZyme hits | ~26 K |
| Intergenic distances | ~38 M |

---

## 1. Build the Database

### Environment setup (one-time)
```bash
pixi install          # creates .pixi/envs/default/ with Python 3.11, duckdb, mcp
pixi run check        # verify imports
```

### Run the build
```bash
# Recommended: submit as a SLURM job (64 GB RAM, ~30–60 min)
sbatch scripts/build_duckdb.sh

# Or run interactively on a high-memory node
bash scripts/build_duckdb.sh
```

The script reads from `tables/` and writes `functionalDB/function.duckdb`.

### What the build script does
1. Loads `species.csv.gz` — taxonomy + LOCUSTAG identifier
2. Loads `asm_stats.tsv.gz` — assembly stats (TSV, joined to species on `ASMID` to add `LOCUSTAG`); column aliases: `gc_pct → GC_PERCENT`, `total_length_bp → TOTAL_LENGTH` for R-script compatibility
3. Loads all gene-structure tables: `gene_info`, `gene_proteins`, `gene_transcripts`, `gene_exons`, `gene_CDS`, `gene_introns`, `gene_trna`
4. Loads `gene_intergenic_distances` from `All_Taxa/` (~38 M rows)
5. Loads functional annotation tables: `pfam`, `cazy`, `cazy_overview`, `merops`
6. Loads subcellular targeting/secretome tables: `signalp`, `targetp`, `wolfpsort`, `predgpi`, `tmhmm`
7. Loads intrinsic disorder tables: `idp_summary`, `idp`
8. Loads composition tables: `aa_frequency`, `codon_frequency`
9. Creates two analytical views (see §3)

### Data quirks to know
| Issue | Detail |
|---|---|
| `LOCUSTAG` uniqueness | LOCUSTAG is the primary per-genome key and is enforced UNIQUE in the build. `ASMID` is also UNIQUE. |
| `ORDER` is a reserved SQL keyword | All references to the taxonomic ORDER column must be quoted: `"ORDER"` |
| `asm_stats` joins via ASMID | The raw TSV has no LOCUSTAG; it is added at load time by joining with species on ASMID |
| `statement_timeout` not supported | DuckDB does not expose this setting; row-cap via `LIMIT` in the MCP server is the only runaway-query protection |
| `gene_intergenic_distances` is large | 38 M rows; always filter by `species_prefix` before querying |

---

## 2. Join Key Reference

### Genome level
```
species.LOCUSTAG      ↔  asm_stats.LOCUSTAG
species.ASMID         ↔  asm_stats.ASMID           (unique join key)
species.LOCUSTAG      ↔  aa_frequency.species_prefix
species.LOCUSTAG      ↔  codon_frequency.species_prefix
```

### Species → proteins (all functional annotation tables)
```
species.LOCUSTAG  =  gene_proteins.species_prefix
                  =  pfam.species_prefix
                  =  cazy.species_prefix / cazy_overview.species_prefix
                  =  merops.species_prefix
                  =  signalp.species_prefix
                  =  targetp.species_prefix
                  =  wolfpsort.species_prefix
                  =  predgpi.species_prefix
                  =  tmhmm.species_prefix
                  =  idp_summary.species_prefix / idp.species_prefix
                  =  gene_intergenic_distances.species_prefix
```

### Protein → annotation
```
gene_proteins.protein_id  =  pfam.protein_id
                          =  cazy_overview.protein_id
                          =  merops.protein_id
                          =  signalp.protein_id
                          =  targetp.protein_id
                          =  wolfpsort.protein_id
                          =  predgpi.protein_id
                          =  tmhmm.protein_id
                          =  idp_summary.protein_id
```

---

## 3. Analytical Views

### `v_species_summary`
One row per species. Best starting point for cross-table investigation.

Columns: `LOCUSTAG, SPECIES, GENUS, FAMILY, "ORDER", CLASS, SUBPHYLUM, PHYLUM, BUSCO_LINEAGE, TOTAL_LENGTH, GC_PERCENT, contig_count, N50_bp, masked_pct, gene_count, mean_protein_length, pfam_proteins, cazy_proteins, merops_proteins, signalp_count, tmhmm_count, idp_high_disorder`

### `v_protein_annotation`
One row per protein with all annotation flags joined.

Columns: `protein_id, LOCUSTAG, length, pfam_domain_count, has_cazyme, has_merops, has_signal_peptide, targetp_class, wolfpsort_loc, tmhmm_helices, is_gpi_anchored, IDP_fraction`

---

## 4. MCP Server

### Start manually
```bash
PROTEIN_DB_PATH=functionalDB/function.duckdb \
  .pixi/envs/default/bin/python mcp_server.py
```

### Configure in claude-code
Add `.mcp.json` is already present in this directory. When you open claude-code inside this project, the `protein_db` MCP server is available automatically.

To add to your global Claude Desktop config (`~/.claude/claude_desktop_config.json`), copy the entry from `.mcp.json`.

### Tools

| Tool | Purpose |
|---|---|
| `db_execute_sql(sql, limit=1000)` | Run any read-only SELECT. Row cap enforced. The primary investigation tool. |
| `db_list_tables()` | All tables and views with row counts and descriptions |
| `db_describe_table(table_name)` | Column names, types, sample values for taxonomy/categorical columns |
| `db_sample_rows(table_name, n=5)` | Preview N rows |
| `db_pearson_correlation(x_col, y_col, from_sql, group_by?)` | Pearson r ± p-value, optionally stratified by taxonomy |
| `db_schema_overview()` | Full join-key map and table catalogue in Markdown |

### Security model
- DuckDB opened in **read-only** mode (`duckdb.connect(path, read_only=True)`)
- Write keywords (`INSERT`, `CREATE`, `DROP`, …) blocked by regex before execution
- Default row cap: 1,000 rows; maximum: 5,000 rows

---

## 5. Example Agent Questions

```
# Correlation investigation
"What is the Pearson correlation between genome GC content and the fraction
 of intrinsically disordered proteins, stratified by phylum?"

# Domain enrichment
"Which Pfam domains are most enriched in Mucoromycota compared to Ascomycota?"

# Secretome profile
"For species in the Polyporales order, what fraction of proteins have
 a signal peptide, and how does that correlate with CAZyme count?"

# Amino acid composition
"Do Basidiomycota genomes encode proteins with significantly different
 cysteine frequency compared to Ascomycota?"
```

---

## 6. R Script Compatibility

Existing R scripts in `scripts/Rscripts/` connect to `functionalDB/function.duckdb` and reference:
- `asm_stats.LOCUSTAG` — available (added at load time via ASMID join)
- `asm_stats.GC_PERCENT` — aliased from `gc_pct`
- `asm_stats.TOTAL_LENGTH` — aliased from `total_length_bp`

All existing R scripts should work against the newly built database without modification.

---

## 7. File Index

```
protein_mcp/
├── tables/                    # source CSVs (gzipped)
│   ├── species.csv.gz
│   ├── asm_stats.tsv.gz       # TSV not CSV; joins via ASMID
│   ├── gene_*.csv.gz
│   ├── pfam.csv.gz
│   ├── cazy*.csv.gz
│   ├── merops.csv.gz
│   ├── signalp.signal_peptide.csv.gz
│   ├── targetP.csv.gz
│   ├── wolfpsort.csv.gz
│   ├── predgpi.csv.gz
│   ├── tmhmm.csv.gz
│   ├── idp*.csv.gz
│   ├── aa_freq.csv.gz
│   ├── codon_freq.csv.gz
│   └── All_Taxa/
│       ├── aa_freq.csv.gz
│       ├── codon_freq.csv.gz
│       └── gene_intergenic_distances.csv.gz   # 38M rows
├── scripts/
│   ├── build_duckdb.sh        # SLURM build script
│   └── Rscripts/              # existing exploration scripts
├── functionalDB/
│   └── function.duckdb        # built by build_duckdb.sh
├── mcp_server.py              # FastMCP stdio server
├── pixi.toml                  # reproducible Python environment
└── .mcp.json                  # claude-code MCP configuration
```
