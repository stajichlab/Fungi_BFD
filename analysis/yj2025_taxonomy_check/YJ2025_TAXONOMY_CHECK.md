# YJ2025 taxonomy check: "Saccharomycotina sp. YJ-2025a/b" are Pucciniomycotina (Basidiomycota)

Date: 2026-11-18. Status: complete. Fix: `data/curation/overrides.csv` and `samples.csv` (this branch).

## The problem

| ASMID | LOCUSTAG | NCBI name (taxid) | strain | BFD before | BFD after |
|---|---|---|---|---|---|
| GCA_054471075.1_ASM5447107v1 | FF90B97F | Saccharomycotina sp. YJ-2025a (3475847) | BHBT-Y41 | Ascomycota, no class | **Basidiomycota, Cystobasidiomycetes** |
| GCA_054471015.1_ASM5447101v1 | F6EB82EF | Saccharomycotina sp. YJ-2025b (3476451) | S2A6 | Ascomycota, no class | **Basidiomycota, Microbotryomycetes** |

NCBI Taxonomy (efetch, 2026-11-18) places both taxids in "Ascomycota; saccharomyceta; Saccharomycotina; unclassified
Saccharomycotina". Both assemblies were submitted by the Institute of Microbiology, Chinese Academy of Sciences, released
2026-01-12 (BioSamples SAMN52395896, SAMN52392630; BioProjects PRJNA1338122, PRJNA1338043). NCBI's assembly summary
flags both "unverified source organism; genus undefined" (PopSeq `catalog/kingdom/fungi_assemblies.csv.gz`).
BFD copies the NCBI lineage through `ncbi_accessions_taxonomy.csv` (taxonkit), which gives PHYLUM=Ascomycota and no class.

## Evidence

1. **Phylogenomics** (antarctolithica_genomics `analyses/V18_phylogenomics/README.md`, 154 taxa, ML and ASTRAL): both
   genomes sit inside Basidiomycota, next to *Rhodotorula pacifica*, with full support in every tree. That matrix has
   one Pucciniomycotina genome (*R. pacifica*) and one genome each of six other basidiomycete classes. So it shows
   Basidiomycota and, because they join the only Pucciniomycotina taxon rather than Agaricomycotina or
   Ustilaginomycotina, Pucciniomycotina. It cannot resolve the class.
2. **AAI against a Pucciniomycotina panel** (this folder; `run.sh`, SLURM job 32173766, 2026-11-18). One-way AAI from
   DIAMOND 2.1.7 `--sensitive` best hits (evalue <= 1e-5, alignment >= 50% of the query). Panel (72 genomes,
   `outputs/panel.tsv`): one annotated genome per genus for every Pucciniomycotina class in BFD (33 genera besides *Rhodotorula*, 9 classes),
   one per species of *Rhodotorula* (26 genomes, *R. pacifica* included), and 11 Ascomycota/Basidiomycota controls. *R. pacifica* (FBBC614A) was run as a
   third query to calibrate. Full table: `outputs/aai.tsv`.

Best AAI per class (excerpt of `outputs/aai.tsv`):

| reference class | YJ-2025a (FF90B97F) | YJ-2025b (F6EB82EF) |
|---|---|---|
| Cystobasidiomycetes | **58.1** (*Erythrobasidium hasegawianum*); ranks 1-4 all this class | 47.5 |
| Microbotryomycetes | 47.7 (*Glaciozyma antarctica*) | **59.5** (*Leucosporidium creatinivorum*); ranks 1-12 all this class |
| Agaricostilbomycetes | 47.3 | 47.0 |
| Mixiomycetes | 47.2 | 47.4 |
| other Pucciniomycotina classes | 45.2-46.6 | 45.6-46.8 |
| Agaricomycotina / Ustilaginomycotina controls | 43.6-44.9 | 44.5-45.6 |
| Saccharomycotina controls | 40.8-41.2 | 41.1-41.6 |
| each other (YJ-a vs YJ-b) | 46.8 | 47.5 |

Reading:

- Both are Basidiomycota and Pucciniomycotina: the Saccharomycotina controls are the most distant references (about 41%).
- **They are not the same lineage.** YJ-2025a is closest to Cystobasidiomycetes, YJ-2025b to Microbotryomycetes, with
  a gap of about 10 AAI points to the next class in each case. Their AAI to each other (47%) is the between-class level.
  V18 grouped both with *R. pacifica* only because it was the sole Pucciniomycotina genome in that matrix.
- Order and genus are not resolved. Within its class each genome is about equally far from several orders
  (YJ-2025b: Leucosporidiales 59.5, Kriegeriales 56.3, Microbotryales 55.7, Sporidiobolales 55.4; YJ-2025a:
  Erythrobasidiales 58.1, *Symmetrospora* 56.7, Cystobasidiales 56.0). For calibration, *R. pacifica* has AAI 63-77 to
  other *Rhodotorula* species; neither YJ genome reaches 60 to any reference, so neither belongs to a sampled genus.

Limits: one-way AAI loses resolution at deep divergences, so this is a class call, not a tree. The class gap (about
10 points) is far larger than the spread within classes. A tree with these Pucciniomycotina genera would confirm it and
might settle the orders. UNMEASURED: order placement.

## The fix

Ranks set: PHYLUM and CLASS only. SUBPHYLUM stays empty because it is empty for every BFD row (23,680 of 23,680 in
`Fungi_BFD_runs/samples.csv`); consumers derive it from CLASS (for example BFD_dashboard `bin/summary/taxonomy.py`), which
gives Pucciniomycotina for both. ORDER, FAMILY and GENUS stay empty (not resolved). SPECIES, SPECIES_IN and NCBI_TAXONID
keep the NCBI values, so the genomes stay joinable to NCBI; the misleading "Saccharomycotina sp." name remains until
NCBI changes it. BUSCO_LINEAGE (dikarya) and TRANSL_TABLE (1) are right for Basidiomycota and are unchanged.

`data/curation/overrides.csv` gets four rows (ASMID, FIELD, VALUE), so `scripts/create_samples_file.py` reproduces the
fix on the next regeneration. The two rows of the committed `samples.csv` are edited to match.

## Still to do (owner)

1. Merge this branch.
2. Apply the same two-row change to the run copy `Fungi_BFD_runs/samples.csv` (lines for FF90B97F and F6EB82EF:
   column 7 PHYLUM `Ascomycota` -> `Basidiomycota`, column 9 CLASS empty -> `Cystobasidiomycetes` / `Microbotryomycetes`).
3. Rebuild `Fungi_BFD_runs/tables/species.parquet` and `samples.parquet`, then BFD.duckdb (BFD_dashboard
   `workflows/build_bfd_db.sh`) and atlas labels (`bin/atlas/build_labels.py`).
4. Optionally tell NCBI / the submitter; a draft note is in the BFD_dashboard hand-off report of 2026-11-18.

## Reproduce

```bash
mkdir -p ~/projects/BFD_dashboard/bench/yj2025_taxonomy/scripts && cd ~/projects/BFD_dashboard/bench/yj2025_taxonomy
cp <repo>/analysis/yj2025_taxonomy_check/*.py <repo>/analysis/yj2025_taxonomy_check/run.sh scripts/
HERE=$PWD/scripts sbatch --export=ALL scripts/run.sh        # about 10 min on 32 cores
```

Test of the summariser: `python3 -m pytest analysis/yj2025_taxonomy_check/test_aai_summary.py` (synthetic hits).

| file | what |
|---|---|
| `select_panel.py` | picks the reference genomes from BFD tables and finds their proteomes |
| `run.sh` | SLURM job: DIAMOND databases, 213 pairwise searches, summary |
| `aai_summary.py` | one-way AAI per query x reference |
| `outputs/panel.tsv`, `outputs/nprot.tsv`, `outputs/aai.tsv` | panel, proteome sizes, results |
