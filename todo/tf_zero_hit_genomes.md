# T-038: 16 genomes with zero TF-domain hits

Added 2026-11-06. Status: open.

## Finding

`tables/tf_inventory.parquet` (standalone MERGE_TF, job 31806285, 2026-11-06)
covers 23,394 of the 23,410 Pfam `domtblout.gz` files. The other 16 genomes
have zero hits to the 25 TF Pfam domains in
`nextflow/assets/fungal_tf_pfam_domains.csv`.

This is not a key or prefix problem. Every file locus tag is unique, and every
`species_prefix` in the table matches a file.

These 16 files also have few Pfam hits overall (1 to 3,246 hit lines). One
reference genome (F01BDC31) has 9,367 lines. That is one data point, not a
measured norm.

Hypothesis (not verified): thin or failed gene predictions, not genomes that
truly lack TFs. Every fungal genome is expected to have Zn(II)2Cys6 and C2H2
TFs. 14 of 16 are Agaricomycetes (mostly ectomycorrhizal Boletaceae, *Amanita*,
*Tricholoma*). Several assembly names suggest metagenome or herbarium-specimen
sources. *Hamiltosporidium* is a microsporidian with a reduced proteome.

## To check for each genome

1. Total number of predicted proteins (`<LOCUSTAG>.proteins.fa` in input/pep,
   or `tables/gene_info.parquet` after the next MERGE_GENE_STATS).
2. BUSCO scores: genome mode (`tables/busco_genome.parquet`) and protein mode
   (`results/genome_stats/BUSCO_protein/`).
3. Compare 1 and 2 with genomes of the same genus or family.
4. Decide for each genome: re-predict, suppress, or keep with a flag.

## Genomes

| LOCUSTAG | Species | Strain | ASMID | Pfam hit lines | Proteins with Pfam hit | TF hits |
|---|---|---|---|---|---|---|
| F00081D8 | Amanita vidua | AV01 | GCA_059575585.1_UBL_Avid_1 | 3246 | 1389 | 0 |
| F1E1EE23 | aff. Porphyrellus sp. UT-M0000367 | UT-M0000367 | GCA_041722675.1_BD885 | 114 | 60 | 0 |
| F2E2A438 | Amanita subjunquillea | LH00001 | GCA_020011035.1_ASM2001103v1 | 3208 | 861 | 0 |
| F636EC73 | Calvatia gigantea | CGMCC 5.9 | GCA_025093545.1_ASM2509354v1 | 199 | 105 | 0 |
| F77192FB | Boletus austroedulis | NY1393612 | GCA_038089435.1_REH8497-Boletus-austroedulis | 272 | 203 | 0 |
| FA58BA6D | Boletaceae sp. KM206737 | KM206737 | GCA_041722175.1_MAD-m71-Xerocomus-sp | 107 | 86 | 0 |
| FAB5C0F4 | Lactarius quietus | S23C | GCA_015025155.1_Lacqui1 | 3009 | 2839 | 0 |
| FB348870 | Amanita bisporigera | Ingham01 | GCA_001983365.1_ASM198336v1 | 487 | 352 | 0 |
| FBDE1BFE | Russula nigricans | | GCA_964656355.1_gfRusNigr1.hap2.1 | 629 | 554 | 0 |
| FBFB4B38 | Hamiltosporidium tvaerminnensis | FI-OER-3-3 | GCA_022605425.2_FIOER33_v3 | 73 | 53 | 0 |
| FC0DB305 | Tricholoma bakamatsutake | SF-Tf05 | GCA_029959215.1_ASM2995921v1 | 331 | 305 | 0 |
| FCF521AA | Boletaceae sp. KM206743 | KM206743 | GCA_041722155.1_MAD-m63-Boletellus-sp | 327 | 283 | 0 |
| FCF75579 | Massospora cicadina | MCPNR19 | GCA_022478985.1_UCR_MCPNR19_1.0 | 1 | 1 | 0 |
| FD6195C5 | Tricholoma bakamatsutake | 20200826-1 | GCA_045514295.1_ASM4551429v1 | 290 | 265 | 0 |
| FD8664C9 | Austroboletus festivus | MICH11833 | GCA_038091375.1_MICH11833-Porphyrellus-festivus | 327 | 266 | 0 |
| FDAED364 | Boletaceae sp. KM206711 | KM206711 | GCA_041722215.1_MAD-m47-Boletaceae | 201 | 151 | 0 |
