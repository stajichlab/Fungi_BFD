# T-038: gene prediction failed for genomes with good assemblies

Added 2026-11-06. Status: in progress. Cause: not verified.

## How it was found

`tables/tf_inventory.parquet` (standalone MERGE_TF, job 31806285, 2026-11-06)
covers 23,394 of the 23,410 Pfam `domtblout.gz` files. The other 16 genomes
have zero hits to the TF Pfam domains in
`nextflow/assets/fungal_tf_pfam_domains.csv`. This is not a key or prefix
problem: every file locus tag is unique, and every `species_prefix` in the
table matches a file.

## The 16 zero-TF genomes

All 16 have protein BUSCO complete <= 0.7% (fungi_odb12), so gene prediction
failed in every one. Sources: `input/pep/*.proteins.fa.gz`,
`results/genome_stats/BUSCO_{genome,protein}/`, `tables/asm_stats.parquet`.
The full table is in `Fungi_BFD_runs/BFD_function/run_merge_tf.sh` (comment
block at the end).

| Group | Genomes | Meaning |
|---|---|---|
| A | F00081D8, F2E2A438, F636EC73, FAB5C0F4, FBDE1BFE, FC0DB305, FD6195C5 | Good assembly (gBUSCO 97.4-99.8%), failed annotation |
| B | F1E1EE23, F77192FB, FB348870, FDAED364, FCF521AA, FA58BA6D, FD8664C9 | Fragmented or partial assembly (gBUSCO 6.7-74.9%, N50 2.8-14.9 kb) |
| C | FBFB4B38 (microsporidian), FCF75579 (Massospora, 1.49 Gb) | Special cases |

## Scope: 52 genomes

Scan of 23,433 genomes with both BUSCO results (suppressed genomes excluded):

| Genome BUSCO | pBUSCO < 10% | 10-50% | 50-80% | >= 80% |
|---|---|---|---|---|
| >= 90% | 26 | 26 | 117 | 15,599 |
| < 90% | 145 | 422 | 5,366 | 1,732 |

52 genomes have gBUSCO >= 90% and pBUSCO < 50%: Agaricomycetes 28,
Sordariomycetes 13, Pezizomycetes 5, Tremellomycetes 3, Eurotiomycetes 2,
Leotiomycetes 1. The worst are mostly ectomycorrhizal (Amanita x7, Lactarius,
Lactifluus, Russula, Tricholoma, boletes, Tuber). Group A is in the 52; groups
B and C are not (gBUSCO < 90%). The 145 genomes with gBUSCO < 90% and
pBUSCO < 10% mix this failure with poor assemblies; BUSCO alone cannot
separate the two causes.

## Stage trace: Amanita vidua F00081D8 (+3 Group A genomes for GeneMark)

1. No RNA-seq or PASA (no `training/` folder in any of the 16). Augustus and
   SNAP were trained from BUSCO models.
2. Training data was fine: 998 complete single-copy BUSCOs validated;
   Augustus trained on 936 normal multi-exon genes.
3. Raw GeneMark (from `genemark_store/`, self-trained per genome, EVM weight 4,
   the highest) was already poor:

   | Genome | GeneMark genes | Assembly | Median aa (GeneMark -> final) |
   |---|---|---|---|
   | Calvatia gigantea | 1,648 | 40 Mb | 328 -> 300 |
   | Russula nigricans | 4,241 | 45 Mb | 237 -> 126 |
   | Tricholoma bakamatsutake SF-Tf05 | 12,473 | 142 Mb | 271 -> 168 |
   | Amanita vidua | 9,625 | 56 Mb | 290 -> 260 |

4. EVM made 9,235 models; the TE filter removed 3,379 (37%) of them as
   transposable elements; 5,834 remained.
5. Final proteins have no internal stops and start with Met, but 20-40% are
   < 100 aa (reference F01BDC31: 3.9%). They are fragments.

Intermediate files (raw Augustus, EVM input/output) were deleted by the
`predict_misc/` trim, so steps 3-4 could not be traced fully.

Hypothesis (not verified): GeneMark-ES self-training gives poor models on these
genomes, EVM trusts GeneMark most, and the TE filter removes many of the rest.
Ectomycorrhizal genomes are large and TE-rich, which fits.

## Diagnostic re-run (set up 2026-11-06, by another session)

Folder `Fungi_BFD_runs/do_annotation_busco_gap/`: isolated `genome_annotation/`,
`genemark_store/` and ab initio stores; no RNA-seq; no model sharing or reuse;
fresh GeneMark; `predict_keep_misc: true` (branch `predict-keep-misc`, worktree
`Fungi_BFD_wt_keepmisc`). First run: `samples_busco_gap_diag3.csv`
(Russula nigricans, Amanita vidua, Tuber melanosporum). All 52:
`samples_busco_gap.csv`.

## Next steps

1. Run the 3-genome diagnostic. Run protein BUSCO on each stage (Augustus,
   GeneMark, EVM, post-TE-filter, final) to find where genes are lost.
2. Test variants: lower GeneMark EVM weight (`predict_evm_weights_norna`), and
   TE filter off (funannotate `--repeat_filter`; no pipeline param yet).
3. Apply the fix to the 52, then decide groups B and C (re-predict, suppress,
   or flag).

## The 52 genomes (gBUSCO >= 90%, pBUSCO < 50%), sorted by protein BUSCO

| LOCUSTAG | Species | Strain | ASMID | Genome BUSCO C % | Protein BUSCO C % | Class |
|---|---|---|---|---|---|---|
| FD6195C5 | Tricholoma bakamatsutake | 20200826-1 | GCA_045514295.1_ASM4551429v1 | 99.8 | 0.1 | Agaricomycetes |
| FC0DB305 | Tricholoma bakamatsutake | SF-Tf05 | GCA_029959215.1_ASM2995921v1 | 99.7 | 0.1 | Agaricomycetes |
| FAB5C0F4 | Lactarius quietus | S23C | GCA_015025155.1_Lacqui1 | 98.0 | 0.2 | Agaricomycetes |
| FE1ED680 | Lactarius vividus | 141 | GCA_021527035.1_Lacviv1 | 98.8 | 0.2 | Agaricomycetes |
| FBDE1BFE | Russula nigricans |  | GCA_964656355.1_gfRusNigr1.hap2.1 | 99.2 | 0.2 | Agaricomycetes |
| FA685757 | Amanita molliuscula | LH00009 | GCA_023014405.1_ASM2301440v1 | 96.0 | 0.3 | Agaricomycetes |
| F3C410A9 | Amanita rimosa | LH00003 | GCA_023014395.1_ASM2301439v1 | 98.3 | 0.3 | Agaricomycetes |
| F2E2A438 | Amanita subjunquillea | LH00001 | GCA_020011035.1_ASM2001103v1 | 97.4 | 0.3 | Agaricomycetes |
| F00081D8 | Amanita vidua | AV01 | GCA_059575585.1_UBL_Avid_1 | 98.1 | 0.4 | Agaricomycetes |
| FA9FBCE9 | Entoloma clypeatum |  | GCA_963565915.2_gfEntClyp1.2 | 98.4 | 0.4 | Agaricomycetes |
| F2C21852 | Tuber magnatum | 091103-1 | GCA_003182015.1_ASM318201v1 | 95.1 | 0.4 | Pezizomycetes |
| FC60ED31 | Pulveroboletus aff. annulatus | UT-M0000386 | GCA_038094595.1_BD784-Pulveroboletus-sp1 | 98.4 | 0.5 | Agaricomycetes |
| F7E3887E | Amanita fuliginea | LH00002 | GCA_023014435.1_ASM2301443v1 | 98.5 | 0.6 | Agaricomycetes |
| F2288138 | Lactarius akahatsu | QP | GCA_021524915.1_Lacaka1 | 98.8 | 0.6 | Agaricomycetes |
| F636EC73 | Calvatia gigantea | CGMCC 5.9 | GCA_025093545.1_ASM2509354v1 | 98.0 | 0.7 | Agaricomycetes |
| FB2B1CC7 | Tuber brumale | 140227-1 | GCA_014065205.1_ASM1406520v1 | 98.5 | 0.8 | Pezizomycetes |
| F5668708 | Austroboletus festivus | G4425 | GCA_038088875.1_TH9250-Austroboletus-festivus | 91.4 | 0.9 | Agaricomycetes |
| F0A03C30 | Amanita exitialis | LH00005 | GCA_023014455.1_ASM2301445v1 | 98.1 | 1.4 | Agaricomycetes |
| F7BFE534 | Lanmaoa asiatica | LA01 | GCA_020310465.1_ASM2031046v1 | 98.9 | 1.7 | Agaricomycetes |
| F09DACD7 | Lactarius evosmus |  | GCA_951394315.1_gfLacEvos1.1 | 99.1 | 2.0 | Agaricomycetes |
| F09BA05D | Russula ochroleuca | Prilba | GCA_015178965.1_Rusoch1 | 98.5 | 2.1 | Agaricomycetes |
| FFBE8F54 | Tuber melanosporum | Mel28 | GCF_000151645.1_ASM15164v1 | 98.2 | 2.6 | Pezizomycetes |
| F236FCE0 | Tuber calosporum | MG102 | GCA_003316355.1_ASM331635v1 | 96.4 | 3.0 | Pezizomycetes |
| F25E1C95 | Lactifluus vellereus |  | GCA_984589725.1_gfLacVell1.hap1.1 | 98.4 | 4.9 | Agaricomycetes |
| FF9A3D86 | Gautieria morchelliformis | GMNE.BST | GCA_015179125.1_Gaumor1_1 | 96.7 | 6.0 | Agaricomycetes |
| F01E5862 | Fusarium circinatum | S | GCA_013168835.2_ASM1316883v2 | 90.5 | 9.2 | Sordariomycetes |
| FE33BEFC | Macrolepiota albuminosa | Ter1-01 | GCA_047678695.1_ASM4767869v1 | 94.0 | 11.7 | Agaricomycetes |
| F364B67D | Leucocoprinus sp. KB1807203 | KB1807203 | GCA_019327925.1_ASM1932792v1 | 98.3 | 13.5 | Agaricomycetes |
| F997B409 | Rhodocollybia maculata |  | GCA_965654425.1_gfRhoMacu1.hap1.1 | 98.4 | 14.0 | Agaricomycetes |
| F3234657 | Boletellus aff. lepidospora | UT-M0000424 | GCA_038094175.1_BD832-Boletellus-aff-lepidosporta | 95.2 | 18.1 | Agaricomycetes |
| F37DB1A3 | Beauveria bassiana | BbFZ51 | GCA_044509585.1_ASM4450958v1 | 98.4 | 22.5 | Sordariomycetes |
| FA14FB9E | Verticillium dahliae | V13 | GCA_051257765.1_ASM5125776v1 | 96.2 | 25.0 | Sordariomycetes |
| FFC954C5 | Monilinia fructigena | KZ-E9_2025 | GCA_052815175.1_IPBB_Mfgen_1 | 96.7 | 28.3 | Leotiomycetes |
| FC6905A2 | Beauveria pseudobassiana | RGM 2184 | GCA_022985055.1_ASM2298505v1 | 98.3 | 29.5 | Sordariomycetes |
| FB7DB937 | Verticillium dahliae | Vd-653 | GCA_050229775.1_ASM5022977v1 | 96.1 | 30.2 | Sordariomycetes |
| FA1C6A97 | Cryptococcus neoformans | VYD146 | GCA_011057555.1_ASM1105755v1 | 94.7 | 33.5 | Tremellomycetes |
| FE67B7D3 | Pyricularia oryzae | N06047 | GCA_911175055.1_N06047_contigs | 94.1 | 33.7 | Sordariomycetes |
| FD011935 | Cryptococcus neoformans | VYD147 | GCA_011057515.1_ASM1105751v1 | 95.3 | 33.8 | Tremellomycetes |
| F808427B | Metarhizium pingshaense | WP08 | GCA_049901535.1_ASM4990153v1 | 98.7 | 36.2 | Sordariomycetes |
| F8EAB1C9 | Phlegmacium triumphans |  | GCA_965275855.1_gfCorTriu1.hap1.1 | 99.4 | 38.4 | Agaricomycetes |
| F9DD1BDF | Microsporum canis | BMU 11032 | GCA_026259295.1_ASM2625929v1 | 98.0 | 38.7 | Eurotiomycetes |
| F91F5D20 | Lactarius indigo | MG109 | GCA_003313985.1_ASM331398v1 | 91.3 | 39.5 | Agaricomycetes |
| FA8D96BC | Verticillium dahliae | HOMCF | GCA_050229575.1_ASM5022957v1 | 98.0 | 42.4 | Sordariomycetes |
| F1B5F3CC | Verticillium dahliae | Vd-H5 | GCA_050229755.1_ASM5022975v1 | 98.3 | 43.6 | Sordariomycetes |
| FD4D4E23 | Trichophaea hybrida | UTF0779 | GCA_015178995.1_Trihyb1 | 94.9 | 43.9 | Pezizomycetes |
| F178DE30 | Verticillium dahliae | Vd-R1 | GCA_050229925.1_ASM5022992v1 | 98.5 | 44.7 | Sordariomycetes |
| FB678174 | Verticillium dahliae | I1V | GCA_051257705.1_ASM5125770v1 | 96.3 | 45.2 | Sordariomycetes |
| F7B7E80F | Trichophyton interdigitale | IGIB-BK-CI04 | GCA_037576225.1_ASM3757622v1 | 98.7 | 47.2 | Eurotiomycetes |
| F81DC707 | Cryptococcus neoformans x Cryptococcus deneoformans | CBS 132 | GCA_006992865.1_ASM699286v1 | 95.5 | 48.3 | Tremellomycetes |
| FFDE7030 | Rosellinia necatrix | R10 | GCA_021209715.1_ASM2120971v1 | 98.4 | 48.7 | Sordariomycetes |
| FE37F312 | Verticillium dahliae | VT-2A | GCA_050229655.1_ASM5022965v1 | 98.3 | 49.1 | Sordariomycetes |
| F3E6C2A1 | Macrolepiota albuminosa | NeiJiang_2021 | GCA_027475525.1_NeiJiang_Talb_v1.0_genomic | 90.7 | 49.6 | Agaricomycetes |
