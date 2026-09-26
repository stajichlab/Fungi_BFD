## Table M1. Locus rule for one PASA model per locus (step B with each rule's own PASA evidence)

Holdout-chromosome gffcompare, CDS level. Values are locus Sn / Pr (%), predicted genes.

| Genome | tx_strand | cds_strand | cds_blind |
|---|---|---|---|
| N. crassa OR74A | 63.8 / 71.8 (3447) | 63.5 / 70.9 (3480) | 63.5 / 70.8 (3478) |
| A. nidulans FGSC A4 | 54.7 / 55.7 (4911) | 54.7 / 55.1 (4975) | 54.7 / 55.1 (4975) |
| B. cinerea B05.10 | 82.9 / 81.4 (5711) | 83.2 / 80.3 (5811) | 83.2 / 80.2 (5811) |
| C. neoformans H99 | 80.5 / 81.4 (2691) | 79.8 / 78.8 (2754) | 79.9 / 79.0 (2753) |
| S. commune H4-8 | 25.0 / 36.0 (4984) | 25.1 / 35.6 (5060) | 25.1 / 35.6 (5059) |

## Table M2. Training-set construction, same fixed EVM evidence (training effect only)

Locus Sn / Pr (%); training models in parentheses. old = funannotate 41a2fd7; oldnew = old getBestModel + new selection; tx = new getBestModel (tx_strand, complete-first) + new selection; busco = BUSCO-forced training.

| Genome | old | oldnew | tx | busco |
|---|---|---|---|---|
| N. crassa OR74A | 63.3 / 71.6 (418) | 63.5 / 71.8 (418) | 64.3 / 72.6 (465) | 66.0 / 74.1 |
| A. nidulans FGSC A4 | 54.6 / 55.7 (1103) | 54.7 / 55.8 (1103) | 54.6 / 55.8 (1122) | 54.8 / 56.8 |
| B. cinerea B05.10 | 83.3 / 82.2 (1229) | 83.1 / 82.0 (1229) | 83.3 / 82.2 (1284) | 82.0 / 82.0 |
| C. neoformans H99 | 80.1 / 80.9 (814) | 80.1 / 81.0 (814) | 79.9 / 80.8 (832) | 78.1 / 80.6 |
| S. commune H4-8 | 29.2 / 44.2 (552) | 24.5 / 35.5 (87) | 24.9 / 36.0 (92) | 37.4 / 49.2 |

## Table M3. Ranking inside a locus: complete-first vs guarded (complete only if CDS ≥ 80% of the locus's longest)

| Genome | Evidence | complete-first locus Sn / Pr (train models) | guarded locus Sn / Pr (train models) |
|---|---|---|---|
| N. crassa OR74A | fixed | 64.3 / 72.6 (465) | 63.9 / 71.9 (441) |
| N. crassa OR74A | own | 63.8 / 71.8 (465) | 63.6 / 71.6 (441) |
| A. nidulans FGSC A4 | fixed | 54.6 / 55.8 (1122) | 54.6 / 55.8 (1111) |
| A. nidulans FGSC A4 | own | 54.7 / 55.7 (1122) | 54.7 / 55.8 (1111) |
| B. cinerea B05.10 | fixed | 83.3 / 82.2 (1284) | 83.3 / 82.1 (1251) |
| B. cinerea B05.10 | own | 82.9 / 81.4 (1284) | 83.3 / 82.0 (1251) |

## Table M4. Single-exon training genes (R6 option b): se (on) vs tx2 (off), same code, fixed evidence

Exact CDS-chain match to RefSeq protein-coding mRNAs on holdout chromosomes. Sn per RefSeq gene, Pr per predicted model.

| Genome | RefSeq single / multi genes | Single Sn | Single Pr | Multi Sn | Multi Pr |
|---|---|---|---|---|---|
| N. crassa OR74A | 877 / 3035 | 49.4 → 54.7 | 65.9 → 60.6 | 59.9 → 59.0 | 65.4 → 67.8 |
| A. nidulans FGSC A4 | 661 / 4353 | 71.3 → 78.5 | 58.6 → 55.6 | 47.7 → 47.4 | 50.6 → 51.8 |
| B. cinerea B05.10 | 1236 / 4484 | 64.2 → 71.8 | 74.0 → 73.2 | 79.5 → 79.2 | 77.2 → 79.0 |
| C. neoformans H99 | 77 / 2649 | 50.6 → 55.8 | 34.5 → 32.8 | 74.2 → 74.2 | 76.2 → 76.6 |
| S. commune H4-8 | 1237 / 5988 | 7.6 → 7.3 | 29.0 → 28.3 | 23.7 → 23.8 | 30.5 → 30.6 |

## Table M5. PASA training vs BUSCO training, divergent reads (N. crassa OR74A; RNA-seq from strain HJDF, median read identity 96.7%)

| Arm | PASA input | Single-exon training | Training | Fixed evidence locus Sn / Pr | Own evidence locus Sn / Pr |
|---|---|---|---|---|---|
| tx | rc1 (old minimap2 parser) | no | PASA | 64.3 / 72.6 | 63.8 / 71.8 |
| se | rc1 | yes | PASA | 65.1 / 73.6 | n/a |
| txR1 | R1 | no | PASA | 64.4 / 72.7 | 65.4 / 73.1 |
| txR1R2 | R1 + R2 | no | PASA | 64.6 / 73.0 | 65.8 / 73.6 |
| txID90 | R1 + gmap, 90% identity | no | PASA | 64.5 / 72.6 | 64.6 / 71.9 |
| txRel | R1, relaxed validation | no | PASA | 64.5 / 72.6 | 65.8 / 72.9 |
| seR1R2 | R1 + R2 | yes | PASA | 66.4 / 74.9 | 67.1 / 75.0 |
| busco | rc1 (evidence only) | n/a | BUSCO | 66.0 / 74.1 | n/a |
| buscoR1R2 | R1 + R2 (evidence only) | n/a | BUSCO | 66.1 / 74.2 | 66.8 / 74.4 |

## Table M6. Genome with few complete PASA models (S. commune H4-8: 93 complete of 822 PASA models, 30 keepers; divergent reads)

| Arm | Locus Sn / Pr | Intron-chain Sn / Pr | Exon Sn / Pr | Proteome BUSCO C (%) |
|---|---|---|---|---|
| busco | 37.4 / 49.2 | 42.3 / 48.6 | 69.8 / 79.6 | 97.5 |
| old | 29.2 / 44.2 | 34.0 / 43.8 | 56.4 / 78.1 | 86.9 |
| txR1R2 | 28.8 / 39.2 | 32.5 / 39.0 | 54.2 / 75.5 | 85.0 |
| tx | 24.9 / 36.0 | 28.1 / 36.2 | 47.0 / 73.8 | 79.4 |

## Table M7. Training set against RefSeq: old pipeline vs new selection (no filterGeneMark keeper set; upper bound, see text)

| Genome | Old: models, exact % | New: models, exact % | Exact RefSeq genes old → new | Redundant old → new |
|---|---|---|---|---|
| N. crassa OR74A | 2889, 38.9 | 1515, 78.7 | 1124 → 1192 | 46 → 0 |
| A. nidulans FGSC A4 | 5007, 50.1 | 4123, 62.1 | 2508 → 2559 | 37 → 4 |
| B. cinerea B05.10 | 5637, 75.5 | 4866, 89.4 | 4255 → 4348 | 59 → 0 |

## Table M8. getBestModel ranking against RefSeq (PASA models passed to EVM): exact CDS / exact intron chain

| PASA source | Genome | complete-first | structure-first | guarded |
|---|---|---|---|---|
| rc1 | N. crassa OR74A | 1549 / 2060 | 1530 / 2225 | 1536 / 2183 |
| rc1 | A. nidulans FGSC A4 | 3224 / 3210 | 3202 / 3221 | 3222 / 3221 |
| rc1 | B. cinerea B05.10 | 5455 / 5030 | 5419 / 5147 | 5447 / 5113 |
| R1 | N. crassa OR74A | 2207 / 3020 | 2174 / 3255 | 2187 / 3206 |
| R1 | A. nidulans FGSC A4 | 3258 / 3234 | 3238 / 3248 | 3258 / 3249 |
| R1 | B. cinerea B05.10 | 5767 / 5398 | 5726 / 5515 | 5761 / 5475 |
| R1R2 | N. crassa OR74A | 2303 / 3060 | 2275 / 3259 | 2290 / 3229 |
| R1R2 | A. nidulans FGSC A4 | 3259 / 3235 | 3239 / 3248 | 3259 / 3249 |
