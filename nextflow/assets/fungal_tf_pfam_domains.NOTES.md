# fungal_tf_pfam_domains.csv: notes

Columns: `pfam_acc`, `pfam_name` (Pfam short name, `#=GF ID`), `tf_family` (label).
The CSV cannot hold comments. Notes are here.
Check the list with `nextflow/bin/check_tf_pfam_list.py`.

## Release checked

Pfam 38.2. File: `/srv/projects/db/pfam/current/Pfam-A.hmm.dat`.

## Change: fix wrong accessions

| Old row | Pfam 38.2 name of that accession | Action |
|---|---|---|
| PF00446 `SBP` | GnRH | Replaced by PF03110 `SBP` (SBP zinc finger). It has 0 hits in BFD proteomes (plant-type domain). |
| PF11955 `WOPR` | PORR | Removed. |
| PF04844 `DNA_binding_1` | Ovate (plant family) | Removed. Replaced by PF04383 `KilA-N`, the APSES DNA-binding domain. |
| PF11951 `zf-C6H2` | Fungal_trans_2 | Kept. Name and family label corrected. The hits are real. |
| PF08880 `Rrn7` | QLQ (SWI/SNF-related domain) | Removed. Not a TF DNA-binding domain. |

Pfam 38.2 has no WOPR model.
WOPR proteins (for example Wor1, Ryp1) are not counted now.
A custom HMM is needed to count them.

## Effect on outputs (BFD counts as of Pfam 38.2 hits)

- `tf_inventory`: PF08880 rows disappear (19,290 proteins in 18,376 genomes).
- `tf_inventory`: PF04383 rows appear (about 75,540 hits in 23,033 genomes).
- `tf_inventory`: PF11951 rows stay. Their `pfam_name` and `tf_family` change.
- `tf_inventory`: PF00446, PF11955, PF04844 had no hits. PF03110 has 0 hits. No row change.

## Rebuild

Existing tables do not change until MERGE_TF runs again.
`merge_if_stale` compares file mtimes only. The CSV is not an input of the task.
Run MERGE_TF with `--force_merges MERGE_TF`.
Then rebuild tables and views made from `tf_inventory` (for example `tf_genome`, `v_tf_family_counts`)
and reload the DuckDB views (`nextflow/bin/build_BFD_duckDB.sh`).

## Added families used in the BFD TF survey

Rows added for the TF survey (BFD_dashboard `config/tf_family_map.tsv`).
The `tf_family` text is the `family` column of that file.
`pfam_name` is the Pfam 38.2 ID from the dat file.

| pfam_acc | pfam_name | tf_family |
|---|---|---|
| PF02319 | WHD_E2F_TDP | E2F_TDP |
| PF01167 | Tub | Tub |
| PF02257 | RFX_DNA_binding | RFX_DNA_binding |
| PF03106 | WRKY | WRKY |
| PF04516 | CP2 | CP2 |
| PF00808 | CBFD_NFYB_HMF | CBFD_NFYB_HMF |
| PF00649 | Copper-fist | Copper-fist |
| PF01285 | TEA | TEA |
| PF01017 | STAT_alpha | STAT |
| PF02864 | STAT_bind | STAT |
| PF21354 | STAT_linker | STAT |

Adding rows changes the size of `tf_inventory`.
The TF survey builder already counted these families from the full `pfam` table.
`tf_inventory` in BFD does not have them yet.
PF00808 is a histone-fold domain. It also occurs in proteins that are not CBF/NF-Y TFs.
Its counts are an upper bound (23,334 genomes in the survey).
This block is a separate commit. Drop that commit to keep only the accession fixes.
