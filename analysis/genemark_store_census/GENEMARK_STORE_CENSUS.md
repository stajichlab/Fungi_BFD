# genemark_store census (2026-11-25)

Script: `scripts/one-off/genemark_store_census.py --runs Fungi_BFD_runs`. Data: `census_20261125.tsv`, log: `census_20261125.log`.
Sizes are stat() apparent bytes (GB = 1e9). Read-only.

| store | species dirs | .genemark.gtf GB | .genemark.mod GB | .other.gff3 GB |
|---|---:|---:|---:|---:|
| do_annotation_fusarium_luffae/genemark_store (shared) | 23,831 | 262.7 | 22.4 | 0.07 |
| do_annotation/genemark_store | 1,126 | 11.1 | 1.0 | 0.02 |
| do_annotation_big_drops_step2 | 81 | 0.6 | 0.1 | 0.00 |
| do_annotation_busco_gap_earlgrey | 51 | 0.8 | 0.0 | 0.00 |
| others (backup_pre_wave3..., busco_gap, rosellinia_step2) | 10 | <0.1 | <0.1 | 0 |
| TOTAL | 25,099 | 275.3 | 23.5 | 0.09 |

- `Fungi_BFD_runs/genemark_store` and `do_annotation_{wave1,wave2,wave3,fusarium_luffae}/genemark_store` resolve to one real directory (counted once).
- 24,890 plain GTFs, 0 `.gtf.gz`, 24,100 `.mod`. 209 shared-store folders have no plain GTF.
- The GTF is 92% of the bytes. The `.other.gff3` files are almost all 0 bytes.
- Estimated saving: 275.3 GB / about 10.8 (gzip ratio from a 12-file sample, 91.3 MB -> 8.5 MB; 3 more real folders gave 10.7x; 1 genome gave 11.7x) = about 25 GB left, about 250 GB freed. This is an estimate from 16 files, not a measurement.
