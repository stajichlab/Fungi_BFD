# PASA + MariaDB troubleshooting

`FUNANNOTATE_TRAIN`/`FUNANNOTATE_UPDATE` run PASA with `--PASACONF` pointed at a
per-task MariaDB instance started inside the funannotate container itself
(baked-in `mariadbd`, no separate sidecar container as of 2026-09-09 — see
`nextflow/modules/funannotate/predict/FUNANNOTATE_TRAIN/main.nf`). MariaDB has
been a recurring source of pipeline-wide failures because the assumptions
PASA's schema and Debian's default configs make about MySQL/MariaDB are
decades out of date. This doc collects the known failure modes and fixes so
they don't get rediscovered from scratch each time. The config file all of
this applies to is `nextflow/assets/pasa_conf/my.cnf`.

## `character-set-server`/`collation-server` must be `latin1` (added 2026-09-09)

**Symptom:**
```
DBD::mysql::st execute failed: Can't drop database '<tag>_pasa'; database doesn't exist ...
ERROR 1071 (42000) at line 151 in file: '.../schema/cdna_alignment_mysqlschema':
Specified key was too long; max key length is 3072 bytes
CMD: .../bin/mysql -ufunannotate -pPASAfun ... -e 'source .../cdna_alignment_mysqlschema' failed.
```
`funannotate train` exits 1 during `Running PASA alignment step`, after
MariaDB itself has already started successfully (socket up, "ready for
connections" logged) — so this is NOT a repeat of the container-startup
issues below, it's a distinct, later-stage failure.

**Cause:** PASA's bundled schema (`cdna_alignment_mysqlschema`, part of the
PASA install, not something this repo controls) was written assuming MySQL
5.x's old default charset, `latin1` (1 byte/char). This MariaDB 11.8 build
defaults to `utf8mb4` (4 bytes/char) when no server-level charset is
configured. Every `VARCHAR(N)` index in the schema becomes 4x wider under
`utf8mb4` than it was sized for, and at least one composite/VARCHAR key
exceeds InnoDB's 3072-byte max key length as a result. This is a well-known
class of issue for old MySQL-era tools (like PASA) run against modern
MySQL/MariaDB installs.

**Fix:** force the server back to the byte width PASA's schema expects, in
`nextflow/assets/pasa_conf/my.cnf`'s `[mysqld]` section:
```
character-set-server	= latin1
collation-server	= latin1_swedish_ci
```
Confirmed against the failed run at
`do_annotation_asco/work/funannotate/05/320769eda949fca77e3d742eaabae9`
(`Alternaria_alternata_DZ`, 2026-09-09).

## Other MariaDB gotchas already fixed in `my.cnf` (context, not new)

These are documented inline as comments in `my.cnf` itself; summarized here
so they're findable from one place:

- **No `user = ...` directive.** mariadbd runs inside the funannotate
  container as whichever OS user submitted the Nextflow job, never as root,
  so there's no privileged-to-unprivileged user to switch to. A hardcoded
  `user = jstajich` broke the pipeline for every other user (mariadbd tries
  to setuid to an account that isn't the one that launched it).
- **`bind-address = 127.0.0.1` only.** This instance exists only for the
  lifetime of one `FUNANNOTATE_TRAIN`/`UPDATE` task, talking only to the PASA
  process launched by the same task in the same network namespace. Binding
  wider than loopback was confirmed exposing mysqld to the whole cluster
  network fabric with no corresponding need.
  - `!include /etc/mysql/mariadb.cnf` / `!includedir /etc/mysql/conf.d/` are
  commented out. Vestigial from a stock Debian/Ubuntu system layout that
  never exists in this per-task scratch/conda context; newer
  `my_print_defaults` treats a missing `includedir` target as **fatal**
  ("Stopped processing the 'includedir' directive ... Program aborted!")
  rather than silently skipping it, which broke every
  `pasa_mysql=true` task after a mysql-libs/mysql-common version bump.

## Earlier container-startup issue (superseded design, kept for context)

`nextflow/docs/mariadb-10.3.9.def` (an Apptainer/Singularity definition for a
separate MariaDB *sidecar* container, since replaced by mariadb baked
directly into the funannotate container) documents that
`apptainer instance start <sif> <name> [args...]` only runs `[args...]` if
the image defines a `%startscript` — Docker's `ENTRYPOINT`/`CMD` only
populate `%runscript` (used by `apptainer run`), which `instance start`
never consults. Without a `%startscript`, `instance start` silently starts an
empty instance and never launches `mysqld_safe`, and PASA fails with
`Can't connect to MySQL server ... (111)`. Not applicable to the current
baked-in design, but keep in mind if sidecar-container MariaDB is ever
revisited.

## If a new MariaDB/PASA failure shows up

1. Check whether the server actually started (`ready for connections` in
   `.command.log`/`.command.err`) before assuming this is a repeat of the
   startup issues above — a clean startup followed by a schema/query error
   (like the key-length issue) is a different class of problem.
2. Check `<workdir>/train_local/training/pasa/pasa-assembly.log` for the
   actual PASA/mysql error, not just funannotate's own summary line — the
   funannotate log truncates the underlying `mysql`/`mysqld` error.
3. Add the fix (and the failure signature) to this file.
