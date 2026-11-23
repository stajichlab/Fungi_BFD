# Ubiquitous Language

Domain glossary for the Fungi_BFD RNA-seq acquisition and funannotate annotation pipeline. Terms here are how the team should refer to concepts in code, commits, docs, and conversation.

## RNA-seq acquisition

| Term | Definition | Aliases to avoid |
| ---- | ---------- | ---------------- |
| **SRA query** | The cached per-species CSV of candidate RNA-seq runs (`rnaseq_reads/sra_query/<species_tag>.sra_query.csv`) produced by `SRA_QUERY`/`SRA_QUERY_BATCH` | query cache, sra csv |
| **Accession** | A single SRA run identifier (e.g. `SRR2093556`) listed in the SRA query | run, SRR, sample |
| **Layout** | The SRA-declared read structure of an accession — `PAIRED` or `SINGLE` (column 6 of the SRA query) | pairing, type |
| **Normalized reads** | The final per-species read files consumed downstream: `<species_tag>_norm_R1/R2/SE.fastq.gz` after truncation, bbnorm, and fastp | processed reads, cleaned reads |
| **PE fetch** | The `SRA_FETCH` process: downloads paired-end accessions and emits `_norm_R1` + `_norm_R2` (with a zero-byte `_norm_SE` stub) | SRA_FETCH, paired fetch |
| **SE fetch** | The `SRA_FETCH_SE` process: downloads single-end accessions and emits `_norm_SE` (with zero-byte `_norm_R1`/`_norm_R2` stubs) | SRA_FETCH_SE, single fetch |
| **Trinity assembly** | Trinity genome-guided transcriptome (`rnaseq_data/<species_tag>.trinity-GG.fasta`) built from normalized reads. Not the genome assembly (`ASMID`, `asm_stats`) | assembly (alone), trinity, transcriptome |

## Read-file states

| Term | Definition | Aliases to avoid |
| ---- | ---------- | ---------------- |
| **Real reads** | A `_norm_*.fastq.gz` file with actual sequence content (non-zero bytes) | data file, populated |
| **Stub** | A deliberately zero-byte `_norm_*.fastq.gz` written so every species carries the same 4-tuple output shape (present but empty) | placeholder, empty file, dummy |
| **Missing output** | A declared output file absent from `rnaseq_reads/` entirely — distinct from a stub, which exists but is empty | missing file, gap |
| **Complete store** | A species whose `rnaseq_reads/` holds *all three* declared outputs (R1, R2, SE), whether real or stub — the condition for `storeDir` to skip re-execution | cached, done |

## Caching & re-execution

| Term | Definition | Aliases to avoid |
| ---- | ---------- | ---------------- |
| **storeDir** | Nextflow directive pointing a process's outputs at a persistent directory; the process is skipped only when *every* declared output already exists there | cache dir, output dir |
| **Cache miss** | Nextflow deciding to (re)run a task; here caused by an *incomplete store*, not a changed input | re-trigger, invalidation |
| **Incomplete store** | A `storeDir` folder missing one or more declared outputs for a species, forcing re-execution even though other outputs are present and valid | partial cache, stale cache |
| **Schema drift** | A mismatch between what a process now declares (e.g. the added `_norm_SE` output, or the 6-column SRA query) and what older on-disk artifacts contain | version skew, format change |
| **Backfill** | Retroactively writing the missing artifact (e.g. a zero-byte `_norm_SE` stub) so an incomplete store becomes complete and stops re-triggering | patch, fill-in |

## Blacklist & routing

| Term | Definition | Aliases to avoid |
| ---- | ---------- | ---------------- |
| **Blacklist** | `rnaseq_blacklist.csv` — per-accession override table (column 4 = action) consulted during fetch | skiplist, exclusions |
| **Action** | The per-accession override in the blacklist: `skip`, `rename_headers`, or `SE_trinity` | flag, override, strategy |
| **SE_trinity** | An action marking an accession that SRA labels `PAIRED` but is actually single-end, routing it to SE fetch | mislabeled-SE, se flag |
| **SE candidate** | An accession detected at fetch time as single-end (yielded `_1` but no `_2`); logged to `rnaseq_se_candidates.csv` for review before adding to the blacklist | se guess, unpaired candidate |
| **Blacklist candidate** | An accession whose download failed outright (both pfd and EBI fallback); logged to `rnaseq_blacklist_candidates.csv` for review | failed download, dead accession |
| **EBI fallback** | The aria2 download of `_1`/`_2` FASTQs directly from EBI FTP when `parallel-fastq-dump` produces mismatched read names | ftp fallback, ena fallback |

## Relationships

- A **species_tag** has exactly one **SRA query**, which lists zero or more **accessions**.
- Each **accession** has a **layout** and may carry a blacklist **action**.
- **PE fetch** and **SE fetch** each emit the same 4-tuple (`species_tag`, R1, R2, SE); the process that doesn't own a given orientation writes **stubs** for it.
- A species is skipped by **storeDir** only when its store is **complete** — all three `_norm_*` files exist (real or stub).
- A **missing output** (not a stub) is what produces the **incomplete store** → **cache miss** → re-run.
- **Schema drift** (the `_norm_SE` output added after some reads were built) is the source of the current **incomplete stores**; **backfill** resolves them.

## Example dialogue

> **Dev:** "Pseudocercospora_eumusae has real R1 and R2 in `rnaseq_reads/`, so why is **PE fetch** re-running it instead of skipping?"

> **Domain expert:** "Because **storeDir** only skips on a **complete store**. The process declares three outputs now — R1, R2, and `_norm_SE` — and that species has no SE file at all. That's a **missing output**, so it's an **incomplete store** and Nextflow calls a **cache miss**."

> **Dev:** "But those reads were fine. Nothing about the input changed."

> **Domain expert:** "Right — this is **schema drift**, not stale data. The `_norm_SE` output was added after the reads were built. The fix is a **backfill**: write a zero-byte SE **stub** so the store reads as complete."

> **Dev:** "And when it does re-run, it dies. The accession is single-end but going through **PE fetch**?"

> **Domain expert:** "Yes. Its **SRA query** is the old 5-column format with no **layout**, and it has no `SE_trinity` **action**, so **SE fetch** never claims it. It falls into **PE fetch**, `--split-files` yields only `_1`, it gets logged as an **SE candidate**, the **EBI fallback** 404s, and the paired enforcement step exits with no reads to work on."

## Flagged ambiguities

- **"cached" / "it should just use the cache"** — conflates two different conditions. Nextflow's task cache (work-dir hashing) is not in play here; `storeDir` is, and it keys purely on *presence of declared outputs*, not input hashes. Say **complete store** vs **incomplete store**, not "cached" vs "not cached."
- **"empty" / "missing"** — a **stub** (zero-byte file that exists) and a **missing output** (no file) behave oppositely for `storeDir`: a stub *satisfies* the store, a missing file *breaks* it. Keep these distinct; don't call a stub "missing" or a missing file "empty."
- **"SE"** was used for three related-but-distinct things: an accession's **layout** (`SINGLE`), the **SE_trinity** blacklist action for mislabeled-paired data, and the **SE fetch** process. Qualify which one is meant.
- **"placeholder"** — used loosely for both the deliberate output-shape **stub** and the to-be-written **backfill** file. They're the same artifact (zero-byte SE), but name the intent: a **stub** is what a fetch process writes for the orientation it doesn't own; a **backfill** is one written after the fact to repair an incomplete store.
