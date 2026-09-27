---
topic: ani-method-suitability
description: Suitability of ANI tools (skani / fastANI / mash / sourmash / LZ-ANI) for genome-level all-vs-all ANI on fungal-scale genomes, incl. containerized evaluation results.
created: 2026-08-18
last_updated: 2026-08-18
status: active
---

# ANI Method Suitability

> Opened 2026-08-18 while evaluating LZ-ANI 1.2.3 as a fifth `--ani_method` for
> the `compare_ANI` workflow. Tracks tool-level fitness (semantics of the output,
> memory/thread scaling, failure modes, output-format contracts) for the identity
> methods in the ANI pipeline.

## Findings

### F-009 — LZ-ANI 1.2.3 is unsuitable as a genome-level fungal ANI method (silent no-output + OOM; scaffold-level only in its working mode)

**Status:** established (2026-08-18, 4-core / 503 GB host, real fungal scaffolds)

**Claim:** In an evaluation driven by custom comparison tooling, the containerized
LZ-ANI 1.2.3 (`docker://quay.io/biocontainers/lz-ani:1.2.3--h9ee0642_0`) could not
produce genome-level ANI for real fungal genomes:

- **Correct genome-level mode (`--multisample-fasta false`):** silently writes **no
  output** — even on single-genome runs (exit 0, normal logs, ~4.9 GB RSS).
- **Default mode (`--multisample-fasta true`):** emits a clean 3-column table but
  **scaffold-level** (each contig treated as a sample; `NW_022983500.1`-style IDs),
  so it is not comparable to genome-level skani/fastANI/mash ANI.
- **Memory/thread scaling:** genome-level runs OOM-kill (`exit 137`) at `-t >= 4`
  even for 1–2 genomes; peak RSS ~5 GB per 1–2 fungal genomes (37.8/32.2/52.3 MB
  scaffolds).
- **Output contract:** `--out-format` is ignored on small inputs — emits a sectioned
  `[lz_similarities]` file indexed by sample number instead of a table.

**Why it matters:** silent-no-output on real data is the failure class that cannot be
caught downstream (no missing TSV to notice when a tool exits 0). LZ-ANI was
de-scoped: compare_ANI now supports `skani | mash | sourmash | fastani` only. The
method-agnostic comparison tooling (`scripts/compare_ani_methods.py`) is unaffected
and remains useful for cross-method agreement/timing (e.g. `--methods skani,fastani`).

**Evidence:** `analysis/ani_method_evaluation/LZANI_PERFORMANCE.md` (measurements,
thread/CPU/memory probes, reproduction recipe, both container modes); test artifacts
in `/scratch/jstajich/27545314/opencode/lzani_real/` (genome_list.txt, out_2.tsv,
s1/s2 synthetic probes).

**Not yet done:** a systematic comparison of the *remaining* four methods'
agreement on a shared clade (only skani production results + a fastANI prefilter
design exist); the LZ-ANI evaluation was tooling-motivated, not a full method
benchmark.
