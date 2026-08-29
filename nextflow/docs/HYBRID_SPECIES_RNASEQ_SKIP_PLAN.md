# Plan: composite parent-transcript RNA-seq evidence for hybrid-cross species

Status: IMPLEMENTED 2026-08-28 (reviewed by an independent Opus pass, revised
per maintainer direction toward the composite-evidence strategy, refined via
a grilling session that resolved the remaining parameters, then implemented
and validated via `nextflow -preview` + a partial `-stub-run`). The
`samples.csv HYBRID` column mentioned in earlier sections below was added
and then reverted the same day (hybrids are <0.5% of rows; `hybrid_parentage.csv`
membership alone is the runtime signal — see "Detection" below) — where a
section still shows the column, it's describing that superseded intermediate
state, not current behavior. Bounded change (the flow it touches —
`FUNANNOTATE_RNASEQ.nf` — already exists; see
`DATA_FLOW_ANI_REPRESENTATIVE_RNASEQ_TRAINING.md` for the full context this
plan builds on).

**Current recommendation, final (all open parameters resolved by grilling,
2026-08-28):**
1. RNA-seq evidence: build a **composite transcript FASTA** per hybrid
   species from its parent species' *already-built* Trinity assemblies
   (free — no new SRA queries, no new Trinity runs), feed it to
   `FUNANNOTATE_TRAIN` via `--trinity` with no reads, under a new relaxed
   `composite` PASA tier (provisional thresholds, shipped now, calibrated
   later). True skip (ab-initio-only) is the fallback only when no
   parent/genus-mate Trinity exists at all. **On by default** once
   implemented — no master flag, just an after-the-fact audit report.
2. Ab-initio params: **no code change needed.** `pick_representative_strain.py`
   already groups strains strictly by the exact `SPECIES` string
   (`by_species = defaultdict(list)`, keyed on the raw field,
   `pick_representative_strain.py:198-202`), so a hybrid never pools with its
   parent species or with a different hybrid cross today — reuse already
   happens only within an identical hybrid-cross group, which is exactly the
   scope the maintainer wants. (An earlier draft of this plan proposed
   excluding hybrids from reuse entirely, i.e. even within their own cross
   group — that was a misreading of the maintainer's intent, corrected below.)

## Problem

`samples.csv` currently has 90 rows across 13 distinct interspecific/hybrid
"species" values (any row whose `SPECIES` field contains a standalone `x`
token, e.g. `Saccharomyces cerevisiae x Saccharomyces eubayanus`, or the
formal nothospecies form `Saccharomyces x bayanus`):

```
37  Saccharomyces cerevisiae x Saccharomyces eubayanus
20  Saccharomyces eubayanus x Saccharomyces uvarum
 9  Saccharomyces cerevisiae x Saccharomyces kudriavzevii
 7  Saccharomyces x bayanus
 5  Cryptococcus neoformans x Cryptococcus deneoformans
 3  Saccharomyces cerevisiae x S. eubayanus x S. kudriavzevii x S. uvarum
 2  Saccharomyces cerevisiae x Saccharomyces kudriavzevii x Saccharomyces uvarum
 2  Saccharomyces cerevisiae x Saccharomyces eubayanus x Saccharomyces uvarum
 1  Saccharomycopsis fibuligera x Saccharomycopsis cf. fibuligera
 1  Saccharomyces cerevisiae x Saccharomyces uvarum
 1  Saccharomyces cerevisiae x Saccharomyces kudriavzevii x Saccharomyces eubayanus
 1  Ogataea polymorpha x Ogataea parapolymorpha
 1  Fusarium meridionale x Fusarium asiaticum
```

Today's pipeline (Stage 4 of the data-flow doc) treats each of these 13
values as an ordinary species: it picks one strain as the RNA-seq
representative, builds one shared `rnaseq_data/<species_tag>.trinity-GG.fasta`
from that strain's own reads/genome, and fans that single assembly out to
*every* sibling strain's `FUNANNOTATE_TRAIN` (e.g. all 37
`Saccharomyces_cerevisiae_x_Saccharomyces_eubayanus` strains PASA-align
against one shared Trinity).

**Why that doesn't make biological sense here.** Interspecific hybrids
(especially the well-studied *Saccharomyces* lager/cider hybrids) are mosaic:
each isolate can carry a different chromosome-copy-number balance between
parental subgenomes, different introgression breakpoints, and even different
ploidy — two isolates sharing the same nominal "species" label (same pair of
parents) are not interchangeable the way two strains of a true biological
species are. Sharing one representative's Trinity-GG assembly across the
group assumes strain-to-strain transcriptome similarity that the biology
doesn't support.

**Correction (external review, 2026-08-28):** the originally-cited empirical
failure does **not** actually demonstrate this. `Saccharomyces_x_bayanus_FM677`
(`work/funannotate/dc/6f90730bd25a9653c1da5760003e86/`) — genome-guided
Trinity aligning only ~11 reads across 6 clusters, 0 transcripts — is almost
certainly a **read-acquisition failure**: the per-species SRA query for this
narrow nothospecies taxid simply returned next to nothing, not a case of
reads existing but mis-anchoring to a mosaic genome. The `misc/poor_trinity/`
artifacts corroborate this reading — `Saccharomyces_x_bayanus.trinity-GG.fasta`
(252 B) and `Saccharomyces_eubayanus_x_Saccharomyces_uvarum.trinity-GG.fasta`
(54 KB) are near-empty inputs, not large-but-mis-anchored ones. The
mosaic-genome argument is still a valid reason not to *share* one
representative's Trinity across differently-recombined siblings, but it is a
distinct claim from "hybrid RNA-seq doesn't align" — this plan rests on two
separate justifications (sparse/unreliable per-hybrid-taxid SRA data, and
sharing-across-mosaic-siblings being unsound), not one confirmed failure.

Every sibling strain still pays a real `FUNANNOTATE_TRAIN` job (Stage 4 fans
out to every strain regardless of `--predict_scope` — see
`HOW_SIBLINGS_ARE_TRANSCRIPT_TRAINED.md`) even when the shared evidence is
this thin, so the cost is not confined to the one representative that builds
the Trinity — for the 37-strain `S. cerevisiae x S. eubayanus` group, a
useless representative assembly means up to 37 wasted PASA-alignment jobs on
top of it.

## RNA-seq evidence strategy: composite parent-transcript FASTA

### Why this beats a blanket skip

Every parent species of every composite-cross hybrid in this dataset already
has its own real, non-hybrid Trinity assembly, built from that parent's own
samples.csv rows and RNA-seq:

```
rnaseq_data/Saccharomyces_cerevisiae.trinity-GG.fasta      (exists)
rnaseq_data/Saccharomyces_eubayanus.trinity-GG.fasta       (exists)
rnaseq_data/Saccharomyces_kudriavzevii.trinity-GG.fasta    (exists)
rnaseq_data/Saccharomyces_uvarum.trinity-GG.fasta          (exists)
rnaseq_data/Cryptococcus_neoformans.trinity-GG.fasta       (exists)
rnaseq_data/Cryptococcus_deneoformans.trinity-GG.fasta     (exists)
rnaseq_data/Fusarium_meridionale.trinity-GG.fasta          (exists)
rnaseq_data/Fusarium_asiaticum.trinity-GG.fasta            (exists)
rnaseq_data/Ogataea_polymorpha.trinity-GG.fasta            (exists)
rnaseq_data/Ogataea_parapolymorpha.trinity-GG.fasta        (exists)
```

Building a composite transcript-evidence set per hybrid species is therefore
essentially free — no new SRA queries, no new Trinity assembly, just
concatenating files the pipeline has already built. It's also arguably more
scientifically sound than the original per-hybrid-representative design ever
was: PASA aligns transcripts locus-by-locus against the target genome, so a
composite of the *real parent transcriptomes* should align each hybrid
subgenome block to whichever parent it actually descends from — better
grounded than sharing one hybrid isolate's own admixed Trinity assembly
across differently-recombined siblings (the mosaic-genome objection above,
applied constructively instead of as a reason to skip).

### Parent-name resolution — `hybrid_parentage.csv`, single runtime lookup

**Resolved (grilling, 2026-08-28), generalized beyond the original
nothospecies-only lookup table.** Rather than three live-runtime paths
(string-split, else lookup table, else genus fallback), commit **one**
narrow/long CSV covering every hybrid `species_tag` in the dataset, and make
it the sole runtime source of truth:

```csv
hybrid_species_tag,parent_species
Saccharomyces_cerevisiae_x_Saccharomyces_eubayanus,Saccharomyces cerevisiae
Saccharomyces_cerevisiae_x_Saccharomyces_eubayanus,Saccharomyces eubayanus
Saccharomyces_x_bayanus,Saccharomyces cerevisiae
Saccharomyces_x_bayanus,Saccharomyces eubayanus
...
```

One row per (hybrid, parent) pair — handles 2-way through N-way crosses
uniformly (the 4-way `S. cerevisiae x S. eubayanus x S. kudriavzevii x
S. uvarum` cross is just 4 rows, no ragged/wide columns, no fixed max-arity).
CSV, not YAML — matches every other lookup/override file already used in
this pipeline (`rnaseq_representative_override.csv`, `rnaseq_blacklist.csv`,
`abinitio_reuse_assignments.csv`), loaded with the same `readLines()`/
`split(',')` idiom already in `utils.nf`, no new parser dependency.

**Population**: `scripts/bootstrap_hybrid_metadata.py` scans `samples.csv`'s
`SPECIES` column (read-only — no persisted marker column; see "Detection"
below) and splits each hybrid `SPECIES` value on `x`/`×`, writing a row per
resolved parent for the 12 string-decomposable combos (verified directly
against the current dataset — `SPECIES` splits cleanly into full parent
binomials for 12 of the 13 hybrid combos). The one exception, the formal
nothospecies shorthand `Saccharomyces x bayanus` (`bayanus` alone isn't a
resolvable second parent token from the string) — literature-documented as
the lager-yeast hybrid lineage (historically `S. pastorianus`/
`S. carlsbergensis`), descending from *S. cerevisiae × S. eubayanus* — is
looked up from a small hardcoded table in that same script (one entry
today), sourced from that outside knowledge rather than inferred. Re-running
the script is idempotent and `--no-clobber` by default (existing rows for a
species already present are left untouched, so a manually-corrected entry
never gets silently overwritten); `--lint` reports drift without writing.
After bootstrap, the string-split logic isn't part of the runtime path at
all — it only runs inside that script.

Handling at runtime (`FUNANNOTATE_RNASEQ.nf`), in order:
1. `hybrid_species_tag` has ≥1 row in `hybrid_parentage.csv` → use those
   parent species' Trinity files directly.
2. No rows (a hybrid added to `samples.csv` without a matching
   `hybrid_parentage.csv` update — `scripts/bootstrap_hybrid_metadata.py
   --lint` catches this, nonzero exit if any hybrid species has no coverage)
   → **genus-wide fallback**: concatenate every non-hybrid
   `rnaseq_data/<Genus>_*.trinity-GG.fasta` under the hybrid's `GENUS` field.
3. Fallback produces nothing either (no genus-mates with Trinity) → true
   skip (ab-initio-only), the last-resort path.

### Building the composite — real channel dependency, not a file-existence check

Must not become a filesystem `.exists()` probe at script-generation time —
that risks a race if a hybrid's composite-build step is scheduled before its
parent species' own `RNASEQ_PREPARE` completes in the same run (they're
independent branches of the same DAG with no inherent ordering guarantee
otherwise). Instead: key off `RNASEQ_PREPARE.out.shared` (already emits
`tuple(species_tag, path(trinity_fa))` for every species processed this run)
collected into an offline map (mirrors how `repr_ch` already does an
offline-map lookup against `abinitioReuseMap`), so Nextflow's own DAG
scheduling — not convention — guarantees parents are built first when they're
part of the same run, and correctly falls back to reading an already-on-disk
parent Trinity from a *prior* run otherwise (`file()` check only as the
not-in-this-run fallback, same pattern `RNASEQ_PREPARE` itself already uses
for pre-existing training output).

New lightweight process, e.g. `BUILD_HYBRID_COMPOSITE_TRINITY`: input
`(hybrid_species_tag, list_of_parent_trinity_paths)`; output
`rnaseq_data/<hybrid_species_tag>.composite-parents.trinity-GG.fasta`
(simple `cat`). `storeDir`'d the same way `RNASEQ_PREPARE` is, keyed on the
hybrid `species_tag` — inherits the same storeDir-invalidation gap noted
under "Cleanup required" below, so the same quarantine/invalidation caveat
applies here too once `DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md`'s Option 2
sentinel exists.

### Feeding it to `FUNANNOTATE_TRAIN` — new readless branch

Existing `FUNANNOTATE_TRAIN/main.nf` has a "PASA only, no reads" branch, but
it still references `R1_REAL`/`R2_REAL` in its `--left_norm`/`--right_norm`
args (a latent bug — worth a separate one-line fix regardless of this plan:
those should be omitted entirely when there are no reads, not passed as
empty strings). The composite-evidence case needs a genuinely new branch:
`--trinity <composite-parents.trinity-GG.fasta>` with **no**
`--left_norm`/`--right_norm`/`--single_norm` at all.

### Detection — `hybrid_parentage.csv` membership, no samples.csv column

Original design considered adding a `HYBRID` boolean column to `samples.csv`
(extra taxonomy columns are already accepted-but-optional per
`assets/schema_input.json`) so hybrid status would be reviewable, committed
data rather than a live runtime regex. **Reverted** (2026-08-28): interspecific
hybrids are <0.5% of `samples.csv` rows (90/23,683) — not worth a whole extra
column on every row for that few. Since `hybrid_parentage.csv`'s row
membership already had to exist as authoritative, reviewable data for the
parent-resolution mechanism above, it can *also* serve as the hybrid-status
signal directly — no second, redundant place to keep in sync. Final design:

- `loadHybridParentage()` (`utils.nf`) is the single source of truth:
  `hybridParentage.containsKey(species_tag)` — presence as a key IS "this
  species is a hybrid." No column, no per-row flag anywhere in `samples.csv`.
- `scripts/bootstrap_hybrid_metadata.py` still uses a regex
  (`/(?:\sx\s|×)/i`) internally, but only to *scan* `samples.csv` (read-only)
  when generating/re-linting `hybrid_parentage.csv` — the regex never runs at
  pipeline execution time, and `samples.csv` itself is never written by it.
  `--lint` mode reports (nonzero exit) any hybrid species with zero
  `hybrid_parentage.csv` coverage, without writing anything — the drift
  check the `HYBRID`-column design was trying to get, without the column.

### Where it plugs in — a branch ahead of SRA acquisition, not inside `rnaseqPolicyFor`

Original design proposed folding the hybrid rule into
`DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` Option 3's proposed
`rnaseqPolicyFor(out)` generalization of `pasaTierFor`, reasoning that one
policy function beats two independent decision mechanisms. **Superseded by
implementation reality**: a hybrid `species_tag` must never even reach the
normal RNA-seq representative pick (`repr_ch`) or trigger `RNASEQ_PREPARE`
in the first place (see below) — that decision happens *before* any
per-strain tiering logic like `pasaTierFor`/`rnaseqPolicyFor` would ever run,
so folding hybrid handling into that function wasn't actually possible
without restructuring it to also gate the representative-pick step, which
would have made it do two unrelated jobs. What's actually implemented:

```groovy
// FUNANNOTATE_RNASEQ.nf, top of the run_sra_fetch block:
def genome_branched = predict_genome_ch.branch { out, asmid, species, strain, locustag, busco, hlen, ttable, genome_fa, taxonid ->
    hybrid: hybridParentage.containsKey(species.replaceAll(/\s+/, '_'))
    normal: true
}
def normalGenomeCh = genome_branched.normal   // unchanged existing pipeline
def hybridGenomeCh  = genome_branched.hybrid  // composite-evidence path below
```

`hybridGenomeCh` never enters `sra_input`/`SRA_QUERY_BATCH`/`repr_ch`/
`RNASEQ_PREPARE` at all. Its own composite-build path (parent lookup →
`BUILD_HYBRID_COMPOSITE_TRINITY` → empty-reads placeholders → per-strain
`(out, ..., trinity_fa, pasa_tier)` tuples, `pasa_tier` computed inline as
`'composite'`/`'skip'` by composite success, not by `pasaTierFor`) produces
the *same tuple shape* `train_input` already uses, and is `mix()`ed into
`train_input` right before the existing `ani_skip`/`has_rnaseq`/`no_rnaseq`
`.branch{}` call — so `pasaTierFor` itself is completely untouched, and every
downstream staleness-filtering/`FUNANNOTATE_TRAIN`-invocation/caching rule
already in place applies to hybrids with zero duplicated logic. One audit
trail (the `pasa_tier` value each row already carries into `FUNANNOTATE_TRAIN`),
just produced by two different code paths feeding the same channel, not by
two competing policy functions.

This still needs a companion fix to `repr_ch` (`FUNANNOTATE_RNASEQ.nf:217-233`):
a hybrid `species_tag` must never trigger a normal `RNASEQ_PREPARE` build
(picking one hybrid strain as "the representative" and running
`funannotate train` on its own reads/genome) — implemented above via the
`genome_branched` split: hybrid `species_tag`s never enter `repr_ch` at all,
routing to `BUILD_HYBRID_COMPOSITE_TRINITY` instead.

### New PASA tier: `composite`

Cross-species identity within a genus is looser than the existing `relaxed`
tier's thresholds (85% avg-id / 70% aligned, tuned for 90-97% ANI intra-
species divergence). *Saccharomyces* sensu stricto members run roughly
80-90% nucleotide identity to each other — a `composite` tier needs its own,
more permissive `--pasa_min_avg_per_id`/`--pasa_min_pct_aligned` pair.

**Resolved (grilling, 2026-08-28):** ship a provisional guess now (~70%
avg-id / ~55% aligned) rather than blocking on empirical calibration first.
A too-strict threshold just degrades to the ab-initio-only fallback — safe,
and consistent with treating these as lower-priority genomes — so there's no
correctness risk in shipping a guess and tuning later against real PASA
output on one hybrid genome.

### Rollout

**Resolved (grilling, 2026-08-28):** on by default once implemented — no
master flag gating it. The risk profile here differs from a pure "skip
RNA-seq" toggle: worst case is the composite evidence gets rejected by the
`composite` PASA tier and the strain falls to ab-initio-only, which is
already today's behavior for these strains in practice (thin/misaligned
shared evidence). A flag would add a step without a matching safety gain.

Still emit an after-the-fact audit report via `collectFile()` (not an ad-hoc
write at the branch point — non-deterministic/`-resume`-unsafe under
Nextflow's execution model), e.g. `misc/hybrid_species_rnaseq_composite.tsv`
(`species_tag, out, asmid, strain, evidence_source` where `evidence_source`
is one of `parents` / `genus_fallback` / `none`) — matches the project's
existing convention of always leaving a reviewable trail (`repr_assignments.tsv`,
`rnaseq_blacklist_candidates.csv`, `rnaseq_se_candidates.csv`).

Manual escape hatch unchanged: `rnaseq_representative_override.csv` still
lets a maintainer force a specific hybrid strain to build its own independent
RNA-seq (Option 4 style, own reads against own genome) instead of using the
composite, for a specific well-studied group deemed worth the extra cost.

### Cleanup required before/alongside this change

`rnaseq_data/` already holds hybrid Trinity artifacts built under today's
un-gated behavior, some already known-bad (`misc/poor_trinity/
Saccharomyces_x_bayanus.trinity-GG.fasta` at 252 B,
`Saccharomyces_eubayanus_x_Saccharomyces_uvarum.trinity-GG.fasta` at 54 KB).
`RNASEQ_PREPARE`'s `storeDir` has no representative-identity invalidation
sentinel (`DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md` Option 2, RNA-seq side
still unimplemented) — any path that re-enters `RNASEQ_PREPARE` for a hybrid
`species_tag` (there shouldn't be one after this change, but old cached
artifacts don't know that) would silently re-serve the stale/near-empty
artifact. Quarantine (don't delete) existing
`rnaseq_data/<hybrid species_tag>.trinity-GG.fasta` files as part of landing
this change.

## Ab-initio parameter reuse: no code change needed

**Resolved (grilling, 2026-08-28), correcting an earlier draft.** An earlier
revision of this plan proposed excluding `HYBRID=true` rows from
`pick_representative_strain.py`'s candidate pool entirely — i.e. every one of
the 90 hybrid strains, even the 37 nominally-identical `S. cerevisiae x
S. eubayanus` strains, would train ab-initio parameters fully independently.
That over-read the maintainer's intent.

What the maintainer actually wants: hybrids should never reuse parameters
**across** hybrid crosses or with true parent species (never pool
`S. cerevisiae x S. eubayanus` with `S. cerevisiae x S. kudriavzevii`, or
with plain `S. cerevisiae`), but reuse **within** an identical hybrid-cross
group (the 37 `S. cerevisiae x S. eubayanus` strains sharing one
representative's trained parameters among themselves) is fine and desirable.

Checked directly against the code: `pick_representative_strain.py:198-202`
groups strains with `by_species = defaultdict(list)`, keyed on the raw
`SPECIES` string, verbatim. A hybrid's `SPECIES` value
(`"Saccharomyces cerevisiae x Saccharomyces eubayanus"`) is already a
distinct dictionary key from its parent (`"Saccharomyces cerevisiae"`) and
from every other hybrid combination. **The existing mechanism already
implements exactly the desired scope, with zero code changes required** —
each hybrid cross gets its own representative pick and its own
`reuse_eligible`/`is_representative` assignments, scoped to just that cross,
today, unmodified.

The only thing left is unrelated to *this* axis: the RNA-seq/PASA-evidence
axis above (composite parent transcripts) is what actually needed new
mechanism, precisely because the *transcript-alignment* case doesn't have an
equivalent "already scoped correctly by construction" property — Stage 4's
representative-pick-and-share logic had to be told not to treat a hybrid
`species_tag` as an ordinary species needing its own `RNASEQ_PREPARE` build.

## Still open

1. **`composite` tier threshold calibration** — provisional values ship now;
   validate against at least one real hybrid genome's PASA output when
   convenient, not blocking.
2. **Sequencing relative to `DIVERGENT_REPRESENTATIVE_RNASEQ_PLAN.md`
   Option 1** (representative-path alignment floor) — still relevant: Option
   1 is a general safety net for any representative whose shared evidence
   doesn't actually align, hybrid or not. This plan's composite strategy
   reduces (but doesn't eliminate — the composite could still align poorly
   to a particular mosaic genome) the odds of hitting that failure mode for
   hybrids specifically. Worth deciding implementation order explicitly.
3. **Genus-level fallback quality** — for genera with many unrelated species
   under one `GENUS` value, concatenating *all* of them could dilute the
   composite with irrelevant transcripts. Not a concern for the current
   dataset (every affected genus here is small/tight), but worth a sanity
   check (e.g. cap by ANI-relatedness to the hybrid genome, or just to the
   genuinely close species) if this generalizes to a genus with broader
   membership.
