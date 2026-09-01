# How siblings get transcript-trained even with `predict_scope: representative_only`

Context: running `funannotate.nf` with `params_predict_representatives.yaml`
(`predict_scope: "representative_only"`) still runs `FUNANNOTATE_TRAIN` for
*every* strain of a species (e.g. all ~53 `Rhizopus_arrhizus` strains), not
just the ANI-selected representative (`Rhizopus_arrhizus_Y5`). This is
expected behavior, not a bug. Below is why, and how the single shared
`rnaseq_data/<species>.trinity-GG.fasta` gets built without contention
between strains.

## Two separate gates, and `predict_scope` only controls one of them

`predict_scope` is read in exactly one place:
`subworkflows/local/FUNANNOTATE_PREDICTION.nf:45`. It gates only the
**gene-calling/predict step** (`FUNANNOTATE_PREDICT`) and the ab-initio
param backfill (`BACKFILL_ABINITIO_PARAMS`). It is never referenced in
`subworkflows/local/FUNANNOTATE_RNASEQ.nf`, which owns `FUNANNOTATE_TRAIN`.

So `--predict_scope representative_only` restricts *predict* to
representatives, but *train* runs on every strain regardless — that's by
design, per the header comment in `FUNANNOTATE_RNASEQ.nf:1-9`: "every other
strain of that species runs `FUNANNOTATE_TRAIN --trinity` against the
archived assembly."

Concretely, for each strain the code just checks
(`FUNANNOTATE_RNASEQ.nf:236-239`) whether that strain already has a
non-empty `funannotate_train.pasa.gff3`, or whether its RNA-seq/genome
inputs are stale — nothing there checks `is_representative` or
`predict_scope`. So every strain that has reads and no existing PASA output
gets trained, even though only the representative strain will actually get
predicted in a `representative_only` pass.

Representative assignment itself comes from
`genome_annotation/_reuse_assignments/abinitio_reuse_assignments.csv`
(one row per strain, `is_representative` True for exactly one strain per
species — picked by ANI + BUSCO completeness/N50 in
`PICK_REPRESENTATIVE_STRAIN`).

## Who builds `rnaseq_data/<species>.trinity-GG.fasta` — no race

Only one job ever builds it. In `FUNANNOTATE_RNASEQ.nf:184-192`:

```groovy
def repr_ch = assembly_with_reads
    .groupTuple(by: 0)                       // collapse ALL strains of the species into one row
    .map { species_tag, outs, asmids, ... ->
        def repIdx = outs.findIndexOf { out -> abinitioReuseMap[out]?.is_representative }
        def i = repIdx >= 0 ? repIdx : 0      // fallback to first strain if no rep assigned
        tuple(species_tag, outs[i], ...)      // emit exactly ONE tuple per species
    }
```

`groupTuple` merges every strain of a species into a single channel item
*before* `RNASEQ_PREPARE` is called, and the `.map` picks out just the
`is_representative` row from inside that group (falling back to index 0 if
`abinitioReuseMap` has no assignment for that species — e.g.
`--run_ani_reuse false`). So `RNASEQ_PREPARE` — which runs
`funannotate train --stop_after_trinity` and archives the Trinity-GG
FASTA — only ever fires once per species, on the representative's reads.
There's no contention between strains for that file.

That single Trinity FASTA is then joined back to every strain by
`species_tag` (`.combine(shared_ch, by: 0)` at line 215), which is why all
the siblings' `FUNANNOTATE_TRAIN` calls end up using the same archived
assembly — they're all reusing the representative's Trinity build, just
running their own genome-specific PASA alignment against it.

## Bottom line

Nothing is broken. `predict_scope: representative_only` does exactly what
it says — restricting the expensive predict/backfill step to the
representative strain in that pass. The Trinity build is correctly
single-sourced from the one representative. The wide `FUNANNOTATE_TRAIN`
fan-out across all strains is intentional pre-staging: by the time
`params_predict_all.yaml` runs for the siblings, their PASA training will
already be cached and predict can proceed immediately without waiting on
`train` serially per strain.

If that fan-out is *not* desired during a representative-only pass (e.g.
to save compute/time until you're sure you'll predict the siblings later),
that would require adding a `predict_scope`-aware filter into
`FUNANNOTATE_RNASEQ.nf`'s `train_todo`/`branched` logic — no such option
exists currently.
