# Preflight check for broken `training` symlinks before funannotate predict

**Status**: done 2026-09-10 (registry T-034; code in `modules/funannotate/utils.nf` `trainingSymlinkPreflight`, smoke test `nextflow/training_preflight_smoke.nf`). Written as an idea first. Triggered by a real incident: dangling `training`
symlinks in `genome_annotation/<SPECIES>/training` (meant to point at the
project-local `genome_annotation_training/<SPECIES>/training`) silently
blocked `predict` from finding training data for an unknown but non-trivial
number of species. Fixed reactively via
`find . -maxdepth 1 -xtype l -exec rm {} \;` inside `genome_annotation`, which
strips only the dangling symlinks. See `.living/learnings.md` 2026-09-10 entry
for the full incident writeup.

## Problem

Nothing currently checks, before `predict` runs, whether `training` in a
species' `genome_annotation/<SPECIES>/` directory is:
- a symlink at all (vs. a real directory, or missing entirely),
- pointing at the correct target
  (`genome_annotation_training/<SPECIES>/training`),
- resolving to an existing directory with real training files inside.

A broken or wrong-target link doesn't surface as a clear "missing training
data" error — funannotate/predict either proceeds without training data or
fails downstream in a way that doesn't point back at the symlink. That makes
this class of failure cryptic and easy to miss at BFD's genome count.

## Proposed check (not yet implemented)

Before (or as an early step of) the predict stage, for every species with a
`genome_annotation/<SPECIES>/` entry:

1. If `genome_annotation_training/<SPECIES>/training` exists (i.e. training
   data is expected for this species):
   - `genome_annotation/<SPECIES>/training` must be a symlink (`test -L`).
   - It must resolve (`test -e`, not just `-L`) and point at
     `genome_annotation_training/<SPECIES>/training` specifically (not just
     any existing target).
   - If any of the above fails — dangling symlink, wrong target, a plain
     directory instead of a symlink, or the symlink missing entirely — flag
     the species explicitly as a config/corruption problem (e.g. a
     `predict_training_symlink_broken.tsv` report, mirroring the existing
     `predict_skipped_too_small.tsv` / `predict_blocked_awaiting_representative.tsv`
     pattern) rather than letting predict silently continue.
2. Consider a periodic repo-wide `find genome_annotation -maxdepth 2 -xtype l`
   sweep as a standalone health check, independent of the predict gate, since
   this class of breakage is invisible until predict quietly underperforms.

## Open questions

- Should this be a Nextflow process (own report file, participates in the
  DAG/caching) or a standalone pre-flight script run before launching the
  pipeline?
- After removing the broken links in this incident, do the affected species
  need their symlinks regenerated pointing at
  `genome_annotation_training/<SPECIES>/training`, or does predict retrain
  from scratch and repopulate correctly on its own? Not yet checked.

**Tags**: symlink, genome_annotation, genome_annotation_training, funannotate,
predict, training-data, preflight-check, data-integrity
