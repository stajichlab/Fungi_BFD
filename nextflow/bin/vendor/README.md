# Vendored scripts

## `filterIntronsFindStrand.pl`

Source: [`Gaius-Augustus/BRAKER`](https://github.com/Gaius-Augustus/BRAKER)
`scripts/filterIntronsFindStrand.pl`, fetched 2026-08-12 from the `master`
branch. Artistic License (header preserved in the file, full license text at
[opensource.org/licenses/artistic-license.php](http://www.opensource.org/licenses/artistic-license.php)).
Authors: Simone Lange & Katharina J. Hoff.

Not modified from upstream. Used by the GeneMark-ET intron-hints recipe
validated in `nextflow/docs/GENEMARK_RUN_DESIGN.md`'s "ET mode" section and
`analysis/genemark_run_validation/et_eval2/` — checks each intron hint's
splice-site boundaries against the genome sequence to assign strand and
canonicality, a step `gmes_petap.pl --ET`'s branch-point training silently
fails without (unstranded hints produce no usable branch-point signal).

Re-fetch from upstream if BRAKER releases a newer version; this copy is
pinned to the 2026-08-12 `master` snapshot, not a tagged release (BRAKER
doesn't tag this script independently of the whole pipeline).
