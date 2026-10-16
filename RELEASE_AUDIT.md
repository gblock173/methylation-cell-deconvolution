# Initial release audit

Author: Yujin Jeon

This is a fresh code-only Git history. Original history, biological data,
reference matrices, cohort metadata, sample-level results, notebooks, images,
logs, scheduler submissions and machine-specific configuration are excluded.

The release contains the portable baseline NNLS implementation, parameterized
atlas builder, balanced marker selection, five reusable fragmentomics scripts,
synthetic generation and tests, MIT license and privacy checks.

Verification:
- Baseline extraction, solver, R2 and inference methods have identical syntax
  trees after removing location information from the comparison.
- Twenty-five synthetic fit comparisons produced identical numeric results.
- Tiny synthetic atlas and feature-selection comparisons matched the baseline.
- Nine synthetic tests passed: mixtures, zero signal, missing data, minimum
  coverage, interval boundaries, reference order, CLI and tiny atlas building.
- Seven synthetic privacy detector probes passed. The package built and
  installed offline, and the installed command line interface passed its check.
- All five optional fragmentomics script help commands passed.
- The 24-file release passed a comparison against 508 source identifier tokens,
  known private terms, paths, server addresses, contact and credential patterns.
- Exact file allowlisting excludes all raw and derived data and binary files.
- Git author and committer names are Yujin Jeon; email fields are empty.

All example values are generated mathematically. User-supplied files and
outputs can contain sensitive information and are excluded from the allowlist.
The legacy atlas CpG slice convention is retained; see README.md.

The checks describe this initial release. Future edits and any change from
private to public visibility require a fresh content and history review.
