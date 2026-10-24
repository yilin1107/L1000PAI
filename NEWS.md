# L1000PAI 1.0.0

First public release. It contains the 0.3 development series below.

## Candidate revision 2 (unreleased)

- Sample-level scoring rejects duplicated gene identifiers instead of selecting
  different rows in single-sample and cohort runs. Fix gene mapping upstream.
- Corrected `method = "ssgsea"` to use within-sample rank weights, with a
  rank-weighted reference walk and shift/scale-invariance regression tests.
  Scores from the former raw-expression-weighted implementation must be
  recomputed; they are not numerically interchangeable with this correction.
- Added sample/cohort and gene-row-order invariance checks with a fixed gene
  universe and random seed. Cross-platform predictive validity is not implied.
  Canonical gene ordering can change historical seeded null draws even when
  raw rank-mean scores are unchanged; recompute previous sample-level outputs.
- Added explicit handling of constant samples and degenerate enrichment sets.
- Paired GSEA tables containing missing NES retain unknown scores and states
  without aborting estimable pairs; unknown directions are not discordance.
- The current data provider validates a requested cached library before reading
  it. The companion dependency for release installation/CI is pinned to v1.0.0.
- This candidate has not been published or tagged. Candidate package versions
  remain 1.0.0; use the delivery SHA256 manifest to distinguish revisions.

## Signature libraries

- The L1000 libraries are no longer bundled. With `L1000PAIdata` (>= 1.0.0),
  `L1000PAIdata::download_l1000pai_data()` downloads the `crispr`, `chem` and
  kinase/GPCR libraries from their public source and accepts them only when
  their content matches the reference release used for the L1000PAI analyses.
- `get_l1000_library()` explains how to download a library that is not yet
  available locally.
- Tests and vignette chunks that need a library are skipped when it is not
  available.

## Fixes

- Declared `fgsea` in `Imports` (used for conservative completion of
  prespecified terms) and qualified `stats::setNames()`.
- Documented `score_samples_l1000pai()` and `sample_score_table()`.

## Packaging

- Added `inst/CITATION`, `URL`, `BugReports` and `Remotes` fields.
- Removed the GSE47042 worked example from the README and vignette; its
  scripts and data were not part of the package.

---

# L1000PAI 0.3.0.9000 (development)

## Versioned Chem-family taxonomy

- Added a frozen Chem-family taxonomy v2 with exact coverage of all 5,425
  paired Chem terms and curated many-to-many memberships for the 80 terms used
  in the historical ten-case mechanism analysis.
- Separated strict direct-mechanism families from pathway-context families,
  retained the legacy annotation API unchanged, and added
  `get_chem_family_catalog()`, `annotate_chem_families()`, and
  `verify_chem_family_catalogs()`.
- Added explicit canonical analysis units for Tanespimycin/17-AAG and the
  Ixazomib/MLN-9708 drug program, and restored the omitted `PLX-4720` paired
  term without outcome-based family selection.

## Paired statistical inference

- Added explicit directional intersection-union tests for paired Up/Down
  signatures. Activation and inhibition hypotheses are adjusted together, and
  `padj_coherent` is now preferred for strict paired classification.
- Added `score_joint_activity()` for pair-specific, covariance-aware GLS
  scoring from synchronized Up/Down null draws. It reports empirical latent
  activity evidence separately from strict conjunction evidence.
- Preserved `activity_score`, `p_combined`, and `padj_combined` for backward
  compatibility, with `activity_score_legacy` added as an explicit audit field.
- Added `score_col` support to `classify_activity()` and `plot_l1000_pai()` so
  GLS Z scores can be used without overwriting legacy NES-difference scores.
- Added a default `require_coherent` guard so explicitly selecting a global GLS
  or legacy p-value cannot bypass the directional conjunction gate.
- Updated annotation summaries to prefer GLS Z scores and coherent p-values
  when present, while recording the score source used.
- Named the covariance-based heterogeneity tail probability
  `p_heterogeneity_asymptotic` to make its calibration explicit.
- Added an explicit `rng_seed` option for reproducible GSEA execution.
- Added opt-in `keep_all_terms` and `complete_prespecified_terms` controls for
  fixed-universe audits. Non-estimable or size-ineligible prespecified arms are
  retained conservatively as p=1 and NES=0 without relaxing gene-set-size
  thresholds; these rows cannot create a directional formal call.
- Added validation for invalid p-values, duplicate directional signatures, and
  incomplete pairs.

## Statistical scope

- GLS p-values are labelled global/latent evidence and never replace the
  directional conjunction gate by default.
- Joint GLS scoring requires synchronized null draws. The package does not
  estimate null covariance from observed perturbagen results or silently assume
  independent Up/Down arms.

---

# L1000PAI 0.2.0.9000

## Project structure

- Separated the core R package from the new sibling `L1000PAI-analysis`
  workspace for case configuration, workflows, provenance, results and
  manuscript sources.
- Added an allowlist-based paired package check script so Windows builds do not
  traverse the historical multi-gigabyte `submission/` tree.
- Added cross-platform core-package CI and explicit build-artifact exclusions.
- Declared compatibility with `L1000PAIdata >= 0.1.0` to prevent silently using
  an older, structurally incompatible companion package.

## Development status

- The historical `submission/` directory is frozen and remains in its original
  location pending checksum-verified migration.
- Version 0.2.0 remains a development version until the full benchmark and data
  governance work is complete.

---

# L1000PAI 0.1.1

## Added

- Added support for both **paired** and **single-signature** L1000 libraries in `run_l1000_pai()`.
- Added automatic library-mode detection via `detect_library_mode()`.
- Added `summarize_single_signatures()` and `classify_single_signatures()` for single-signature enrichment workflows.
- Added `resolve_score_direction()` to support library-aware paired-score direction selection.
- Added helper functions for `ExpressionSet` workflows:
  - `inspect_eset()`
  - `looks_like_raw_counts()`
  - `choose_feature_id_col()`
  - `run_limma_from_eset()`
  - `add_symbol_to_limma_table()`
  - `collapse_limma_to_genes()`
  - `eset_to_rank_input()`
- Added result-table utilities for downstream integration and plotting:
  - `strip_direction_suffix()`
  - `add_cluster_column()`
  - `merge_l1000_results()`
  - `prepare_cluster_result()`
- Added unified `term` column to paired result tables to facilitate downstream merging with single-signature outputs.
- Added `result_table()` and `result_mode()` accessors for standardized extraction of final results.

## Changed

- Updated paired-score calculation to support both:
  - `NES_Down - NES_Up`
  - `NES_Up - NES_Down`
- `run_l1000_pai(score_direction = "auto")` now resolves paired-score direction based on bundled library type:
  - `single_drug`, `chem`, `ligand`, `oe`, and `kinase_gpcr` default to `up_minus_down`
  - `crispr` defaults to `down_minus_up`
- Improved `collapse_pairs()` to recognize more directional suffix styles, including:
  - `-up` / `-dn`
  - `_Up` / `_Down`
  - ` Up` / ` Down`
- Updated `classify_activity()` so directional consistency is evaluated according to the selected score direction.
- Standardized GSEA outputs so downstream paired workflows consistently use a `term` column even when the backend returns `ID`.
- Expanded vignette and README coverage for:
  - bundled libraries from `L1000PAIdata`
  - library switching across `crispr`, `chem`, `ligand`, `oe`, `single_drug`, and `kinase_gpcr`
  - paired vs single-signature interpretation

## Fixed

- Fixed failure to recognize paired signatures when terms used space-separated directional suffixes such as `"CDK6 Down"` or `"WNT7A Up"`.
- Fixed incompatibility between paired result tables and downstream workflows expecting a `term` column.
- Fixed backend-specific GSEA issues by exposing `nPerm` for the `"DOSE"` backend.
- Improved overlap checking between ranked genes and `TERM2GENE` libraries before running GSEA.
- Improved robustness of feature-to-symbol mapping when working with `ExpressionSet` objects lacking explicit `SYMBOL` columns.
- Improved handling of tied preranked statistics through optional stable tie-breaking in `rank_input`.

## Documentation

- Added and expanded vignette examples for:
  - paired-signature libraries
  - single-signature libraries
  - bundled library retrieval from `L1000PAIdata`
  - switching across perturbation-resource types
- Added README examples showing one-step and stepwise workflows.
- Clarified interpretation of paired `Activated` / `Inhibited` versus single-signature `Positive_enrichment` / `Negative_enrichment`.

---

# L1000PAI 0.1.0

## Added

- Initial public release of `L1000PAI`.
- Added core workflow functions:
  - `rank_genes()`
  - `read_l1000_gmt()`
  - `run_l1000_gsea()`
  - `collapse_pairs()`
  - `classify_activity()`
  - `run_l1000_pai()`
- Added support for bundled library retrieval via `get_l1000_library()`.
- Added companion-package integration with `L1000PAIdata`.
- Added S4 result container `L1000PAIResult`.
- Added accessors for rank, GSEA, result table, and metadata extraction.
- Added testthat-based unit tests for ranking, GMT parsing, GSEA, paired summarization, and workflow execution.
- Added vignette and README documentation.

## Notes

- `L1000PAI` is designed for perturbation activity inference from ranked transcriptomic signatures using paired and single-signature L1000-style libraries.
- The package is intended for use with differential-expression outputs from workflows such as `limma` and `DESeq2`.
