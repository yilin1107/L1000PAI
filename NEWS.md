# L1000PAI 1.1.0

## One-call report: `pai_deg_report()`

- `pai_deg_report()` turns a differential-expression table into a PDF report
  of up to 19 pages written for non-specialists, an Excel workbook of every table
  behind it (with a column dictionary) and an RDS file of the R objects. The
  gene column (symbols, Ensembl or Entrez IDs) and the ranking statistic
  (limma `t`, DESeq2 `stat`, or sign(log fold change) x -log10 P for edgeR)
  are found automatically; duplicate genes keep the largest |statistic|.
- Both L1000 libraries are scored at the package settings; pathways come from
  `pai_enrich()` or, when no term is significant, `pai_gsea_pathways()`, and
  are shown in depth with every enrichment display of the package
  (`pai_enrichplot()`, `pai_gobubble()`, `pai_term_cnetplot()`,
  `pai_gocircle()`, `pai_gochord()`, `pai_gobar()`, `pai_govenn()`,
  `pai_gocluster()`, `pai_goheat()`, `pai_armmap()`, `pai_hubplot()`,
  `pai_emapplot()`, `pai_chordplot()`, `pai_enrich_cnetplot()`;
  `enrichment = "core"` keeps three pages); the top knockout and compound
  calls are shown with `pai_evidence_card()` and `pai_cnetplot()`;
  compound mechanisms and evidence tiers from taxonomy 3.0.1.
- Seven quality checks (OK, Check, Caution): gene coverage, number of calls,
  Discordant share, a proliferation signal in the input (MSigDB Hallmark E2F
  targets and G2M checkpoint, shipped with the package under CC BY 4.0),
  common-hit compounds (a reference of how each compound usually ranks in 12
  public GEO perturbation contrasts, shipped in `inst/extdata/report` with the
  contrasts listed), and called genes absent from the input.
- New imports: openxlsx, zip and DBI.
- The patchwork namespace is now imported, so `pai_gocircle()`, `pai_gochord()`,
  `pai_goheat()` and the other plots that combine panels with `|` or `/` work
  when they are the first plot of a session (they failed with "Can't find
  method for generic `|`" unless patchwork had been loaded).
- `get_pai_genesets()` builds the GO gene sets with one SQL query per ontology
  on the org.Hs.eg.db tables. The files are byte-identical to the earlier
  `select()`-based build, which needed about 6 GB of memory; the new build
  takes about 15 seconds and under 1 GB.

## Chem compound annotation: taxonomy 3.0.1

- Every Chem library compound now has a mechanism call built from public
  databases and matched by structure (LINCS 2020 InChIKey to ChEMBL): primary
  target genes and action from ChEMBL drug mechanisms, the Guide to
  PHARMACOLOGY primary targets, the LINCS 2020 compound annotation or ChEMBL
  bioactivity, with an evidence tier (`A` curated, `B` bioactivity, `C`/`D`
  earlier label only, `U` unknown, `X` not a small molecule), a confidence
  flag (`corroborated`, `single`, `conflict`), three class levels and the
  record identifiers behind the call. 3,283 compounds have an A or B call.
- `get_l1000_annotation_catalog("chem")`, `annotate_l1000_result()`,
  `build_l1000_annotation_tables()` and `write_l1000_annotation_tables()` use
  the 3.0.1 call by default; chem summaries are now `broad_class_summary`,
  `family_class_summary` and `target_class_summary`. The earlier catalog is
  available with `version = "legacy"`.
- `get_chem_family_catalog()`, `annotate_chem_families()` and
  `verify_chem_family_catalogs()` default to taxonomy 3.0.1. The v2 families
  are kept and every membership is graded (`A` rule met and two sources agree;
  `B` curated v2 member not contradicted; `C` single source; `U` contradicted).
  The default `evidence_tier = c("A", "B")` keeps all curated v2 members and
  adds corroborated new ones. `MICROTUBULE_TARGETING` joins the two microtubule
  families, whose split the databases do not record. New components:
  `term_annotation` and `evidence`. The frozen v2 taxonomy remains available
  with `version = "2.0.0"`.
- Version 3.0.1 applies author decisions on top of the database rules
  (`get_chem_family_catalog("curation_log")`). Mubritinib, recorded as an
  ERBB2 inhibitor in ChEMBL and LINCS, is assigned to mitochondrial complex I
  (Baccelli et al., Cancer Cell 2019; Stephenson et al., eLife 2020) and
  removed from `ERBB_TKI_CORE`. Gene symbols that are aliases of several genes
  are no longer resolved to the first one; LINCS `ND1` is read as `MT-ND1`.
- The 3.0.1 data files are released under CC BY-SA 4.0 (adapted from ChEMBL
  and the Guide to PHARMACOLOGY); see
  `inst/extdata/annotation_catalogs/v3/LICENSE_DATA.md`. No LINCS source text
  is redistributed.

## Plotting

- New plotting functions draw an L1000PAI result in one visual language: colour
  means direction only, and explanations are drawn inside the figure. Every
  function returns a ggplot or patchwork object, attaches the values drawn as
  `attr(p, "plot_data")` and a recommended size used by `pai_save()`. See
  `?pai_plotting`.
- Targets: `pai_dotplot()`, `pai_armplot()`, `pai_gseaplot()`,
  `pai_evidence_card()`, `pai_compareplot()`, `pai_cnetplot()`,
  `pai_familyplot()` and the one-page `pai_report()`.
- Pathways: `pai_enrich()` tests the Activated and the Inhibited targets
  separately against GO, KEGG, Reactome, the Elsevier Pathway Collection or a
  GMT file, with the scored targets as background; `pai_gsea_pathways()` ranks
  all scored targets by activity score. Plots: `pai_enrichplot()`,
  `pai_pathwayplot()`, the GOplot-style `pai_gobubble()`, `pai_gobar()`,
  `pai_gocircle()`, `pai_gochord()`, `pai_goheat()`, `pai_govenn()` and
  `pai_gocluster()`, and the two-arm displays `pai_armmap()`, `pai_ladder()`
  and `pai_term_cnetplot()`. Pathway plots take the result of `pai_enrich()`.
- Hubs: `pai_enrich_hubs()` groups redundant terms into modules;
  `pai_hubplot()`, `pai_chordplot()`, `pai_enrich_cnetplot()` and
  `pai_emapplot()` draw them.
- Style: `theme_pai()`, `pai_use_palette()`, `pai_style_guide()` and
  `pai_save()`; fonts are set with `options(pai.font = , pai.mono = )`.
- Save PDF files with `pai_save()` or another cairo device
  (`ggsave(..., device = cairo_pdf)`): the labels use arrows, a dagger and the
  greater-or-equal sign, which the base `pdf()` device replaces with dots.
- `n` means the top n everywhere. For `pai_evidence_card()` it is the number of
  rows in total, shared between the two directions in proportion to their
  numbers of calls.
- The z-score of the GOplot-style displays defaults to a balanced z-score that
  allows for the different sizes of the two call lists (`z = "goplot"` gives
  GOplot's).

## Dependencies

- `ggplot2` (>= 4.0.0), `patchwork`, `scales`, `ggnewscale`, `ggrepel`,
  `ggforce`, `igraph`, `GO.db` and `digest` are now in `Imports`; all of them
  are already installed with `clusterProfiler`.

---

# L1000PAI 1.0.0

First public release. It contains the 0.3 development series below.

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
