
<!-- README.md is generated from README.Rmd. Please edit that file -->

# L1000PAI

`L1000PAI` is an R package for **direction-aware perturbation activity
inference** from ranked transcriptomic signatures using L1000-style
perturbation-signature libraries.

The package supports:

- **paired libraries** with `Up` / `Down` signatures
- **single-signature libraries** with one consensus signature per
  perturbagen
- **one-step analysis** through `run_l1000_pai()`
- **L1000 signature libraries** downloaded and verified by the companion
  package `L1000PAIdata`
- **ExpressionSet-to-rank pipeline** for direct use after transcriptomic
  preprocessing
- **library-aware score direction handling** for drug-like versus
  target-activity libraries

## Installation

For a published 1.0.0 release, use the versioned tags below. An
unpublished candidate can instead be installed from its two source
tarballs (data package first); retain the accompanying SHA256 manifest
to identify the exact revision.

``` r
# install.packages("remotes")
remotes::install_github("yilin1107/L1000PAIdata@v1.0.0")
remotes::install_github("yilin1107/L1000PAI@v1.0.0")

# Download and verify the CRISPR knockout and chemical perturbation libraries
# (about 32 MB, once per machine):
L1000PAIdata::download_l1000pai_data(c("crispr", "chem"))
```

`L1000PAI` depends on Bioconductor packages (`clusterProfiler`, `limma`,
`Biobase`, `AnnotationDbi`, `org.Hs.eg.db`, `fgsea`). If they are
missing, install them first with `BiocManager::install()`.

The signature libraries are not distributed with either package.
`L1000PAIdata` fetches them from their public source and accepts a file
only when its content matches the reference release used for the
L1000PAI analyses; see the [`L1000PAIdata`
README](https://github.com/yilin1107/L1000PAIdata) for sources, terms of
use and data citations.

Public downloads cover CRISPR, Chem, and the kinase/GPCR Up, Down and
combined files (five library IDs). Ligand, overexpression and
single-drug libraries need verified local copies;
`download_l1000pai_data("all")` does not provide those three local-only
resources.

## Main ideas

`L1000PAI` supports two related perturbation-analysis settings.

### 1. Paired libraries

For paired libraries, perturbagen activity is summarized using a
direction-aware score:

- `NES_Down - NES_Up`
- or `NES_Up - NES_Down`

The direction is controlled by `score_direction` and can be resolved
automatically according to the library type.

Paired inference separates two statistical questions:

- `padj_coherent` is the strict directional intersection-union test:
  both Up and Down arms must support the same activity direction.
- `activity_gls_z` and `p_gls_global` quantify a shared latent activity
  when synchronized Up/Down null draws are supplied.

The GLS p-value is deliberately not used as a replacement for
`padj_coherent`. A strong single arm may support latent activity without
proving that both arms are coherent. Legacy `activity_score`,
`p_combined`, and `padj_combined` columns remain available for
comparisons.

Typical defaults are:

- `crispr` -\> `down_minus_up`
- `single_drug`, `chem`, `ligand`, `oe`, `kinase_gpcr` -\>
  `up_minus_down`

### 2. Single-signature libraries

For single-signature libraries, the package summarizes and classifies
enrichment directly from a single `NES` per perturbagen.

### Sample-level scores and candidate migration

`score_samples_l1000pai()` is a separate within-sample estimator, not
the preranked group-contrast GSEA test above. It requires unique gene
identifiers, fixed preprocessing and the same measured gene universe
when comparing scores. Duplicated gene identifiers and constant samples
are rejected explicitly. The ssGSEA option uses rank-weighted
enrichment, including deterministic tied rank handling. Results from the
earlier raw-expression-weighted implementation must be recomputed; do
not pool scores from the two implementations. Canonical gene ordering
also changes the sequence of seeded null draws for some earlier inputs,
so recalculate prior sample-level null-based outputs. Within-sample
calculation does not by itself establish cross-platform predictive
validity, drug sensitivity or clinical benefit. This software is for
research.

## Core workflow

The main one-step wrapper is `run_l1000_pai()`.

### Minimal paired example with a user-supplied TERM2GENE table

``` r
library(L1000PAI)

de_tbl <- data.frame(
  gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
  stat = c(3, 2, 1, -1, -2, -3),
  stringsAsFactors = FALSE
)

term2gene_paired <- data.frame(
  term = c(
    "DrugA_Up", "DrugA_Up",
    "DrugA_Down", "DrugA_Down",
    "DrugB_Up", "DrugB_Up",
    "DrugB_Down", "DrugB_Down"
  ),
  gene = c(
    "GeneE", "GeneF",
    "GeneA", "GeneB",
    "GeneA", "GeneB",
    "GeneE", "GeneF"
  ),
  stringsAsFactors = FALSE
)

res <- run_l1000_pai(
  x = de_tbl,
  gene_col = "gene",
  stat_col = "stat",
  term2gene = term2gene_paired,
  library_mode = "paired",
  score_direction = "down_minus_up",
  min_gs_size = 2,
  max_gs_size = 500,
  pvalue_cutoff = 1,
  seed = TRUE,
  combine_p_method = "max",
  activity_threshold = 0.5,
  classify_p_value_threshold = 0.25
)

res
```

### Extract results

``` r
result_mode(res)
head(gsea_table(res))
head(result_table(res))
pai_metadata(res)
```

## Stepwise workflow

``` r
ranks <- rank_genes(de_tbl, gene_col = "gene", stat_col = "stat")
ranks
```

``` r
gsea_res <- run_l1000_gsea(
  ranks = ranks,
  term2gene = term2gene_paired,
  min_gs_size = 2,
  max_gs_size = 500,
  pvalue_cutoff = 1,
  seed = TRUE,
  by = "fgsea"
)

head(gsea_res)
```

``` r
pair_res <- collapse_pairs(
  gsea_table = gsea_res,
  combine_p_method = "max",
  score_direction = "down_minus_up"
)

pair_res
```

``` r
classified_res <- classify_activity(
  pair_table = pair_res,
  activity_threshold = 0.5,
  p_value_threshold = 0.25,
  score_direction = "down_minus_up"
)

classified_res
```

## Covariance-aware paired scoring

`score_joint_activity()` estimates a separate 2-by-2 null covariance for
every perturbagen. Its `joint_null` input must contain synchronized
draws with this schema:

``` text
perturbagen | permutation_id | NES_Up | NES_Down
```

For example, after generating the null draws with a common permutation
index:

``` r
pair_gls <- score_joint_activity(
  pair_table = pair_res,
  joint_null = joint_null,
  score_direction = "down_minus_up",
  min_null_draws = 100,
  ridge = 0
)

pair_gls[, c(
  "perturbagen", "activity_score", "activity_gls_z",
  "p_gls_global", "padj_coherent", "heterogeneity_Q"
)]
```

Do not construct `joint_null` from the observed perturbagen rows. It
must come from common gene-label, sample-label, block-respecting, or
residual permutations. The default fgsea summary does not expose
synchronized null draws, so the package will not silently assume zero
Up/Down correlation.

Strict calls with a GLS score threshold can be made without changing the
conjunction evidence column:

``` r
classified_gls <- classify_activity(
  pair_gls,
  score_col = "activity_gls_z",
  p_col = "padj_coherent",
  activity_threshold = 1,
  p_value_threshold = 0.05,
  score_direction = "down_minus_up"
)
```

## ExpressionSet workflow

`L1000PAI` also provides helpers for moving from an `ExpressionSet` to a
ranked gene signature suitable for L1000PAI.

### Inspect an `ExpressionSet`

``` r
inspect_eset(eset)
```

### Convert an `ExpressionSet` to `rank_input`

``` r
res_eset <- eset_to_rank_input(
  eset = eset,
  group_col = "source_name_ch1",
  control_label = "MCF7_unstimulated",
  case_label = "MCF7_Nutlin stimulated",
  stat_col = "t",
  use_voom = "auto"
)

head(res_eset$rank_input)
head(res_eset$tt_gene[, c("SYMBOL", "logFC", "adj.P.Val", "t")], 20)
```

### Run L1000PAI directly from `rank_input`

``` r
res_pai <- run_l1000_pai(
  x = res_eset$rank_input,
  gene_col = "gene",
  stat_col = "stat",
  library = "chem",
  library_mode = "auto",
  score_direction = "auto",
  by = "DOSE",
  nPerm = 1000,
  min_gs_size = 5,
  max_gs_size = 500,
  pvalue_cutoff = 1,
  seed = TRUE
)

head(result_table(res_pai), 20)
```

## Signature libraries from `L1000PAIdata`

`L1000PAI` can analyze any user-supplied `TERM2GENE` table. The L1000
libraries are provided by the companion package `L1000PAIdata`:

``` r
available_l1000_libraries()
L1000PAIdata::download_l1000pai_data("all")   # crispr, chem and kinase/GPCR
L1000PAIdata::list_l1000pai_data_files()
```

The `ligand`, `oe` and `single_drug` libraries have no confirmed public
source and are not downloaded; a local copy whose content matches the
reference release can be added with
`L1000PAIdata::register_l1000pai_data_file()`.

## Switching across libraries

The same ranked transcriptomic input can be analyzed against different
libraries by switching the `library` argument.

``` r
res_crispr <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "crispr",
  library_mode = "auto",
  score_direction = "auto"
)

res_chem <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "chem",
  library_mode = "auto",
  score_direction = "auto"
)

res_kg <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "kinase_gpcr",
  library_mode = "auto",
  score_direction = "auto"
)
```

Individual libraries can also be retrieved directly:

``` r
gmt_crispr <- get_l1000_library("crispr")
gmt_chem <- get_l1000_library("chem")
gmt_kg <- get_l1000_library("kinase_gpcr")
```

## Paired versus single-signature modes

`L1000PAI` can automatically detect library structure:

``` r
term2gene_single <- data.frame(
  term = c("TP53", "TP53", "RB1", "RB1"),
  gene = c("GeneA", "GeneB", "GeneE", "GeneF"),
  stringsAsFactors = FALSE
)

detect_library_mode(term2gene_paired)
detect_library_mode(term2gene_single)
```

When `library_mode = "auto"`, `run_l1000_pai()` chooses:

- **paired mode** for `Up` / `Down` libraries
- **single mode** for consensus-signature libraries

## Result-table utilities

To make downstream merging and plotting easier, paired result tables
also include a unified `term` column equal to `perturbagen`.

Utility helpers are also available for compareCluster-like workflows:

``` r
tbl <- result_table(res_pai)

tbl2 <- add_cluster_column(tbl)
head(tbl2[, c("Cluster", "term")])
```

``` r
final_merged_df <- prepare_cluster_result(
  df = single_drug_df,
  other = multi_drug_df
)
```

These helpers strip common directional suffixes such as:

- `-up` / `-dn`
- `_Up` / `_Down`
- `Up` / `Down`

## GSEA backend notes

`run_l1000_pai()` supports both:

- `by = "fgsea"`
- `by = "DOSE"`

Recommended usage:

- use `by = "fgsea"` for typical runs
- use `by = "DOSE", nPerm = 1000` when needed for backend compatibility
- avoid manually passing `nPerm` to `fgsea` mode

## Interpretation

For paired libraries, the interpretation of `Activated` / `Inhibited`
depends on score direction.

For drug-like libraries such as `chem` or `single_drug`:

- `Activated` usually indicates **concordant perturbation**
- `Inhibited` usually indicates **reverse perturbation**

For target-oriented libraries such as `crispr`, the same labels are
typically interpreted more directly as inferred target-state
relationships.

For single-signature libraries, default classes include:

- `Positive_enrichment`
- `Negative_enrichment`
- `Not_significant`

## Versioned Chem-family taxonomy

Chem family taxonomy v2 is packaged alongside, and does not replace, the
legacy one-row-per-compound annotation catalog. The exact-term registry
covers all 5,425 paired Chem perturbagens; 80 terms have mechanism-level
curation in this release. Direct mechanisms (`direct_moa`) are kept
separate from broader pathway context (`pathway_context`).

``` r
families <- get_chem_family_catalog("family_registry", version = "2.0.0")

chem_long <- annotate_chem_families(
  chem_result,
  family_level = "direct_moa",
  evidence_tier = c("A", "B")
)
```

`annotate_chem_families()` returns a long many-to-many mapping and
performs no family aggregation or enrichment. Downstream analyses should
use the frozen `analysis_unit_id` to give aliases or members of one drug
program a single vote. Unreviewed library terms remain explicit rather
than receiving inferred mechanisms.

## Suggested downstream analyses

Typical downstream uses include:

- ranking candidate perturbagens by `activity_gls_z`, legacy
  `activity_score`, or `NES`
- filtering strict paired calls on `padj_coherent`
- separating `Activated` and `Inhibited` perturbagens
- extracting perturbagen genes for downstream ORA or GSEA
- comparing direct DEG-based enrichment with perturbation-inferred
  regulator-level enrichment
- switching across the `crispr`, `chem` and `kinase_gpcr` libraries
- visualizing KEGG pathways with `pathview`

## Citation

``` r
citation("L1000PAI")
```

When you use the L1000 libraries, also cite the data sources listed in
the [`L1000PAIdata` README](https://github.com/yilin1107/L1000PAIdata).

## License

The package code is released under the MIT license. The signature
libraries are not part of the package and remain under the terms of
their providers.
