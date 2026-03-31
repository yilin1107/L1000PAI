
<!-- README.md is generated from README.Rmd. Please edit that file -->

# L1000PAI

`L1000PAI` is an R package for **direction-aware perturbation activity
inference** from ranked transcriptomic signatures using L1000
perturbation-signature libraries.

The package supports:

- **paired libraries** with `Up` / `Down` signatures
- **single-signature libraries** with one consensus signature per
  perturbagen
- **one-step analysis** through `run_l1000_pai()`
- **bundled perturbation libraries** via the companion data package
  `L1000PAIdata`

## Installation

``` r
# install.packages("devtools")

# install the main package
# devtools::install("path/to/L1000PAI")
```

To use bundled perturbation libraries, also install the companion data
package:

``` r
# devtools::install("path/to/L1000PAIdata")
```

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
  min_gs_size = 2,
  max_gs_size = 500,
  pvalue_cutoff = 1,
  seed = TRUE,
  combine_p_method = "max",
  activity_threshold = 0.5,
  classify_p_value_threshold = 0.25
)
```

    ## 

``` r
res
```

    ## L1000PAIResult
    ##   Ranked genes: 6 
    ##   GSEA rows: 4 
    ##   Final rows: 2 
    ##   Library: user_supplied_term2gene 
    ##   Mode: paired

### Extract results

``` r
result_mode(res)
```

    ## [1] "paired"

``` r
head(gsea_table(res))
```

    ##                  term Description setSize enrichmentScore       NES    pvalue
    ## DrugA_Up     DrugA_Up    DrugA_Up       2              -1 -2.163617 0.0910596
    ## DrugB_Down DrugB_Down  DrugB_Down       2              -1 -2.163617 0.0910596
    ## DrugA_Down DrugA_Down  DrugA_Down       2               1  2.143611 0.1098361
    ## DrugB_Up     DrugB_Up    DrugB_Up       2               1  2.143611 0.1098361
    ##             p.adjust    qvalue rank                     leading_edge
    ## DrugA_Up   0.1098361 0.1098361    3  tags=100%, list=50%, signal=75%
    ## DrugB_Down 0.1098361 0.1098361    3  tags=100%, list=50%, signal=75%
    ## DrugA_Down 0.1098361 0.1098361    2 tags=100%, list=33%, signal=100%
    ## DrugB_Up   0.1098361 0.1098361    2 tags=100%, list=33%, signal=100%
    ##            core_enrichment
    ## DrugA_Up             GeneF
    ## DrugB_Down           GeneF
    ## DrugA_Down     GeneA/GeneB
    ## DrugB_Up       GeneA/GeneB

``` r
head(result_table(res))
```

    ##   perturbagen  term_Up Description_Up setSize_Up enrichmentScore_Up    NES_Up
    ## 1       DrugA DrugA_Up       DrugA_Up          2                 -1 -2.163617
    ## 2       DrugB DrugB_Up       DrugB_Up          2                  1  2.143611
    ##   pvalue_Up p.adjust_Up rank_Up                  leading_edge_Up
    ## 1 0.0910596   0.1098361       3  tags=100%, list=50%, signal=75%
    ## 2 0.1098361   0.1098361       2 tags=100%, list=33%, signal=100%
    ##   core_enrichment_Up  term_Down Description_Down setSize_Down
    ## 1              GeneF DrugA_Down       DrugA_Down            2
    ## 2        GeneA/GeneB DrugB_Down       DrugB_Down            2
    ##   enrichmentScore_Down  NES_Down pvalue_Down p.adjust_Down rank_Down
    ## 1                    1  2.143611   0.1098361     0.1098361         2
    ## 2                   -1 -2.163617   0.0910596     0.1098361         3
    ##                  leading_edge_Down core_enrichment_Down combined_stat
    ## 1 tags=100%, list=33%, signal=100%          GeneA/GeneB            NA
    ## 2  tags=100%, list=50%, signal=75%                GeneF            NA
    ##   p_combined activity_score padj_combined is_significant direction_consistent
    ## 1  0.1098361       4.307228     0.1098361           TRUE                 TRUE
    ## 2  0.1098361      -4.307228     0.1098361           TRUE                 TRUE
    ##   activity_class
    ## 1      Activated
    ## 2      Inhibited

``` r
pai_metadata(res)
```

    ## $input_type
    ## [1] "data.frame"
    ## 
    ## $library
    ## [1] "user_supplied_term2gene"
    ## 
    ## $data_package
    ## [1] "L1000PAIdata"
    ## 
    ## $requested_library_mode
    ## [1] "paired"
    ## 
    ## $detected_library_mode
    ## [1] "paired"
    ## 
    ## $library_mode_used
    ## [1] "paired"
    ## 
    ## $n_ranked_genes
    ## [1] 6
    ## 
    ## $n_term2gene_rows
    ## [1] 8
    ## 
    ## $gsea_params
    ## $gsea_params$exponent
    ## [1] 1
    ## 
    ## $gsea_params$min_gs_size
    ## [1] 2
    ## 
    ## $gsea_params$max_gs_size
    ## [1] 500
    ## 
    ## $gsea_params$eps
    ## [1] 1e-10
    ## 
    ## $gsea_params$pvalue_cutoff
    ## [1] 1
    ## 
    ## $gsea_params$p_adjust_method
    ## [1] "BH"
    ## 
    ## $gsea_params$seed
    ## [1] TRUE
    ## 
    ## $gsea_params$by
    ## [1] "fgsea"
    ## 
    ## $gsea_params$verbose
    ## [1] FALSE
    ## 
    ## 
    ## $pairing_params
    ## $pairing_params$combine_p_method
    ## [1] "max"
    ## 
    ## 
    ## $classify_params
    ## $classify_params$activity_threshold
    ## [1] 0.5
    ## 
    ## $classify_params$classify_p_value_threshold
    ## [1] 0.25
    ## 
    ## $classify_params$classify_p_col
    ## NULL
    ## 
    ## $classify_params$single_nes_threshold
    ## [1] 1

## Stepwise workflow

``` r
ranks <- rank_genes(de_tbl, gene_col = "gene", stat_col = "stat")
ranks
```

    ## GeneA GeneB GeneC GeneD GeneE GeneF 
    ##     3     2     1    -1    -2    -3

``` r
gsea_res <- run_l1000_gsea(
  ranks = ranks,
  term2gene = term2gene_paired,
  min_gs_size = 2,
  max_gs_size = 500,
  pvalue_cutoff = 1,
  seed = TRUE
)

head(gsea_res)
```

    ##                  term Description setSize enrichmentScore       NES    pvalue
    ## DrugA_Up     DrugA_Up    DrugA_Up       2              -1 -2.163617 0.0910596
    ## DrugB_Down DrugB_Down  DrugB_Down       2              -1 -2.163617 0.0910596
    ## DrugA_Down DrugA_Down  DrugA_Down       2               1  2.143611 0.1098361
    ## DrugB_Up     DrugB_Up    DrugB_Up       2               1  2.143611 0.1098361
    ##             p.adjust    qvalue rank                     leading_edge
    ## DrugA_Up   0.1098361 0.1098361    3  tags=100%, list=50%, signal=75%
    ## DrugB_Down 0.1098361 0.1098361    3  tags=100%, list=50%, signal=75%
    ## DrugA_Down 0.1098361 0.1098361    2 tags=100%, list=33%, signal=100%
    ## DrugB_Up   0.1098361 0.1098361    2 tags=100%, list=33%, signal=100%
    ##            core_enrichment
    ## DrugA_Up             GeneF
    ## DrugB_Down           GeneF
    ## DrugA_Down     GeneA/GeneB
    ## DrugB_Up       GeneA/GeneB

``` r
pair_res <- collapse_pairs(
  gsea_table = gsea_res,
  combine_p_method = "max"
)

pair_res
```

    ##   perturbagen  term_Up Description_Up setSize_Up enrichmentScore_Up    NES_Up
    ## 1       DrugA DrugA_Up       DrugA_Up          2                 -1 -2.163617
    ## 2       DrugB DrugB_Up       DrugB_Up          2                  1  2.143611
    ##   pvalue_Up p.adjust_Up rank_Up                  leading_edge_Up
    ## 1 0.0910596   0.1098361       3  tags=100%, list=50%, signal=75%
    ## 2 0.1098361   0.1098361       2 tags=100%, list=33%, signal=100%
    ##   core_enrichment_Up  term_Down Description_Down setSize_Down
    ## 1              GeneF DrugA_Down       DrugA_Down            2
    ## 2        GeneA/GeneB DrugB_Down       DrugB_Down            2
    ##   enrichmentScore_Down  NES_Down pvalue_Down p.adjust_Down rank_Down
    ## 1                    1  2.143611   0.1098361     0.1098361         2
    ## 2                   -1 -2.163617   0.0910596     0.1098361         3
    ##                  leading_edge_Down core_enrichment_Down combined_stat
    ## 1 tags=100%, list=33%, signal=100%          GeneA/GeneB            NA
    ## 2  tags=100%, list=50%, signal=75%                GeneF            NA
    ##   p_combined activity_score padj_combined
    ## 1  0.1098361       4.307228     0.1098361
    ## 2  0.1098361      -4.307228     0.1098361

``` r
classified_res <- classify_activity(
  pair_table = pair_res,
  activity_threshold = 0.5,
  p_value_threshold = 0.25
)

classified_res
```

    ##   perturbagen  term_Up Description_Up setSize_Up enrichmentScore_Up    NES_Up
    ## 1       DrugA DrugA_Up       DrugA_Up          2                 -1 -2.163617
    ## 2       DrugB DrugB_Up       DrugB_Up          2                  1  2.143611
    ##   pvalue_Up p.adjust_Up rank_Up                  leading_edge_Up
    ## 1 0.0910596   0.1098361       3  tags=100%, list=50%, signal=75%
    ## 2 0.1098361   0.1098361       2 tags=100%, list=33%, signal=100%
    ##   core_enrichment_Up  term_Down Description_Down setSize_Down
    ## 1              GeneF DrugA_Down       DrugA_Down            2
    ## 2        GeneA/GeneB DrugB_Down       DrugB_Down            2
    ##   enrichmentScore_Down  NES_Down pvalue_Down p.adjust_Down rank_Down
    ## 1                    1  2.143611   0.1098361     0.1098361         2
    ## 2                   -1 -2.163617   0.0910596     0.1098361         3
    ##                  leading_edge_Down core_enrichment_Down combined_stat
    ## 1 tags=100%, list=33%, signal=100%          GeneA/GeneB            NA
    ## 2  tags=100%, list=50%, signal=75%                GeneF            NA
    ##   p_combined activity_score padj_combined is_significant direction_consistent
    ## 1  0.1098361       4.307228     0.1098361           TRUE                 TRUE
    ## 2  0.1098361      -4.307228     0.1098361           TRUE                 TRUE
    ##   activity_class
    ## 1      Activated
    ## 2      Inhibited

## Companion data package: `L1000PAIdata`

`L1000PAI` can analyze user-supplied `TERM2GENE` tables, but it also
supports bundled libraries distributed through the companion package
`L1000PAIdata`.

Available bundled library names include:

``` r
available_l1000_libraries()
```

    ## [1] "crispr"           "chem"             "ligand"           "oe"              
    ## [5] "single_drug"      "kinase_gpcr"      "kinase_gpcr_up"   "kinase_gpcr_down"

Typical options include:

- `"crispr"`
- `"chem"`
- `"ligand"`
- `"oe"`
- `"single_drug"`
- `"kinase_gpcr"`

## Switching across libraries

The same ranked transcriptomic input can be analyzed against different
bundled libraries simply by switching the `library` argument.

``` r
res_crispr <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "crispr",
  library_mode = "auto"
)

res_chem <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "chem",
  library_mode = "auto"
)

res_ligand <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "ligand",
  library_mode = "auto"
)

res_single_drug <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "single_drug",
  library_mode = "auto"
)

res_kg <- run_l1000_pai(
  x = tt,
  gene_col = "gene",
  stat_col = "t",
  library = "kinase_gpcr",
  library_mode = "auto"
)
```

If `L1000PAIdata` is installed, individual libraries can also be
retrieved directly:

``` r
gmt_crispr <- get_l1000_library("crispr")
gmt_chem <- get_l1000_library("chem")
gmt_ligand <- get_l1000_library("ligand")
gmt_single_drug <- get_l1000_library("single_drug")
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
```

    ## [1] "paired"

``` r
detect_library_mode(term2gene_single)
```

    ## [1] "single"

When `library_mode = "auto"`, `run_l1000_pai()` chooses:

- **paired mode** for `Up` / `Down` libraries
- **single mode** for consensus-signature libraries

## Interpretation

For paired libraries, perturbagen activity is summarized using:

$$
\mathrm{activity\_score} = NES_{Down} - NES_{Up}
$$

Default paired classes include:

- `Activated`
- `Inhibited`
- `Discordant`
- `Not_significant`

For single-signature libraries, default classes include:

- `Positive_enrichment`
- `Negative_enrichment`
- `Not_significant`

## Suggested downstream analyses

Typical downstream uses include:

- ranking candidate perturbagens by `activity_score` or `NES`
- filtering on combined adjusted p-values
- separating `Activated` and `Inhibited` perturbagens
- extracting perturbagen genes for downstream ORA or GSEA
- comparing direct DEG-based enrichment with perturbation-inferred
  regulator-level enrichment
- switching across `crispr`, `chem`, `ligand`, `oe`, `single_drug`, and
  `kinase_gpcr` libraries
