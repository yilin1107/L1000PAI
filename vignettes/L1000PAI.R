## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(L1000PAI)

## -----------------------------------------------------------------------------
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

de_tbl
term2gene_paired

## -----------------------------------------------------------------------------
ranks <- rank_genes(
  x = de_tbl,
  gene_col = "gene",
  stat_col = "stat"
)

ranks

## -----------------------------------------------------------------------------
detect_library_mode(term2gene_paired)

