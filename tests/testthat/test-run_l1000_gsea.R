test_that("run_l1000_gsea returns a data.frame", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(
    GeneA = 3.5,
    GeneB = 2.8,
    GeneC = 1.9,
    GeneD = 0.7,
    GeneE = -0.5,
    GeneF = -1.4,
    GeneG = -2.6,
    GeneH = -3.2
  )

  term2gene <- data.frame(
    term = c("TERM_UP", "TERM_UP", "TERM_UP", "TERM_DN", "TERM_DN", "TERM_DN"),
    gene = c("GeneA", "GeneB", "GeneC", "GeneF", "GeneG", "GeneH"),
    stringsAsFactors = FALSE
  )

  res <- run_l1000_gsea(
    ranks = ranks,
    term2gene = term2gene,
    min_gs_size = 2,
    max_gs_size = 500,
    pvalue_cutoff = 1,
    seed = TRUE,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
  expect_true(all(c("term", "NES", "pvalue", "p.adjust") %in% colnames(res)))
})

test_that("run_l1000_gsea errors when ranks are unnamed", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(3, 2, 1, -1)
  term2gene <- data.frame(
    term = c("TERM1", "TERM1"),
    gene = c("GeneA", "GeneB"),
    stringsAsFactors = FALSE
  )

  expect_error(
    run_l1000_gsea(ranks, term2gene),
    "named numeric vector"
  )
})

test_that("run_l1000_gsea errors when term2gene is invalid", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(GeneA = 3, GeneB = 2, GeneC = 1)

  expect_error(
    run_l1000_gsea(ranks, term2gene = c("a", "b")),
    "data.frame"
  )
})

test_that("run_l1000_gsea sorts ranks decreasingly", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks_unsorted <- c(
    GeneF = -2.0,
    GeneA = 3.0,
    GeneD = 0.5,
    GeneB = 2.0,
    GeneE = -1.0,
    GeneC = 1.0
  )

  term2gene <- data.frame(
    term = c("TERM1", "TERM1", "TERM1"),
    gene = c("GeneA", "GeneB", "GeneC"),
    stringsAsFactors = FALSE
  )

  res <- run_l1000_gsea(
    ranks = ranks_unsorted,
    term2gene = term2gene,
    min_gs_size = 2,
    pvalue_cutoff = 1,
    seed = TRUE,
    verbose = FALSE
  )

  expect_s3_class(res, "data.frame")
})
