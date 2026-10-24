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

test_that("run_l1000_gsea initializes RNG state when seed is TRUE", {
  testthat::skip_if_not_installed("clusterProfiler")

  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }

  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    rm(".Random.seed", envir = .GlobalEnv)
  }

  on.exit({
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

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
  expect_true(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
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

test_that("run_l1000_gsea validates an explicit RNG seed", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(GeneA = 3, GeneB = 2, GeneC = 1, GeneD = -1)
  term2gene <- data.frame(
    term = c("TERM1", "TERM1"),
    gene = c("GeneA", "GeneB"),
    stringsAsFactors = FALSE
  )

  expect_error(
    run_l1000_gsea(ranks, term2gene, rng_seed = -1),
    "between 0 and"
  )
  expect_error(
    run_l1000_gsea(ranks, term2gene, rng_seed = 1.5),
    "between 0 and"
  )
})

test_that("run_l1000_gsea validates complete-universe mode", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(GeneA = 3, GeneB = 2, GeneC = 1, GeneD = -1)
  term2gene <- data.frame(
    term = c("TERM1", "TERM1"),
    gene = c("GeneA", "GeneB"),
    stringsAsFactors = FALSE
  )

  expect_error(
    run_l1000_gsea(ranks, term2gene, keep_all_terms = NA),
    "single non-missing logical"
  )
  expect_error(
    run_l1000_gsea(ranks, term2gene, keep_all_terms = 1),
    "single non-missing logical"
  )
  expect_error(
    run_l1000_gsea(ranks, term2gene, complete_prespecified_terms = NA),
    "single non-missing logical"
  )
  expect_error(
    run_l1000_gsea(
      ranks, term2gene, keep_all_terms = TRUE, pvalue_cutoff = 0.05
    ),
    "pvalue_cutoff = 1"
  )
})

test_that("prespecified completion retains size-ineligible terms conservatively", {
  testthat::skip_if_not_installed("clusterProfiler")

  genes <- paste0("Gene", seq_len(30))
  ranks <- stats::setNames(seq(3, -3, length.out = 30), genes)
  term2gene <- data.frame(
    term = c(rep("TERM_A", 10), rep("TERM_B", 10), rep("TERM_SMALL", 2)),
    gene = c(genes[1:10], genes[21:30], genes[14:15]),
    stringsAsFactors = FALSE
  )

  result <- suppressWarnings(run_l1000_gsea(
    ranks,
    term2gene,
    min_gs_size = 3,
    max_gs_size = 500,
    pvalue_cutoff = 1,
    rng_seed = 2026L,
    keep_all_terms = TRUE,
    complete_prespecified_terms = TRUE
  ))
  small <- result[result$term == "TERM_SMALL", , drop = FALSE]

  expect_equal(nrow(small), 1L)
  expect_equal(small$setSize, 2L)
  expect_equal(small$NES, 0)
  expect_equal(small$pvalue, 1)
  expect_equal(small$p.adjust, 1)
  expect_true(is.na(small$rank))
})

test_that("non-estimable fgsea terms receive conservative fallback rows", {
  testthat::skip_if_not_installed("fgsea")

  ranks <- c(GeneA = 3, GeneB = 2, GeneC = 1, GeneD = -1, GeneE = -2)
  term2gene <- data.frame(
    term = c("TERM1", "TERM1"),
    gene = c("GeneA", "GeneC"),
    stringsAsFactors = FALSE
  )
  eligible <- .eligible_l1000_terms(
    names(ranks), term2gene, min_gs_size = 2, max_gs_size = 500
  )
  fallback <- .score_missing_fgsea_terms(
    ranks, term2gene, eligible, "TERM1", exponent = 1
  )

  expect_equal(nrow(fallback), 1L)
  expect_identical(fallback$ID, "TERM1")
  expect_equal(fallback$NES, 0)
  expect_equal(fallback$pvalue, 1)
  expect_equal(fallback$p.adjust, 1)
  expect_equal(fallback$qvalue, 1)
  expect_true(is.finite(fallback$enrichmentScore))
})

test_that("run_l1000_gsea is reproducible with rng_seed", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(
    GeneA = 3.5, GeneB = 2.8, GeneC = 1.9, GeneD = 0.7,
    GeneE = -0.5, GeneF = -1.4, GeneG = -2.6, GeneH = -3.2
  )
  term2gene <- data.frame(
    term = c("TERM_UP", "TERM_UP", "TERM_UP", "TERM_DN", "TERM_DN", "TERM_DN"),
    gene = c("GeneA", "GeneB", "GeneC", "GeneF", "GeneG", "GeneH"),
    stringsAsFactors = FALSE
  )

  first <- suppressWarnings(run_l1000_gsea(
    ranks,
    term2gene,
    min_gs_size = 2,
    pvalue_cutoff = 1,
    rng_seed = 104L
  ))
  stats::runif(10)
  second <- suppressWarnings(run_l1000_gsea(
    ranks,
    term2gene,
    min_gs_size = 2,
    pvalue_cutoff = 1,
    rng_seed = 104L
  ))

  expect_equal(first, second)
})
