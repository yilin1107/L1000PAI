test_that("detect_library_mode identifies paired and single libraries", {
  paired_t2g <- data.frame(
    term = c("TP53_Up", "TP53_Up", "TP53_Down", "TP53_Down"),
    gene = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )

  single_t2g <- data.frame(
    term = c("TP53", "TP53", "RB1", "RB1"),
    gene = c("A", "B", "C", "D"),
    stringsAsFactors = FALSE
  )

  expect_equal(detect_library_mode(paired_t2g), "paired")
  expect_equal(detect_library_mode(single_t2g), "single")
})

test_that("classify_single_signatures classifies positive and negative enrichment", {
  single_tbl <- data.frame(
    term = c("TP53", "RB1", "MYC"),
    NES = c(1.8, -1.6, 0.3),
    pvalue = c(0.01, 0.02, 0.5),
    p.adjust = c(0.02, 0.03, 0.5),
    stringsAsFactors = FALSE
  )

  res <- classify_single_signatures(
    single_table = single_tbl,
    nes_threshold = 1,
    p_value_threshold = 0.05
  )

  expect_equal(as.character(res$activity_class), c(
    "Positive_enrichment",
    "Negative_enrichment",
    "Not_significant"
  ))
})

test_that("run_l1000_pai supports single-signature library mode", {
  testthat::skip_if_not_installed("clusterProfiler")

  de_tbl <- data.frame(
    gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
    stat = c(3, 2, 1, -1, -2, -3),
    stringsAsFactors = FALSE
  )

  single_t2g <- data.frame(
    term = c("TP53", "TP53", "RB1", "RB1"),
    gene = c("GeneA", "GeneB", "GeneE", "GeneF"),
    stringsAsFactors = FALSE
  )

  res <- run_l1000_pai(
    x = de_tbl,
    gene_col = "gene",
    stat_col = "stat",
    term2gene = single_t2g,
    library_mode = "single",
    min_gs_size = 2,
    pvalue_cutoff = 1,
    seed = TRUE,
    single_nes_threshold = 0.5,
    classify_p_value_threshold = 1
  )

  expect_true(methods::is(res, "L1000PAIResult"))
  expect_equal(result_mode(res), "single")
  expect_s3_class(result_table(res), "data.frame")
  expect_true("activity_class" %in% colnames(result_table(res)))
})
