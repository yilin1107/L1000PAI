test_that("run_l1000_pai works with a named rank vector and user term2gene", {
  testthat::skip_if_not_installed("clusterProfiler")

  ranks <- c(
    DrugA_Target1 = 3.5,
    DrugA_Target2 = 2.8,
    DrugA_Target3 = 2.0,
    DrugB_Target1 = -2.5,
    DrugB_Target2 = -2.0,
    DrugB_Target3 = -1.8,
    Noise1 = 0.4,
    Noise2 = -0.3
  )
  names(ranks) <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH")

  term2gene <- data.frame(
    term = c(
      "DrugA_Up", "DrugA_Up", "DrugA_Up",
      "DrugA_Down", "DrugA_Down", "DrugA_Down",
      "DrugB_Up", "DrugB_Up", "DrugB_Up",
      "DrugB_Down", "DrugB_Down", "DrugB_Down"
    ),
    gene = c(
      "GeneD", "GeneE", "GeneF",
      "GeneA", "GeneB", "GeneC",
      "GeneA", "GeneB", "GeneC",
      "GeneD", "GeneE", "GeneF"
    ),
    stringsAsFactors = FALSE
  )

  res <- run_l1000_pai(
    x = ranks,
    term2gene = term2gene,
    min_gs_size = 2,
    max_gs_size = 500,
    pvalue_cutoff = 1,
    seed = TRUE,
    combine_p_method = "max",
    activity_threshold = 0.5,
    classify_p_value_threshold = 0.25
  )

  expect_true(methods::is(res, "L1000PAIResult"))
  expect_s3_class(gsea_table(res), "data.frame")
  expect_s3_class(pair_table(res), "data.frame")
  expect_type(rank_table(res), "double")
  expect_type(pai_metadata(res), "list")
})

test_that("run_l1000_pai works with data.frame input", {
  testthat::skip_if_not_installed("clusterProfiler")

  de_tbl <- data.frame(
    gene = c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF"),
    stat = c(3, 2, 1, -1, -2, -3),
    stringsAsFactors = FALSE
  )

  term2gene <- data.frame(
    term = c(
      "DrugA_Up", "DrugA_Up",
      "DrugA_Down", "DrugA_Down"
    ),
    gene = c(
      "GeneE", "GeneF",
      "GeneA", "GeneB"
    ),
    stringsAsFactors = FALSE
  )

  res <- run_l1000_pai(
    x = de_tbl,
    gene_col = "gene",
    stat_col = "stat",
    term2gene = term2gene,
    min_gs_size = 2,
    pvalue_cutoff = 1,
    seed = TRUE
  )

  expect_true(methods::is(res, "L1000PAIResult"))
  expect_equal(pai_metadata(res)$input_type, "data.frame")
})

test_that("run_l1000_pai errors when neither term2gene nor library is provided", {
  ranks <- c(GeneA = 2, GeneB = 1, GeneC = -1)

  expect_error(
    run_l1000_pai(x = ranks),
    "Either `term2gene` or `library` must be provided"
  )
})

test_that("run_l1000_pai can use bundled library if L1000PAIdata is installed", {
  testthat::skip_if_not_installed("clusterProfiler")
  testthat::skip_if_not_installed("L1000PAIdata")

  lib_tbl <- get_l1000_library("crispr")

  genes <- unique(lib_tbl$gene)
  genes <- genes[!is.na(genes) & genes != ""]
  testthat::skip_if(length(genes) < 50, "Bundled library has too few genes for test.")

  sel_genes <- genes[1:50]
  ranks <- seq(from = 5, to = -5, length.out = length(sel_genes))
  names(ranks) <- sel_genes

  res <- run_l1000_pai(
    x = ranks,
    library = "crispr",
    min_gs_size = 2,
    max_gs_size = 500,
    pvalue_cutoff = 1,
    seed = TRUE,
    classify_p_value_threshold = 1
  )

  expect_true(methods::is(res, "L1000PAIResult"))
  expect_s3_class(gsea_table(res), "data.frame")
})
