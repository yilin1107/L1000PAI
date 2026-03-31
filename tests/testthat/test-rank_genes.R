test_that("rank_genes.data.frame returns named numeric vector", {
  df <- data.frame(
    gene = c("GeneA", "GeneB", "GeneC"),
    stat = c(2.5, -1.2, 0.8)
  )

  res <- rank_genes(df)

  expect_type(res, "double")
  expect_equal(length(res), 3)
  expect_true(!is.null(names(res)))
  expect_equal(names(res), c("GeneA", "GeneC", "GeneB"))
  expect_equal(unname(res), c(2.5, 0.8, -1.2))
})

test_that("rank_genes.data.frame removes NA rows", {
  df <- data.frame(
    gene = c("GeneA", "GeneB", NA, "GeneD"),
    stat = c(2.5, NA, 1.1, -0.3)
  )

  res <- rank_genes(df)

  expect_equal(names(res), c("GeneA", "GeneD"))
  expect_equal(unname(res), c(2.5, -0.3))
})

test_that("rank_genes.data.frame keeps largest absolute stat for duplicated genes", {
  df <- data.frame(
    gene = c("GeneA", "GeneA", "GeneB"),
    stat = c(1.0, -3.0, 2.0)
  )

  res <- rank_genes(df)

  expect_equal(names(res), c("GeneB", "GeneA"))
  expect_equal(unname(res), c(2.0, -3.0))
})

test_that("rank_genes.data.frame errors when required columns are missing", {
  df <- data.frame(a = 1:3, b = 4:6)

  expect_error(rank_genes(df), "gene_col")
  expect_error(rank_genes(df, gene_col = "a"), "stat_col")
})
