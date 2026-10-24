test_that("plot_l1000_pai creates an activity plot", {
  testthat::skip_if_not_installed("ggplot2")

  df <- data.frame(
    term = c("A", "B", "C"),
    perturbagen = c("A", "B", "C"),
    activity_score = c(2, -2, 0.2),
    padj_combined = c(0.01, 0.02, 0.8),
    activity_class = c("Activated", "Inhibited", "Not_significant"),
    stringsAsFactors = FALSE
  )

  p <- plot_l1000_pai(df, type = "activity", library = NULL, label_top = 0)
  expect_s3_class(p, "ggplot")

  positional <- plot_l1000_pai(
    df, "activity", NULL, NULL, NULL, 0.05, 1, 200, 0, NULL
  )
  expect_s3_class(positional, "ggplot")
})

test_that("plot_l1000_pai accepts a GLS score and coherent p-value", {
  testthat::skip_if_not_installed("ggplot2")

  df <- data.frame(
    term = c("A", "B"),
    perturbagen = c("A", "B"),
    activity_score = c(3, -2),
    activity_gls_z = c(1.5, -1.2),
    padj_coherent = c(0.02, 0.04),
    activity_class = c("Activated", "Inhibited"),
    stringsAsFactors = FALSE
  )

  p <- plot_l1000_pai(
    df,
    type = "activity",
    library = NULL,
    score_col = "activity_gls_z",
    label_top = 0
  )

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$x, "activity_gls_z")
  expect_match(p$labels$y, "padj_coherent", fixed = TRUE)
})

test_that("plot_l1000_pai creates a chem category plot", {
  testthat::skip_if_not_installed("ggplot2")

  df <- data.frame(
    term = c("Erlotinib", "Simvastatin"),
    perturbagen = c("Erlotinib", "Simvastatin"),
    activity_score = c(2, 3),
    p_combined = c(0.01, 0.001),
    padj_combined = c(0.02, 0.002),
    is_significant = c(TRUE, TRUE),
    activity_class = c("Activated", "Activated"),
    stringsAsFactors = FALSE
  )

  p <- plot_l1000_pai(df, type = "category", library = "chem", top_n = 5)
  expect_s3_class(p, "ggplot")
})
