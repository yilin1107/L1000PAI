test_that("collapse_pairs supports up_minus_down scoring", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(2.0, -1.5),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(
    gsea_table,
    combine_p_method = "max",
    score_direction = "up_minus_down"
  )

  expect_equal(res$activity_score, 2.0 - (-1.5))
})
