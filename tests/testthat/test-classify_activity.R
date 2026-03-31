test_that("classify_activity supports up_minus_down direction", {
  pair_table <- data.frame(
    perturbagen = "DrugA",
    NES_Up = 2.0,
    NES_Down = -1.5,
    activity_score = 3.5,
    p_combined = 0.01,
    padj_combined = 0.02,
    stringsAsFactors = FALSE
  )

  res <- classify_activity(
    pair_table,
    activity_threshold = 1,
    p_value_threshold = 0.05,
    score_direction = "up_minus_down"
  )

  expect_equal(as.character(res$activity_class), "Activated")
})
