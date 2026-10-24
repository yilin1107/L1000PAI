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
  expect_true(res$is_significant)
  expect_true(res$is_both_arms_significant)

  positional <- classify_activity(pair_table, 1, 0.05, NULL, "up_minus_down")
  expect_equal(as.character(positional$activity_class), "Activated")
})

test_that("classify_activity prefers coherent adjusted p-values by default", {
  pair_table <- data.frame(
    perturbagen = "DrugA",
    NES_Up = -2,
    NES_Down = 2,
    activity_score = 4,
    p_combined = 0.001,
    padj_combined = 0.001,
    padj_both_arms = 0.001,
    padj_coherent = 0.20,
    stringsAsFactors = FALSE
  )

  strict <- classify_activity(pair_table, activity_threshold = 1)
  legacy <- classify_activity(
    pair_table,
    activity_threshold = 1,
    p_col = "padj_combined",
    require_coherent = FALSE
  )

  expect_false(strict$is_significant)
  expect_true(strict$is_both_arms_significant)
  expect_equal(as.character(strict$activity_class), "Not_significant")
  expect_true(legacy$is_selected_significant)
  expect_false(legacy$is_coherent_significant)
  expect_true(legacy$is_significant)
  expect_equal(as.character(legacy$activity_class), "Activated")
})

test_that("classify_activity can threshold on a GLS score while retaining coherent p-values", {
  pair_table <- data.frame(
    perturbagen = c("WeakA", "StrongI"),
    NES_Up = c(-2, 2),
    NES_Down = c(2, -2),
    activity_score = c(4, -4),
    activity_gls_z = c(0.5, -2),
    padj_both_arms = c(0.01, 0.01),
    padj_coherent = c(0.01, 0.01),
    stringsAsFactors = FALSE
  )

  res <- classify_activity(
    pair_table,
    activity_threshold = 1,
    score_col = "activity_gls_z",
    score_direction = "down_minus_up"
  )

  expect_equal(
    as.character(res$activity_class),
    c("Not_significant", "Inhibited")
  )
  expect_true(all(res$is_significant))
  expect_true(all(res$is_both_arms_significant))
})

test_that("discordance requires both-arm evidence and is not latent significance", {
  pair_table <- data.frame(
    perturbagen = c("Disc", "WeakDisc"),
    NES_Up = c(1, 1),
    NES_Down = c(1, 1),
    activity_score = c(0, 0),
    padj_combined = c(0.01, 0.01),
    padj_both_arms = c(0.01, 0.20),
    padj_coherent = c(1, 1),
    stringsAsFactors = FALSE
  )

  res <- classify_activity(pair_table, activity_threshold = 0)

  expect_equal(
    as.character(res$activity_class),
    c("Discordant", "Not_significant")
  )
  expect_false(any(res$is_significant))
  expect_equal(res$is_both_arms_significant, c(TRUE, FALSE))
  expect_false(any(res$direction_consistent))
})

test_that("strong latent GLS evidence cannot replace the coherent gate", {
  pair_table <- data.frame(
    perturbagen = "SingleArmDominant",
    NES_Up = -10,
    NES_Down = 1,
    activity_score = 11,
    activity_gls_z = 6.36,
    padj_gls_global = 1e-8,
    padj_both_arms = 1,
    padj_coherent = 1,
    stringsAsFactors = FALSE
  )

  res <- classify_activity(
    pair_table,
    p_col = "padj_gls_global",
    score_col = "activity_gls_z",
    activity_threshold = 1,
    p_value_threshold = 0.05,
    score_direction = "down_minus_up"
  )

  expect_true(res$is_selected_significant)
  expect_false(res$is_coherent_significant)
  expect_false(res$is_significant)
  expect_true(res$direction_consistent)
  expect_equal(as.character(res$activity_class), "Not_significant")
})

test_that("classification significance threshold is strict", {
  pair_table <- data.frame(
    NES_Up = -2,
    NES_Down = 2,
    activity_score = 4,
    padj_coherent = 0.05,
    padj_both_arms = 0.05
  )

  res <- classify_activity(
    pair_table,
    activity_threshold = 1,
    p_value_threshold = 0.05
  )

  expect_false(res$is_significant)
  expect_false(res$is_both_arms_significant)
  expect_equal(as.character(res$activity_class), "Not_significant")
})

test_that("classify_activity validates score and p-value columns", {
  pair_table <- data.frame(
    NES_Up = -2,
    NES_Down = 2,
    activity_score = 4,
    padj_coherent = 0.01
  )

  expect_error(
    classify_activity(pair_table, score_col = "missing"),
    "Missing required columns"
  )
  bad_score <- pair_table
  bad_score$activity_score <- "four"
  expect_error(classify_activity(bad_score), "score_col.*numeric")
  expect_error(
    classify_activity(pair_table, p_col = "missing"),
    "p_col.*not found"
  )
  bad_p <- pair_table
  bad_p$padj_coherent <- "small"
  expect_error(classify_activity(bad_p), "p_col.*numeric")
  expect_error(classify_activity(pair_table, activity_threshold = -1), "activity_threshold")
  expect_error(classify_activity(pair_table, p_value_threshold = 1.1), "p_value_threshold")
})
