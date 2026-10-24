.make_joint_null <- function(
    perturbagen,
    covariance,
    score_direction = "down_minus_up"
) {
  base <- sqrt(3 / 2) * rbind(
    c(1, 0),
    c(-1, 0),
    c(0, 1),
    c(0, -1)
  )
  oriented <- base %*% chol(covariance)

  if (score_direction == "down_minus_up") {
    up <- -oriented[, 1]
    down <- oriented[, 2]
  } else {
    up <- oriented[, 1]
    down <- -oriented[, 2]
  }

  data.frame(
    perturbagen = perturbagen,
    permutation_id = seq_len(nrow(oriented)),
    NES_Up = up,
    NES_Down = down,
    stringsAsFactors = FALSE
  )
}

test_that("score_joint_activity reproduces identity-covariance GLS", {
  pair_table <- data.frame(
    perturbagen = "DrugA",
    NES_Up = -2,
    NES_Down = 1,
    activity_score = 3
  )
  joint_null <- .make_joint_null("DrugA", diag(2))

  res <- L1000PAI:::score_joint_activity(
    pair_table,
    joint_null,
    min_null_draws = 4
  )

  expect_equal(res$activity_gls, 1.5, tolerance = 1e-12)
  expect_equal(res$activity_gls_se, sqrt(0.5), tolerance = 1e-12)
  expect_equal(res$activity_gls_z, 2.121320343559642, tolerance = 1e-12)
  expect_equal(res$p_gls_global, 0.2)
  expect_equal(res$padj_gls_global, 0.2)
  expect_equal(res$p_gls_asymptotic, 0.03389485352469, tolerance = 1e-12)
  expect_equal(res$heterogeneity_Q, 0.5, tolerance = 1e-12)
  expect_equal(
    res$p_heterogeneity_asymptotic,
    0.4795001221869535,
    tolerance = 1e-14
  )
  expect_equal(res$null_correlation, 0, tolerance = 1e-12)
  expect_equal(res$null_center_arm1, 0, tolerance = 1e-12)
  expect_equal(res$null_center_arm2, 0, tolerance = 1e-12)
  expect_equal(res$null_variance_arm1, 1, tolerance = 1e-12)
  expect_equal(res$null_variance_arm2, 1, tolerance = 1e-12)
  expect_equal(res$null_covariance_arm12, 0, tolerance = 1e-12)
  expect_equal(res$gls_weight_arm1, 0.5, tolerance = 1e-12)
  expect_equal(res$gls_weight_arm2, 0.5, tolerance = 1e-12)
  expect_equal(res$n_null_draws, 4L)
  expect_equal(res$gls_condition_number, 1, tolerance = 1e-12)
  expect_equal(res$gls_ridge, 0)
  expect_equal(res$gls_covariance_source, "synchronized_joint_null")

  model <- attr(res, "joint_model")
  expect_equal(model$method, "pair_specific_empirical_gls")
  expect_equal(model$score_direction, "down_minus_up")
  expect_equal(model$p_value, "two_sided_empirical_with_finite_sample_correction")
})

test_that("score_joint_activity accounts for positive pair correlation", {
  covariance <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  pair_table <- data.frame(
    perturbagen = "DrugA",
    NES_Up = -2,
    NES_Down = 1
  )

  res <- L1000PAI:::score_joint_activity(
    pair_table,
    .make_joint_null("DrugA", covariance),
    min_null_draws = 4
  )

  expect_equal(res$activity_gls, 1.5, tolerance = 1e-12)
  expect_equal(res$activity_gls_se, sqrt(0.75), tolerance = 1e-12)
  expect_equal(res$activity_gls_z, sqrt(3), tolerance = 1e-12)
  expect_equal(res$p_gls_asymptotic, 0.08326451666355, tolerance = 1e-14)
  expect_equal(res$heterogeneity_Q, 1, tolerance = 1e-12)
  expect_equal(
    res$p_heterogeneity_asymptotic,
    0.3173105078629141,
    tolerance = 1e-14
  )
  expect_equal(res$null_correlation, 0.5, tolerance = 1e-12)
  expect_equal(res$gls_condition_number, 3, tolerance = 1e-12)
})

test_that("score_joint_activity weights unequal arm variances", {
  covariance <- diag(c(1, 4))
  pair_table <- data.frame(
    perturbagen = "DrugA",
    NES_Up = -2,
    NES_Down = 1
  )

  res <- L1000PAI:::score_joint_activity(
    pair_table,
    .make_joint_null("DrugA", covariance),
    min_null_draws = 4
  )

  expect_equal(res$activity_gls, 1.8, tolerance = 1e-12)
  expect_equal(res$activity_gls_se, sqrt(0.8), tolerance = 1e-12)
  expect_equal(res$activity_gls_z, 2.012461179749811, tolerance = 1e-12)
  expect_equal(res$p_gls_asymptotic, 0.04417134490844, tolerance = 1e-12)
  expect_equal(res$heterogeneity_Q, 0.2, tolerance = 1e-12)
  expect_equal(res$null_correlation, 0, tolerance = 1e-12)
  expect_equal(res$gls_condition_number, 4, tolerance = 1e-12)
  expect_equal(res$gls_weight_arm1, 0.8, tolerance = 1e-12)
  expect_equal(res$gls_weight_arm2, 0.2, tolerance = 1e-12)
})

test_that("score direction changes only the input orientation", {
  covariance <- diag(2)
  down_pair <- data.frame(perturbagen = "DrugA", NES_Up = -2, NES_Down = 1)
  up_pair <- data.frame(perturbagen = "DrugA", NES_Up = 2, NES_Down = -1)

  down_res <- L1000PAI:::score_joint_activity(
    down_pair,
    .make_joint_null("DrugA", covariance, "down_minus_up"),
    score_direction = "down_minus_up",
    min_null_draws = 4
  )
  up_res <- L1000PAI:::score_joint_activity(
    up_pair,
    .make_joint_null("DrugA", covariance, "up_minus_down"),
    score_direction = "up_minus_down",
    min_null_draws = 4
  )

  expect_equal(up_res$activity_gls, down_res$activity_gls)
  expect_equal(up_res$activity_gls_se, down_res$activity_gls_se)
  expect_equal(up_res$activity_gls_z, down_res$activity_gls_z)
  expect_equal(up_res$p_gls_global, down_res$p_gls_global)
})

test_that("empirical GLS p-values use finite-sample correction and count ties", {
  theta <- sqrt(3 / 2) / 2
  pair_table <- data.frame(
    perturbagen = "DrugA",
    NES_Up = -theta,
    NES_Down = theta
  )

  res <- L1000PAI:::score_joint_activity(
    pair_table,
    .make_joint_null("DrugA", diag(2)),
    min_null_draws = 4
  )

  expect_equal(res$activity_gls, theta, tolerance = 1e-12)
  expect_equal(res$p_gls_global, 1)
})

test_that("score_joint_activity preserves pair order and adjusts empirical p-values", {
  pair_table <- data.frame(
    perturbagen = c("DrugB", "DrugA"),
    NES_Up = c(-2, 0),
    NES_Down = c(1, 0)
  )
  joint_null <- rbind(
    .make_joint_null("DrugA", diag(2)),
    .make_joint_null("DrugB", diag(2))
  )

  res <- L1000PAI:::score_joint_activity(
    pair_table,
    joint_null,
    min_null_draws = 4
  )

  expect_equal(res$perturbagen, c("DrugB", "DrugA"))
  expect_equal(res$p_gls_global, c(0.2, 1))
  expect_equal(res$padj_gls_global, c(0.4, 1))
})

test_that("a singular covariance requires an explicit ridge", {
  t <- c(-1, 0, 1)
  oriented <- cbind(t, t)
  joint_null <- data.frame(
    perturbagen = "DrugA",
    permutation_id = seq_along(t),
    NES_Up = -oriented[, 1],
    NES_Down = oriented[, 2]
  )
  pair_table <- data.frame(perturbagen = "DrugA", NES_Up = -2, NES_Down = 1)

  expect_error(
    L1000PAI:::score_joint_activity(
      pair_table,
      joint_null,
      min_null_draws = 3
    ),
    "singular null covariance"
  )

  regularized <- L1000PAI:::score_joint_activity(
    pair_table,
    joint_null,
    min_null_draws = 3,
    ridge = 0.1
  )
  expect_true(all(is.finite(regularized$activity_gls_z)))
  expect_equal(regularized$gls_ridge, 0.1)
  expect_equal(regularized$gls_condition_number, 21, tolerance = 1e-10)
})

test_that("score_joint_activity validates synchronized null inputs", {
  pair_table <- data.frame(perturbagen = "DrugA", NES_Up = -2, NES_Down = 1)
  joint_null <- .make_joint_null("DrugA", diag(2))

  expect_error(
    L1000PAI:::score_joint_activity(
      pair_table,
      joint_null,
      min_null_draws = 5
    ),
    "at least 5"
  )

  duplicate <- rbind(joint_null, joint_null[1, ])
  expect_error(
    L1000PAI:::score_joint_activity(pair_table, duplicate, min_null_draws = 4),
    "duplicate perturbagen/permutation"
  )

  missing_id <- joint_null
  missing_id$perturbagen <- "Other"
  expect_error(
    L1000PAI:::score_joint_activity(pair_table, missing_id, min_null_draws = 4),
    "no draws for perturbagen"
  )

  zero_variance <- joint_null
  zero_variance$NES_Down <- 0
  expect_error(
    L1000PAI:::score_joint_activity(pair_table, zero_variance, min_null_draws = 4),
    "zero-variance null covariance"
  )

  nonfinite_pair <- pair_table
  nonfinite_pair$NES_Up <- Inf
  expect_error(
    L1000PAI:::score_joint_activity(nonfinite_pair, joint_null, min_null_draws = 4),
    "finite, non-missing"
  )

  expect_error(
    L1000PAI:::score_joint_activity(pair_table, joint_null, ridge = -1, min_null_draws = 4),
    "ridge"
  )
  expect_error(
    L1000PAI:::score_joint_activity(pair_table, joint_null, min_null_draws = 2),
    "min_null_draws"
  )
  expect_error(
    L1000PAI:::score_joint_activity(
      pair_table,
      joint_null,
      min_null_draws = 4,
      p_adjust_method = "bad"
    ),
    "p_adjust_method"
  )
})
