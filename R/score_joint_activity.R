#' Score paired activity with a synchronized joint null
#'
#' Estimate a covariance-aware generalized least-squares (GLS) activity score
#' for paired Up/Down signatures. The null table must contain synchronized null
#' draws for each perturbagen so that the within-pair covariance is estimated
#' without treating observed biological signals as null observations.
#'
#' This function adds latent/global activity evidence. It does not replace the
#' directional conjunction p-values produced by [collapse_pairs()], which are
#' required when claiming that both arms support a coherent activity state.
#'
#' @param pair_table A paired result data.frame returned by [collapse_pairs()].
#' @param joint_null A data.frame containing synchronized null draws. By
#'   default it must contain `perturbagen`, `permutation_id`, `NES_Up`, and
#'   `NES_Down`.
#' @param score_direction One of `"down_minus_up"` or `"up_minus_down"`.
#' @param perturbagen_col Name of the perturbagen column in both tables.
#' @param draw_col Name of the synchronized null-draw identifier column.
#' @param up_col Name of the Up-score column in `joint_null`.
#' @param down_col Name of the Down-score column in `joint_null`.
#' @param ridge Non-negative relative ridge added to each pair covariance
#'   matrix. The added diagonal value is `ridge * mean(diag(covariance))`.
#'   The default, zero, requires a positive-definite empirical covariance.
#' @param min_null_draws Minimum number of complete synchronized null draws per
#'   perturbagen. Default is `50`.
#' @param p_adjust_method Method passed to [stats::p.adjust()] for the empirical
#'   GLS p-values. Default is `"BH"`.
#'
#' @return `pair_table` with appended GLS estimates, empirical and asymptotic
#'   global p-values, heterogeneity diagnostics, and null-quality metadata.
#'   Existing columns and row order are preserved.
#' @export
score_joint_activity <- function(
    pair_table,
    joint_null,
    score_direction = c("down_minus_up", "up_minus_down"),
    perturbagen_col = "perturbagen",
    draw_col = "permutation_id",
    up_col = "NES_Up",
    down_col = "NES_Down",
    ridge = 0,
    min_null_draws = 50L,
    p_adjust_method = "BH"
) {
  score_direction <- match.arg(score_direction)

  if (!is.data.frame(pair_table)) {
    stop("`pair_table` must be a data.frame.")
  }
  if (!is.data.frame(joint_null)) {
    stop("`joint_null` must be a data.frame of synchronized null draws.")
  }

  pair_required <- c(perturbagen_col, "NES_Up", "NES_Down")
  null_required <- c(perturbagen_col, draw_col, up_col, down_col)
  .require_joint_columns(pair_table, pair_required, "pair_table")
  .require_joint_columns(joint_null, null_required, "joint_null")

  if (!is.numeric(ridge) || length(ridge) != 1L || !is.finite(ridge) || ridge < 0) {
    stop("`ridge` must be a single non-negative finite number.")
  }
  if (!is.numeric(min_null_draws) || length(min_null_draws) != 1L ||
      !is.finite(min_null_draws) || min_null_draws < 3 ||
      min_null_draws != as.integer(min_null_draws)) {
    stop("`min_null_draws` must be a single integer of at least 3.")
  }
  min_null_draws <- as.integer(min_null_draws)

  valid_adjust_methods <- stats::p.adjust.methods
  if (!is.character(p_adjust_method) || length(p_adjust_method) != 1L ||
      !p_adjust_method %in% valid_adjust_methods) {
    stop(
      "`p_adjust_method` must be one of: ",
      paste(valid_adjust_methods, collapse = ", "), "."
    )
  }

  pair_ids <- as.character(pair_table[[perturbagen_col]])
  null_ids <- as.character(joint_null[[perturbagen_col]])
  draw_ids <- as.character(joint_null[[draw_col]])

  if (anyNA(pair_ids) || any(!nzchar(pair_ids))) {
    stop("`pair_table` contains missing or empty perturbagen identifiers.")
  }
  if (anyDuplicated(pair_ids)) {
    stop("`pair_table` must contain at most one row per perturbagen.")
  }
  if (anyNA(null_ids) || any(!nzchar(null_ids)) || anyNA(draw_ids) || any(!nzchar(draw_ids))) {
    stop("`joint_null` contains missing or empty perturbagen/draw identifiers.")
  }

  duplicate_key <- duplicated(data.frame(
    perturbagen = null_ids,
    draw = draw_ids,
    stringsAsFactors = FALSE
  ))
  if (any(duplicate_key)) {
    stop("`joint_null` contains duplicate perturbagen/permutation identifiers.")
  }

  pair_up <- .finite_numeric(pair_table$NES_Up, "pair_table$NES_Up")
  pair_down <- .finite_numeric(pair_table$NES_Down, "pair_table$NES_Down")
  null_up <- .finite_numeric(joint_null[[up_col]], paste0("joint_null$", up_col))
  null_down <- .finite_numeric(joint_null[[down_col]], paste0("joint_null$", down_col))

  missing_ids <- setdiff(pair_ids, unique(null_ids))
  if (length(missing_ids) > 0L) {
    stop(
      "`joint_null` has no draws for perturbagen(s): ",
      paste(utils::head(missing_ids, 10L), collapse = ", "),
      if (length(missing_ids) > 10L) ", ..." else "."
    )
  }

  n_pairs <- nrow(pair_table)
  activity_gls <- activity_gls_se <- activity_gls_z <- rep(NA_real_, n_pairs)
  p_gls_global <- p_gls_asymptotic <- rep(NA_real_, n_pairs)
  heterogeneity_Q <- p_heterogeneity_asymptotic <- rep(NA_real_, n_pairs)
  null_correlation <- condition_number <- rep(NA_real_, n_pairs)
  null_center_arm1 <- null_center_arm2 <- rep(NA_real_, n_pairs)
  null_variance_arm1 <- null_variance_arm2 <- null_covariance_arm12 <- rep(NA_real_, n_pairs)
  gls_weight_arm1 <- gls_weight_arm2 <- rep(NA_real_, n_pairs)
  n_null_draws <- integer(n_pairs)
  near_singular_ids <- character()

  observed_oriented <- .orient_pair_scores(
    up = pair_up,
    down = pair_down,
    score_direction = score_direction
  )

  for (i in seq_len(n_pairs)) {
    idx <- null_ids == pair_ids[[i]]
    oriented_null <- .orient_pair_scores(
      up = null_up[idx],
      down = null_down[idx],
      score_direction = score_direction
    )

    n_draws <- nrow(oriented_null)
    n_null_draws[[i]] <- n_draws
    if (n_draws < min_null_draws) {
      stop(
        "Perturbagen `", pair_ids[[i]], "` has ", n_draws,
        " complete null draws; at least ", min_null_draws, " are required."
      )
    }

    null_center <- colMeans(oriented_null)
    centered_null <- sweep(oriented_null, 2L, null_center, FUN = "-")
    covariance <- stats::cov(centered_null)

    if (any(!is.finite(covariance)) || any(diag(covariance) <= 0)) {
      stop(
        "Perturbagen `", pair_ids[[i]],
        "` has a non-finite or zero-variance null covariance."
      )
    }

    ridge_scale <- mean(diag(covariance))
    covariance_regularized <- covariance + diag(ridge * ridge_scale, nrow = 2L)
    eigenvalues <- eigen(covariance_regularized, symmetric = TRUE, only.values = TRUE)$values

    if (min(eigenvalues) <= 0) {
      stop(
        "Perturbagen `", pair_ids[[i]],
        "` has a singular null covariance; provide more null draws or set `ridge` explicitly."
      )
    }

    condition_number[[i]] <- max(eigenvalues) / min(eigenvalues)
    if (condition_number[[i]] > 1e8) {
      near_singular_ids <- c(near_singular_ids, pair_ids[[i]])
    }

    covariance_inverse <- solve(covariance_regularized)
    loading <- c(1, 1)
    precision_sum <- as.numeric(crossprod(loading, covariance_inverse %*% loading))
    weights <- as.numeric(covariance_inverse %*% loading) / precision_sum
    standard_error <- 1 / sqrt(precision_sum)

    observed_centered <- observed_oriented[i, ] - null_center
    theta <- sum(weights * observed_centered)
    z_value <- theta / standard_error
    null_theta <- as.numeric(centered_null %*% weights)

    residual <- observed_centered - loading * theta
    q_value <- as.numeric(crossprod(residual, covariance_inverse %*% residual))

    activity_gls[[i]] <- theta
    activity_gls_se[[i]] <- standard_error
    activity_gls_z[[i]] <- z_value
    p_gls_global[[i]] <- (1 + sum(abs(null_theta) >= abs(theta))) / (n_draws + 1)
    p_gls_asymptotic[[i]] <- 2 * stats::pnorm(-abs(z_value))
    heterogeneity_Q[[i]] <- q_value
    p_heterogeneity_asymptotic[[i]] <- stats::pchisq(
      q_value,
      df = 1,
      lower.tail = FALSE
    )
    null_correlation[[i]] <- stats::cov2cor(covariance)[1, 2]
    null_center_arm1[[i]] <- null_center[[1]]
    null_center_arm2[[i]] <- null_center[[2]]
    null_variance_arm1[[i]] <- covariance[1, 1]
    null_variance_arm2[[i]] <- covariance[2, 2]
    null_covariance_arm12[[i]] <- covariance[1, 2]
    gls_weight_arm1[[i]] <- weights[[1]]
    gls_weight_arm2[[i]] <- weights[[2]]
  }

  if (length(near_singular_ids) > 0L) {
    warning(
      "Near-singular null covariance (condition number > 1e8) for perturbagen(s): ",
      paste(unique(near_singular_ids), collapse = ", "),
      ". Consider increasing null draws or setting `ridge`.",
      call. = FALSE
    )
  }

  out <- pair_table
  out$activity_gls <- activity_gls
  out$activity_gls_se <- activity_gls_se
  out$activity_gls_z <- activity_gls_z
  out$p_gls_global <- p_gls_global
  out$padj_gls_global <- .adjust_complete_p(p_gls_global, method = p_adjust_method)
  out$p_gls_asymptotic <- p_gls_asymptotic
  out$heterogeneity_Q <- heterogeneity_Q
  out$p_heterogeneity_asymptotic <- p_heterogeneity_asymptotic
  out$null_correlation <- null_correlation
  out$null_center_arm1 <- null_center_arm1
  out$null_center_arm2 <- null_center_arm2
  out$null_variance_arm1 <- null_variance_arm1
  out$null_variance_arm2 <- null_variance_arm2
  out$null_covariance_arm12 <- null_covariance_arm12
  out$gls_weight_arm1 <- gls_weight_arm1
  out$gls_weight_arm2 <- gls_weight_arm2
  out$n_null_draws <- n_null_draws
  out$gls_condition_number <- condition_number
  out$gls_ridge <- rep(ridge, n_pairs)
  out$gls_covariance_source <- rep("synchronized_joint_null", n_pairs)

  attr(out, "joint_model") <- list(
    method = "pair_specific_empirical_gls",
    score_direction = score_direction,
    p_value = "two_sided_empirical_with_finite_sample_correction",
    p_adjust_method = p_adjust_method,
    ridge = ridge,
    min_null_draws = min_null_draws,
    covariance_source = "synchronized_joint_null"
  )

  out
}

.require_joint_columns <- function(x, required, object_name) {
  missing_cols <- setdiff(required, colnames(x))
  if (length(missing_cols) > 0L) {
    stop(
      "Missing required columns in `", object_name, "`: ",
      paste(missing_cols, collapse = ", "), "."
    )
  }
}

.finite_numeric <- function(x, label) {
  if (!is.numeric(x)) {
    stop("`", label, "` must be numeric.")
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`", label, "` must contain only finite, non-missing values.")
  }
  as.numeric(x)
}

.orient_pair_scores <- function(up, down, score_direction) {
  if (score_direction == "down_minus_up") {
    out <- cbind(arm1 = -up, arm2 = down)
  } else {
    out <- cbind(arm1 = up, arm2 = -down)
  }
  storage.mode(out) <- "double"
  out
}

.adjust_complete_p <- function(p, method = "BH") {
  out <- rep(NA_real_, length(p))
  keep <- !is.na(p)
  if (any(keep)) {
    out[keep] <- stats::p.adjust(p[keep], method = method)
  }
  out
}
