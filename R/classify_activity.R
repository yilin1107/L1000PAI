#' Classify perturbation activity states
#'
#' Classify paired L1000 perturbation results into activity states based on
#' activity score, statistical significance, and directional consistency.
#'
#' @param pair_table A data.frame returned by [collapse_pairs()].
#' @param activity_threshold Minimum absolute activity score required for
#'   activation/inhibition calls. Default is `1`.
#' @param p_value_threshold Significance threshold applied to the selected
#'   p-value column. Default is `0.05`.
#' @param p_col Name of the p-value column to use. If `NULL`, the function
#'   will use `"padj_combined"` when available, otherwise `"p_combined"`.
#' @param score_direction One of `"down_minus_up"` or `"up_minus_down"`.
#'
#' @return A data.frame with added columns:
#'   `is_significant`, `direction_consistent`, and `activity_class`.
#' @export
classify_activity <- function(
    pair_table,
    activity_threshold = 1,
    p_value_threshold = 0.05,
    p_col = NULL,
    score_direction = c("down_minus_up", "up_minus_down")
) {
  score_direction <- match.arg(score_direction)

  if (!is.data.frame(pair_table)) {
    stop("`pair_table` must be a data.frame.")
  }

  required_cols <- c("NES_Up", "NES_Down", "activity_score")
  missing_cols <- setdiff(required_cols, colnames(pair_table))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in `pair_table`: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (is.null(p_col)) {
    if ("padj_combined" %in% colnames(pair_table)) {
      p_col <- "padj_combined"
    } else if ("p_combined" %in% colnames(pair_table)) {
      p_col <- "p_combined"
    } else {
      stop("No p-value column found. Please provide `p_col` explicitly.")
    }
  }

  if (!p_col %in% colnames(pair_table)) {
    stop("`p_col` not found in `pair_table`.")
  }

  if (!is.numeric(activity_threshold) || length(activity_threshold) != 1 ||
      is.na(activity_threshold) || activity_threshold < 0) {
    stop("`activity_threshold` must be a single non-negative numeric value.")
  }

  if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1 ||
      is.na(p_value_threshold) || p_value_threshold < 0 || p_value_threshold > 1) {
    stop("`p_value_threshold` must be a single numeric value between 0 and 1.")
  }

  x <- pair_table

  x$is_significant <- !is.na(x[[p_col]]) & x[[p_col]] < p_value_threshold

  if (score_direction == "down_minus_up") {
    activated_pattern <- !is.na(x$NES_Down) & !is.na(x$NES_Up) &
      x$NES_Down > 0 & x$NES_Up < 0

    inhibited_pattern <- !is.na(x$NES_Down) & !is.na(x$NES_Up) &
      x$NES_Down < 0 & x$NES_Up > 0
  } else {
    activated_pattern <- !is.na(x$NES_Down) & !is.na(x$NES_Up) &
      x$NES_Up > 0 & x$NES_Down < 0

    inhibited_pattern <- !is.na(x$NES_Down) & !is.na(x$NES_Up) &
      x$NES_Up < 0 & x$NES_Down > 0
  }

  x$direction_consistent <- activated_pattern | inhibited_pattern

  x$activity_class <- "Not_significant"
  x$activity_class[x$is_significant] <- "Discordant"

  idx_activated <- x$is_significant &
    activated_pattern &
    !is.na(x$activity_score) &
    x$activity_score >= activity_threshold

  x$activity_class[idx_activated] <- "Activated"

  idx_inhibited <- x$is_significant &
    inhibited_pattern &
    !is.na(x$activity_score) &
    x$activity_score <= -activity_threshold

  x$activity_class[idx_inhibited] <- "Inhibited"

  x$activity_class <- factor(
    x$activity_class,
    levels = c("Activated", "Inhibited", "Discordant", "Not_significant")
  )

  x
}
