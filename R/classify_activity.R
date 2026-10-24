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
#'   prefers `"padj_coherent"` when available, followed by
#'   `"padj_combined"` and `"p_combined"`.
#' @param score_col Name of the activity-score column used for the effect-size
#'   threshold. The default preserves the legacy `"activity_score"`. Use
#'   `"activity_gls_z"` after [score_joint_activity()] for a standardized GLS
#'   threshold while retaining a conjunction p-value in `p_col`.
#' @param require_coherent Logical; when `TRUE` (the default) and
#'   `padj_coherent` is available, that directional conjunction must also pass
#'   `p_value_threshold`. This prevents a global or legacy p-value supplied in
#'   `p_col` from replacing evidence that both signature arms are coherent.
#' @param score_direction One of `"down_minus_up"` or `"up_minus_down"`.
#'
#' @return A data.frame with added columns:
#'   `is_selected_significant`, `is_coherent_significant`, `is_significant`,
#'   `direction_consistent`, and `activity_class`.
#' @export
classify_activity <- function(
    pair_table,
    activity_threshold = 1,
    p_value_threshold = 0.05,
    p_col = NULL,
    score_direction = c("down_minus_up", "up_minus_down"),
    score_col = "activity_score",
    require_coherent = TRUE
) {
  score_direction <- match.arg(score_direction)

  if (!is.data.frame(pair_table)) {
    stop("`pair_table` must be a data.frame.")
  }

  if (!is.character(score_col) || length(score_col) != 1L || is.na(score_col)) {
    stop("`score_col` must be a single column name.")
  }

  required_cols <- c("NES_Up", "NES_Down", score_col)
  missing_cols <- setdiff(required_cols, colnames(pair_table))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in `pair_table`: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  if (!is.numeric(pair_table[[score_col]])) {
    stop("`score_col` must identify a numeric column.")
  }

  if (is.null(p_col)) {
    if ("padj_coherent" %in% colnames(pair_table)) {
      p_col <- "padj_coherent"
    } else if ("padj_combined" %in% colnames(pair_table)) {
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
  if (!is.numeric(pair_table[[p_col]])) {
    stop("`p_col` must identify a numeric column.")
  }

  if (!is.numeric(activity_threshold) || length(activity_threshold) != 1 ||
      is.na(activity_threshold) || activity_threshold < 0) {
    stop("`activity_threshold` must be a single non-negative numeric value.")
  }

  if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1 ||
      is.na(p_value_threshold) || p_value_threshold < 0 || p_value_threshold > 1) {
    stop("`p_value_threshold` must be a single numeric value between 0 and 1.")
  }
  if (!is.logical(require_coherent) || length(require_coherent) != 1L ||
      is.na(require_coherent)) {
    stop("`require_coherent` must be a single non-missing logical value.")
  }

  x <- pair_table

  x$is_selected_significant <- !is.na(x[[p_col]]) &
    x[[p_col]] < p_value_threshold

  if ("padj_coherent" %in% colnames(x)) {
    if (!is.numeric(x$padj_coherent)) {
      stop("`padj_coherent` must be numeric when present.")
    }
    x$is_coherent_significant <- !is.na(x$padj_coherent) &
      x$padj_coherent < p_value_threshold
  } else {
    x$is_coherent_significant <- rep(NA, nrow(x))
  }

  x$is_significant <- x$is_selected_significant
  if (isTRUE(require_coherent) && "padj_coherent" %in% colnames(x)) {
    x$is_significant <- x$is_significant & x$is_coherent_significant
  }
  both_arms_p_col <- if ("padj_both_arms" %in% colnames(x)) {
    "padj_both_arms"
  } else {
    p_col
  }
  x$is_both_arms_significant <- !is.na(x[[both_arms_p_col]]) &
    x[[both_arms_p_col]] < p_value_threshold

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
  observed_arms <- !is.na(x$NES_Up) & !is.na(x$NES_Down)
  x$activity_class[x$is_both_arms_significant & observed_arms &
                    !x$direction_consistent] <- "Discordant"

  idx_activated <- x$is_significant &
    activated_pattern &
    !is.na(x[[score_col]]) &
    x[[score_col]] >= activity_threshold

  x$activity_class[idx_activated] <- "Activated"

  idx_inhibited <- x$is_significant &
    inhibited_pattern &
    !is.na(x[[score_col]]) &
    x[[score_col]] <= -activity_threshold

  x$activity_class[idx_inhibited] <- "Inhibited"

  x$activity_class <- factor(
    x$activity_class,
    levels = c("Activated", "Inhibited", "Discordant", "Not_significant")
  )

  x
}
