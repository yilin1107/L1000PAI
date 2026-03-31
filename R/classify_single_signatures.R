#' Classify single-signature enrichment states
#'
#' @param single_table A data.frame returned by [summarize_single_signatures()].
#' @param nes_threshold Minimum absolute NES required for calling positive or
#'   negative enrichment. Default is `1`.
#' @param p_value_threshold Significance threshold applied to the selected
#'   p-value column. Default is `0.05`.
#' @param p_col Name of the p-value column to use. If `NULL`, the function
#'   will use `"p.adjust"` when available, otherwise `"pvalue"`.
#'
#' @return A data.frame with added columns:
#'   `is_significant` and `activity_class`.
#' @export
classify_single_signatures <- function(
    single_table,
    nes_threshold = 1,
    p_value_threshold = 0.05,
    p_col = NULL
) {
  if (!is.data.frame(single_table)) {
    stop("`single_table` must be a data.frame.")
  }

  required_cols <- c("term", "NES")
  missing_cols <- setdiff(required_cols, colnames(single_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in `single_table`: ",
         paste(missing_cols, collapse = ", "))
  }

  if (is.null(p_col)) {
    if ("p.adjust" %in% colnames(single_table)) {
      p_col <- "p.adjust"
    } else if ("pvalue" %in% colnames(single_table)) {
      p_col <- "pvalue"
    } else {
      stop("No p-value column found. Please provide `p_col` explicitly.")
    }
  }

  if (!p_col %in% colnames(single_table)) {
    stop("`p_col` not found in `single_table`.")
  }

  if (!is.numeric(nes_threshold) || length(nes_threshold) != 1 || is.na(nes_threshold) || nes_threshold < 0) {
    stop("`nes_threshold` must be a single non-negative numeric value.")
  }

  if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1 || is.na(p_value_threshold) ||
      p_value_threshold < 0 || p_value_threshold > 1) {
    stop("`p_value_threshold` must be a single numeric value between 0 and 1.")
  }

  x <- single_table
  x$is_significant <- !is.na(x[[p_col]]) & x[[p_col]] < p_value_threshold
  x$activity_class <- "Not_significant"

  idx_pos <- x$is_significant & !is.na(x$NES) & x$NES >= nes_threshold
  idx_neg <- x$is_significant & !is.na(x$NES) & x$NES <= -nes_threshold

  x$activity_class[idx_pos] <- "Positive_enrichment"
  x$activity_class[idx_neg] <- "Negative_enrichment"

  x$activity_class <- factor(
    x$activity_class,
    levels = c("Positive_enrichment", "Negative_enrichment", "Not_significant")
  )

  x
}
