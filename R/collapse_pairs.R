#' Collapse paired Up/Down perturbation signatures
#'
#' Pair L1000 GSEA results with `Up` and `Down` suffixes and summarize them
#' into perturbagen-level activity scores.
#'
#' @param gsea_table A data.frame containing GSEA results. It must contain
#'   at least the columns `term`, `NES`, and `pvalue`.
#' @param combine_p_method Method for combining paired p-values. Default is
#'   `"max"`, which conservatively uses the larger of the paired Up/Down
#'   p-values. `"fisher"` is also supported as an alternative exploratory
#'   option.
#' @param score_direction One of `"down_minus_up"` or `"up_minus_down"`.
#'   `"down_minus_up"` computes `NES_Down - NES_Up`, whereas `"up_minus_down"`
#'   computes `NES_Up - NES_Down`.
#'
#' @return A data.frame with paired `Up`/`Down` statistics, activity score,
#'   and combined p-values.
#' @export
collapse_pairs <- function(
    gsea_table,
    combine_p_method = "max",
    score_direction = c("down_minus_up", "up_minus_down")
) {
  score_direction <- match.arg(score_direction)

  if (!is.data.frame(gsea_table)) {
    stop("`gsea_table` must be a data.frame.")
  }

  required_cols <- c("term", "NES", "pvalue")
  missing_cols <- setdiff(required_cols, colnames(gsea_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in `gsea_table`: ",
         paste(missing_cols, collapse = ", "))
  }

  if (!combine_p_method %in% c("max", "fisher")) {
    stop("`combine_p_method` must be one of: 'max', 'fisher'.")
  }

  term <- as.character(gsea_table$term)

  direction <- ifelse(
    grepl("([[:space:]_.-]+)(Up|up)$", term),
    "Up",
    ifelse(grepl("([[:space:]_.-]+)(Down|down)$", term), "Down", NA_character_)
  )

  if (all(is.na(direction))) {
    stop("No paired Up/Down signatures were found in `gsea_table$term`.")
  }

  perturbagen <- term
  perturbagen <- sub("([[:space:]_.-]+)(Up|up)$", "", perturbagen)
  perturbagen <- sub("([[:space:]_.-]+)(Down|down)$", "", perturbagen)
  perturbagen <- trimws(perturbagen)

  x <- gsea_table
  x$term <- term
  x$direction <- direction
  x$perturbagen <- perturbagen
  x <- x[!is.na(x$direction), , drop = FALSE]

  up_tbl <- x[x$direction == "Up", , drop = FALSE]
  dn_tbl <- x[x$direction == "Down", , drop = FALSE]

  if (nrow(up_tbl) == 0 || nrow(dn_tbl) == 0) {
    stop("Both Up and Down signatures are required to form pairs.")
  }

  keep_cols <- intersect(
    c("perturbagen", "term", "Description", "setSize", "enrichmentScore",
      "NES", "pvalue", "p.adjust", "qvalues", "rank",
      "leading_edge", "core_enrichment"),
    colnames(x)
  )

  up_tbl <- up_tbl[, keep_cols, drop = FALSE]
  dn_tbl <- dn_tbl[, keep_cols, drop = FALSE]

  colnames(up_tbl) <- ifelse(
    colnames(up_tbl) == "perturbagen",
    "perturbagen",
    paste0(colnames(up_tbl), "_Up")
  )
  colnames(dn_tbl) <- ifelse(
    colnames(dn_tbl) == "perturbagen",
    "perturbagen",
    paste0(colnames(dn_tbl), "_Down")
  )

  pair_tbl <- merge(up_tbl, dn_tbl, by = "perturbagen", all = FALSE)

  if (nrow(pair_tbl) == 0) {
    stop("No complete Up/Down pairs were found after merging.")
  }

  p_up <- pmax(pair_tbl$pvalue_Up, .Machine$double.xmin)
  p_dn <- pmax(pair_tbl$pvalue_Down, .Machine$double.xmin)

  if (combine_p_method == "max") {
    pair_tbl$combined_stat <- NA_real_
    pair_tbl$p_combined <- pmax(p_up, p_dn)
  } else {
    pair_tbl$combined_stat <- -2 * (log(p_up) + log(p_dn))
    pair_tbl$p_combined <- stats::pchisq(
      pair_tbl$combined_stat,
      df = 4,
      lower.tail = FALSE
    )
  }

  if (score_direction == "down_minus_up") {
    pair_tbl$activity_score <- pair_tbl$NES_Down - pair_tbl$NES_Up
  } else {
    pair_tbl$activity_score <- pair_tbl$NES_Up - pair_tbl$NES_Down
  }

  pair_tbl$padj_combined <- stats::p.adjust(pair_tbl$p_combined, method = "BH")

  pair_tbl <- pair_tbl[order(pair_tbl$activity_score, decreasing = TRUE), , drop = FALSE]
  rownames(pair_tbl) <- NULL

  pair_tbl
}
