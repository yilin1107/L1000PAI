#' Summarize single-signature enrichment results
#'
#' @param gsea_table A data.frame returned by [run_l1000_gsea()].
#'
#' @return A standardized data.frame for single-signature libraries.
#' @export
summarize_single_signatures <- function(gsea_table) {
  if (!is.data.frame(gsea_table)) {
    stop("`gsea_table` must be a data.frame.")
  }

  required_cols <- c("term", "NES", "pvalue")
  missing_cols <- setdiff(required_cols, colnames(gsea_table))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in `gsea_table`: ",
         paste(missing_cols, collapse = ", "))
  }

  keep_cols <- intersect(
    c("term", "Description", "setSize", "enrichmentScore",
      "NES", "pvalue", "p.adjust", "qvalues", "rank",
      "leading_edge", "core_enrichment"),
    colnames(gsea_table)
  )

  out <- gsea_table[, keep_cols, drop = FALSE]
  out <- out[order(abs(out$NES), decreasing = TRUE), , drop = FALSE]
  rownames(out) <- NULL
  out
}
