#' Rank genes for L1000PAI
#'
#' Convert differential expression results into a named numeric rank vector.
#'
#' @param x A data.frame containing gene identifiers and ranking statistics.
#' @param gene_col Column containing gene identifiers.
#' @param stat_col Column containing ranking statistics.
#'
#' @return A named numeric vector sorted in decreasing order.
#' @export
rank_genes <- function(x, gene_col = "gene", stat_col = "stat") {
  UseMethod("rank_genes")
}

#' @export
rank_genes.data.frame <- function(x, gene_col = "gene", stat_col = "stat") {
  if (!is.data.frame(x)) {
    stop("`x` must be a data.frame.")
  }

  if (!gene_col %in% colnames(x)) {
    stop("`gene_col` not found in `x`.")
  }

  if (!stat_col %in% colnames(x)) {
    stop("`stat_col` not found in `x`.")
  }

  out <- x[, c(gene_col, stat_col), drop = FALSE]
  out <- out[!is.na(out[[gene_col]]) & !is.na(out[[stat_col]]), , drop = FALSE]

  if (nrow(out) == 0) {
    stop("No valid rows remaining after removing NA values.")
  }

  out[[gene_col]] <- as.character(out[[gene_col]])
  out[[stat_col]] <- as.numeric(out[[stat_col]])

  # 若同一 gene 重複，保留絕對值最大的 statistic
  ord_abs <- order(abs(out[[stat_col]]), decreasing = TRUE)
  out <- out[ord_abs, , drop = FALSE]
  out <- out[!duplicated(out[[gene_col]]), , drop = FALSE]

  ord <- order(out[[stat_col]], decreasing = TRUE)
  out <- out[ord, , drop = FALSE]

  stats <- out[[stat_col]]
  names(stats) <- out[[gene_col]]

  stats
}
