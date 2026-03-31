#' Run preranked GSEA on L1000 signatures
#'
#' @param ranks A named numeric vector of ranked genes.
#' @param term2gene A two-column data.frame for enrichment.
#' @param exponent Weighting exponent.
#' @param min_gs_size Minimum gene set size.
#' @param max_gs_size Maximum gene set size.
#' @param eps Boundary for p-value calculation.
#' @param pvalue_cutoff Adjusted p-value cutoff used by GSEA.
#' @param p_adjust_method Multiple-testing adjustment method.
#' @param seed Logical; whether to enable seeded execution in clusterProfiler GSEA.
#' @param by Backend used by clusterProfiler GSEA. One of `"fgsea"` or `"DOSE"`.
#' @param verbose Logical.
#'
#' @return A data.frame of GSEA results.
#' @export
run_l1000_gsea <- function(
    ranks,
    term2gene,
    exponent = 1,
    min_gs_size = 2,
    max_gs_size = 500,
    eps = 1e-10,
    pvalue_cutoff = 1,
    p_adjust_method = "BH",
    seed = FALSE,
    by = "fgsea",
    verbose = FALSE
) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Package `clusterProfiler` is required for `run_l1000_gsea()`. Please install it first.")
  }

  if (!is.numeric(ranks)) {
    stop("`ranks` must be a numeric vector.")
  }

  if (is.null(names(ranks)) || anyNA(names(ranks)) || any(names(ranks) == "")) {
    stop("`ranks` must be a named numeric vector with non-empty gene names.")
  }

  ranks <- ranks[!is.na(ranks)]
  if (length(ranks) == 0) {
    stop("`ranks` contains no valid numeric values after removing NA.")
  }

  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- sort(ranks, decreasing = TRUE)

  if (!is.data.frame(term2gene) || ncol(term2gene) < 2) {
    stop("`term2gene` must be a data.frame with at least two columns.")
  }

  term2gene <- term2gene[, 1:2, drop = FALSE]
  colnames(term2gene) <- c("term", "gene")
  term2gene$term <- as.character(term2gene$term)
  term2gene$gene <- as.character(term2gene$gene)
  term2gene <- term2gene[!is.na(term2gene$term) & !is.na(term2gene$gene), , drop = FALSE]
  term2gene <- term2gene[term2gene$term != "" & term2gene$gene != "", , drop = FALSE]
  term2gene <- unique(term2gene)

  if (nrow(term2gene) == 0) {
    stop("`term2gene` has no valid term-gene pairs.")
  }

  if (!by %in% c("fgsea", "DOSE")) {
    stop("`by` must be either 'fgsea' or 'DOSE'.")
  }

  if (!is.logical(seed) || length(seed) != 1 || is.na(seed)) {
    stop("`seed` must be TRUE or FALSE.")
  }

  if (!is.logical(seed) || length(seed) != 1 || is.na(seed)) {
    stop("`seed` must be TRUE or FALSE.")
  }

  # 某些乾淨 session 下 .Random.seed 尚未建立，
  # 但 clusterProfiler/DOSE 在 seed = TRUE 時會用到它
  if (isTRUE(seed) && !exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    set.seed(1L)
  }

  overlap_n <- sum(names(ranks) %in% term2gene$gene)

  if (overlap_n == 0) {
    stop(
      "No overlapping genes were found between `ranks` and `term2gene`. ",
      "Please check that gene identifiers use the same naming convention."
    )
  }

  gsea_res <- clusterProfiler::GSEA(
    geneList = ranks,
    exponent = exponent,
    minGSSize = min_gs_size,
    maxGSSize = max_gs_size,
    eps = eps,
    pvalueCutoff = pvalue_cutoff,
    pAdjustMethod = p_adjust_method,
    TERM2GENE = term2gene,
    verbose = verbose,
    seed = seed,
    by = by
  )

  res_df <- as.data.frame(gsea_res)

  if (nrow(res_df) == 0) {
    return(data.frame(
      term = character(),
      Description = character(),
      setSize = numeric(),
      enrichmentScore = numeric(),
      NES = numeric(),
      pvalue = numeric(),
      p.adjust = numeric(),
      qvalues = numeric(),
      rank = numeric(),
      leading_edge = character(),
      core_enrichment = character(),
      stringsAsFactors = FALSE
    ))
  }

  if ("ID" %in% colnames(res_df)) {
    colnames(res_df)[colnames(res_df) == "ID"] <- "term"
  }

  res_df
}
