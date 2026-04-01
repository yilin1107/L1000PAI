#' Run L1000 GSEA
#'
#' @param ranks A named numeric rank vector, or a two-column data.frame with
#'   columns `gene` and `stat`.
#' @param term2gene A two-column TERM2GENE data.frame with columns `term` and `gene`.
#' @param exponent Weighting exponent for GSEA.
#' @param min_gs_size Minimum gene set size.
#' @param max_gs_size Maximum gene set size.
#' @param eps Boundary for p-value calculation.
#' @param pvalue_cutoff P-value cutoff.
#' @param p_adjust_method Multiple-testing correction method.
#' @param seed Logical.
#' @param by GSEA backend. One of `"fgsea"` or `"DOSE"`.
#' @param verbose Logical.
#' @param nPerm Number of permutations for the `"DOSE"` backend. If `NULL` and
#'   `by = "DOSE"`, defaults to `1000`.
#'
#' @return A data.frame of GSEA results with a standardized `term` column.
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
    verbose = FALSE,
    nPerm = NULL
) {
  if (is.data.frame(ranks)) {
    if (!all(c("gene", "stat") %in% colnames(ranks))) {
      stop("If `ranks` is a data.frame, it must contain columns `gene` and `stat`.")
    }
    geneList <- ranks$stat
    names(geneList) <- ranks$gene
    geneList <- sort(geneList, decreasing = TRUE)
  } else if (is.numeric(ranks) && !is.null(names(ranks))) {
    geneList <- sort(ranks, decreasing = TRUE)
  } else {
    stop("`ranks` must be a named numeric vector or a data.frame with columns `gene` and `stat`.")
  }

  if (!is.data.frame(term2gene) || ncol(term2gene) < 2) {
    stop("`term2gene` must be a data.frame with at least two columns.")
  }

  term2gene <- term2gene[, 1:2, drop = FALSE]
  colnames(term2gene) <- c("term", "gene")

  overlap_n <- sum(names(geneList) %in% term2gene$gene)
  if (overlap_n == 0) {
    stop(
      "No overlapping genes were found between `ranks` and `term2gene`. ",
      "Please check that gene identifiers use the same naming convention."
    )
  }

  gsea_obj <- if (by == "DOSE") {
    if (is.null(nPerm)) {
      nPerm <- 1000
    }

    clusterProfiler::GSEA(
      geneList = geneList,
      exponent = exponent,
      minGSSize = min_gs_size,
      maxGSSize = max_gs_size,
      eps = eps,
      pvalueCutoff = pvalue_cutoff,
      pAdjustMethod = p_adjust_method,
      TERM2GENE = term2gene,
      verbose = verbose,
      seed = seed,
      by = by,
      nPerm = nPerm
    )
  } else {
    clusterProfiler::GSEA(
      geneList = geneList,
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
  }

  out <- as.data.frame(gsea_obj)

  # ---- 關鍵修正：統一欄名 ----
  if (!"term" %in% colnames(out)) {
    if ("ID" %in% colnames(out)) {
      out$term <- as.character(out$ID)
    } else {
      stop(
        "The GSEA result does not contain a `term` column or an `ID` column. ",
        "Available columns are: ", paste(colnames(out), collapse = ", ")
      )
    }
  }

  # term 放到前面，方便 downstream
  out <- out[, c("term", setdiff(colnames(out), "term")), drop = FALSE]

  out
}
