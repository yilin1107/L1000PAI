#' Inspect an ExpressionSet object
#'
#' Print a concise overview of an [Biobase::ExpressionSet], including assay
#' dimensions, pData columns, sample names, and low-cardinality pData values.
#'
#' @importFrom magrittr %>%
#'
#' @param eset An [Biobase::ExpressionSet] object.
#' @param n Number of rows to print from pData and fData. Default is `10`.
#' @param max_unique Maximum number of unique values to display for a pData
#'   column. Default is `10`.
#'
#' @return Invisibly returns a list with `expr_dim`, `pdata`, and `fdata`.
#' @export
inspect_eset <- function(eset, n = 10, max_unique = 10) {
  if (!methods::is(eset, "ExpressionSet")) {
    stop("`eset` must be an ExpressionSet object.")
  }

  cat("=== ExpressionSet basic info ===\n")
  print(eset)

  cat("\n=== exprs dimension ===\n")
  print(dim(Biobase::exprs(eset)))

  pd <- Biobase::pData(eset)
  fd <- Biobase::fData(eset)

  cat("\n=== pData head ===\n")
  print(utils::head(pd, n))

  cat("\n=== pData column names ===\n")
  print(colnames(pd))

  cat("\n=== sample names ===\n")
  print(Biobase::sampleNames(eset))

  cat("\n=== pData unique values (limited columns only) ===\n")
  for (cc in colnames(pd)) {
    vals <- unique(as.character(pd[[cc]]))
    vals <- vals[!is.na(vals)]
    if (length(vals) <= max_unique) {
      cat("\n[", cc, "]\n", sep = "")
      print(vals)
    }
  }

  if (ncol(fd) > 0) {
    cat("\n=== fData head ===\n")
    print(utils::head(fd, n))
    cat("\n=== fData column names ===\n")
    print(colnames(fd))
  }

  invisible(list(
    expr_dim = dim(Biobase::exprs(eset)),
    pdata = pd,
    fdata = fd
  ))
}


#' Check whether an expression matrix looks like raw counts
#'
#' @param expr_mat A numeric matrix.
#' @param tol Tolerance for integer-likeness. Default is `1e-8`.
#'
#' @return Logical scalar.
#' @export
looks_like_raw_counts <- function(expr_mat, tol = 1e-8) {
  if (!is.matrix(expr_mat)) {
    stop("`expr_mat` must be a matrix.")
  }

  is_integer_like <- all(abs(expr_mat - round(expr_mat)) < tol, na.rm = TRUE)
  all(expr_mat >= 0, na.rm = TRUE) && is_integer_like
}


#' Choose a feature identifier column from fData
#'
#' @param eset An [Biobase::ExpressionSet] object.
#' @param feature_id_priority Character vector of preferred columns.
#'
#' @return A single column name, or `NA_character_` if none found.
#' @export
choose_feature_id_col <- function(
    eset,
    feature_id_priority = c("SYMBOL", "ENTREZID", "ENSEMBL")
) {
  if (!methods::is(eset, "ExpressionSet")) {
    stop("`eset` must be an ExpressionSet object.")
  }

  fd_cols <- colnames(Biobase::fData(eset))
  hit <- feature_id_priority[feature_id_priority %in% fd_cols]

  if (length(hit) == 0) {
    return(NA_character_)
  }

  hit[1]
}


#' Run limma differential expression from an ExpressionSet
#'
#' @param eset An [Biobase::ExpressionSet] object.
#' @param group_col Column name in `pData(eset)` containing group labels.
#' @param control_label Control-group label.
#' @param case_label Case/treatment-group label.
#' @param use_voom One of `"auto"`, `TRUE`, or `FALSE`. If `"auto"`, voom is
#'   used when the expression matrix looks like raw counts.
#'
#' @return A list containing `eset_subset`, `design`, `contrast_mat`,
#'   `fit2`, `tt_all`, and `group`.
#' @export
run_limma_from_eset <- function(
    eset,
    group_col,
    control_label,
    case_label,
    use_voom = "auto"
) {
  if (!methods::is(eset, "ExpressionSet")) {
    stop("`eset` must be an ExpressionSet object.")
  }

  pd <- Biobase::pData(eset)
  if (!group_col %in% colnames(pd)) {
    stop("`group_col` not found in pData(eset).")
  }

  group_raw <- as.character(pd[[group_col]])

  if (!all(c(control_label, case_label) %in% unique(group_raw))) {
    stop("`control_label` and/or `case_label` not found in pData(eset)[[group_col]].")
  }

  keep_samples <- group_raw %in% c(control_label, case_label)
  eset2 <- eset[, keep_samples]
  expr_mat <- Biobase::exprs(eset2)
  pd2 <- Biobase::pData(eset2)

  group <- factor(
    as.character(pd2[[group_col]]),
    levels = c(control_label, case_label)
  )

  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- c("CTRL", "CASE")

  contrast_mat <- limma::makeContrasts(contrasts = "CASE-CTRL", levels = design)

  auto_counts <- looks_like_raw_counts(expr_mat)

  if (identical(use_voom, "auto")) {
    use_voom_final <- auto_counts
  } else {
    use_voom_final <- isTRUE(use_voom)
  }

  if (use_voom_final) {
    v <- limma::voom(expr_mat, design = design, plot = TRUE)
    fit <- limma::lmFit(v, design)
  } else {
    fit <- limma::lmFit(expr_mat, design)
  }

  fit2 <- limma::contrasts.fit(fit, contrast_mat)
  fit2 <- limma::eBayes(fit2)

  tt_all <- limma::topTable(fit2, number = Inf, sort.by = "none")
  tt_all$feature_id <- rownames(tt_all)

  list(
    eset_subset = eset2,
    design = design,
    contrast_mat = contrast_mat,
    fit2 = fit2,
    tt_all = tt_all,
    group = group,
    use_voom = use_voom_final
  )
}


#' Add gene symbols to a limma result table derived from an ExpressionSet
#'
#' @param tt_all A data.frame returned by [run_limma_from_eset()]$tt_all.
#' @param eset An [Biobase::ExpressionSet] object used to generate `tt_all`.
#' @param feature_col Optional feature column in `fData(eset)`. If `NULL`,
#'   [choose_feature_id_col()] is used.
#' @param feature_id_priority Priority order used when `feature_col=NULL`.
#'
#' @return A data.frame with a `SYMBOL` column.
#' @export
add_symbol_to_limma_table <- function(
    tt_all,
    eset,
    feature_col = NULL,
    feature_id_priority = c("SYMBOL", "ENTREZID", "ENSEMBL")
) {
  if (!methods::is(eset, "ExpressionSet")) {
    stop("`eset` must be an ExpressionSet object.")
  }
  if (!is.data.frame(tt_all)) {
    stop("`tt_all` must be a data.frame.")
  }
  if (!"feature_id" %in% colnames(tt_all)) {
    stop("`tt_all` must contain a `feature_id` column.")
  }

  fd <- Biobase::fData(eset)
  fd$feature_id <- rownames(fd)

  if (is.null(feature_col)) {
    feature_col <- choose_feature_id_col(eset, feature_id_priority)
  }

  out <- tt_all

  if (is.na(feature_col)) {
    warning("No supported feature column found; using `feature_id` as SYMBOL.")
    out$SYMBOL <- out$feature_id
    return(out)
  }

  if (!feature_col %in% colnames(fd)) {
    stop("`feature_col` not found in fData(eset).")
  }

  key_vec <- as.character(fd[out$feature_id, feature_col])

  if (feature_col == "SYMBOL") {
    out$SYMBOL <- key_vec
  } else if (feature_col == "ENTREZID") {
    out$SYMBOL <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = key_vec,
      keytype = "ENTREZID",
      column = "SYMBOL",
      multiVals = "first"
    )
  } else if (feature_col == "ENSEMBL") {
    out$SYMBOL <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = key_vec,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
  } else {
    stop("Unsupported `feature_col`: ", feature_col)
  }

  out
}


#' Collapse limma result rows to unique gene symbols
#'
#' @param tt_all A limma result table with a `SYMBOL` column.
#' @param stat_col Column used to prioritize duplicates. Default is `"t"`.
#'
#' @return A gene-collapsed data.frame.
#' @export
collapse_limma_to_genes <- function(tt_all, stat_col = "t") {
  if (!is.data.frame(tt_all)) {
    stop("`tt_all` must be a data.frame.")
  }
  if (!"SYMBOL" %in% colnames(tt_all)) {
    stop("`tt_all` must contain a `SYMBOL` column.")
  }
  if (!stat_col %in% colnames(tt_all)) {
    stop("`stat_col` not found in `tt_all`.")
  }

  tt_all %>%
    dplyr::filter(!is.na(.data$SYMBOL), .data$SYMBOL != "") %>%
    dplyr::arrange(dplyr::desc(abs(.data[[stat_col]]))) %>%
    dplyr::distinct(.data$SYMBOL, .keep_all = TRUE)
}


#' Convert an ExpressionSet to rank_input for L1000PAI
#'
#' @param eset An [Biobase::ExpressionSet] object.
#' @param group_col Column name in `pData(eset)` containing group labels.
#' @param control_label Control-group label.
#' @param case_label Case/treatment-group label.
#' @param feature_col Optional feature column in `fData(eset)`.
#' @param feature_id_priority Priority order used when `feature_col=NULL`.
#' @param stat_col Statistic to use for ranking. Default is `"t"`.
#' @param use_voom One of `"auto"`, `TRUE`, or `FALSE`.
#' @param jitter_ties Logical; whether to add very small offsets to break ties.
#'
#' @return A list containing `rank_input`, `tt_gene`, `tt_all`, `fit2`,
#'   `design`, `contrast_mat`, `gene_universe`, and `use_voom`.
#' @export
eset_to_rank_input <- function(
    eset,
    group_col,
    control_label,
    case_label,
    feature_col = NULL,
    feature_id_priority = c("SYMBOL", "ENTREZID", "ENSEMBL"),
    stat_col = "t",
    use_voom = "auto",
    jitter_ties = TRUE
) {
  limma_res <- run_limma_from_eset(
    eset = eset,
    group_col = group_col,
    control_label = control_label,
    case_label = case_label,
    use_voom = use_voom
  )

  tt_all2 <- add_symbol_to_limma_table(
    tt_all = limma_res$tt_all,
    eset = limma_res$eset_subset,
    feature_col = feature_col,
    feature_id_priority = feature_id_priority
  )

  tt_gene <- collapse_limma_to_genes(tt_all2, stat_col = stat_col)

  gene_universe <- unique(tt_gene$SYMBOL)

  rank_input <- tt_gene %>%
    dplyr::transmute(
      gene = .data$SYMBOL,
      stat = .data[[stat_col]]
    ) %>%
    dplyr::arrange(dplyr::desc(.data$stat))

  if (jitter_ties) {
    rank_input <- rank_input %>%
      dplyr::mutate(stat = .data$stat + seq_along(.data$stat) * 1e-12)
  }

  list(
    rank_input = rank_input,
    tt_gene = tt_gene,
    tt_all = tt_all2,
    fit2 = limma_res$fit2,
    design = limma_res$design,
    contrast_mat = limma_res$contrast_mat,
    gene_universe = gene_universe,
    use_voom = limma_res$use_voom
  )
}
