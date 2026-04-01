#' Run the full L1000PAI workflow
#'
#' A one-step wrapper that performs rank generation (optional), L1000 library
#' loading, preranked GSEA, and either paired or single-signature summarization.
#'
#' @param x Either a named numeric rank vector or an object supported by
#'   [rank_genes()], such as a data.frame.
#' @param gene_col Column containing gene identifiers when `x` is not already
#'   a named numeric vector.
#' @param stat_col Column containing ranking statistics when `x` is not already
#'   a named numeric vector.
#' @param term2gene Optional two-column TERM2GENE data.frame. If supplied,
#'   `library` is ignored.
#' @param library Optional bundled library name to retrieve via
#'   [get_l1000_library()]. Used only when `term2gene` is `NULL`.
#' @param data_package Companion data package name. Default is `"L1000PAIdata"`.
#' @param library_mode One of `"auto"`, `"paired"`, or `"single"`.
#' @param score_direction One of `"auto"`, `"down_minus_up"`, or
#'   `"up_minus_down"`. For paired libraries, `"auto"` resolves according to
#'   the library type.
#' @param exponent Weighting exponent for GSEA.
#' @param min_gs_size Minimum gene set size.
#' @param max_gs_size Maximum gene set size.
#' @param eps Boundary for p-value calculation in GSEA.
#' @param pvalue_cutoff GSEA p-value cutoff.
#' @param p_adjust_method Multiple-testing adjustment method for GSEA.
#' @param seed Logical; whether to enable seeded execution in GSEA.
#' @param by GSEA backend. One of `"fgsea"` or `"DOSE"`.
#' @param nPerm Number of permutations for the `"DOSE"` backend. If `NULL` and
#'   `by = "DOSE"`, defaults to `1000`.
#' @param verbose Logical.
#' @param combine_p_method Method for combining paired p-values in
#'   [collapse_pairs()]. Default is `"max"`.
#' @param activity_threshold Minimum absolute activity score required for
#'   paired activation/inhibition classification.
#' @param classify_p_value_threshold Significance threshold used in classification.
#' @param classify_p_col Optional p-value column used in classification.
#' @param single_nes_threshold Minimum absolute NES required for single-signature
#'   positive/negative enrichment classification. Default is `1`.
#'
#' @return An `L1000PAIResult` object.
#' @export
run_l1000_pai <- function(
    x,
    gene_col = "gene",
    stat_col = "stat",
    term2gene = NULL,
    library = NULL,
    data_package = "L1000PAIdata",
    library_mode = c("auto", "paired", "single"),
    score_direction = c("auto", "down_minus_up", "up_minus_down"),
    exponent = 1,
    min_gs_size = 2,
    max_gs_size = 500,
    eps = 1e-10,
    pvalue_cutoff = 1,
    p_adjust_method = "BH",
    seed = FALSE,
    by = "fgsea",
    nPerm = NULL,
    verbose = FALSE,
    combine_p_method = "max",
    activity_threshold = 1,
    classify_p_value_threshold = 0.05,
    classify_p_col = NULL,
    single_nes_threshold = 1
) {
  library_mode <- match.arg(library_mode)
  score_direction <- match.arg(score_direction)

  # 1) 準備 ranked vector
  if (is.numeric(x) && !is.null(names(x))) {
    ranks <- x
    input_type <- "rank_vector"
  } else {
    ranks <- rank_genes(x, gene_col = gene_col, stat_col = stat_col)
    input_type <- class(x)[1]
  }

  # 2) 準備 library / term2gene
  if (is.null(term2gene)) {
    if (is.null(library)) {
      stop("Either `term2gene` or `library` must be provided.")
    }
    term2gene <- get_l1000_library(library = library, package = data_package)
    library_used <- library
  } else {
    if (!is.data.frame(term2gene) || ncol(term2gene) < 2) {
      stop("`term2gene` must be a data.frame with at least two columns.")
    }
    term2gene <- term2gene[, 1:2, drop = FALSE]
    colnames(term2gene) <- c("term", "gene")
    library_used <- if (is.null(library)) "user_supplied_term2gene" else library
  }

  # 2b) 偵測 library mode
  detected_mode <- detect_library_mode(term2gene)
  mode_used <- if (library_mode == "auto") detected_mode else library_mode

  # 2c) paired score direction
  score_direction_used <- resolve_score_direction(
    library = if (is.null(library)) NULL else library,
    score_direction = score_direction
  )

  # 3) GSEA
  gsea_table <- run_l1000_gsea(
    ranks = ranks,
    term2gene = term2gene,
    exponent = exponent,
    min_gs_size = min_gs_size,
    max_gs_size = max_gs_size,
    eps = eps,
    pvalue_cutoff = pvalue_cutoff,
    p_adjust_method = p_adjust_method,
    seed = seed,
    by = by,
    nPerm = nPerm,
    verbose = verbose
  )

  if (nrow(gsea_table) == 0) {
    stop(
      "No enriched L1000 signatures were returned by GSEA. ",
      "Try adjusting `min_gs_size`, `max_gs_size`, or the input ranking."
    )
  }

  # 4) paired vs single 分支
  if (mode_used == "paired") {
    result_tbl <- collapse_pairs(
      gsea_table = gsea_table,
      combine_p_method = combine_p_method,
      score_direction = score_direction_used
    )

    result_tbl <- classify_activity(
      pair_table = result_tbl,
      activity_threshold = activity_threshold,
      p_value_threshold = classify_p_value_threshold,
      p_col = classify_p_col,
      score_direction = score_direction_used
    )
  } else {
    result_tbl <- summarize_single_signatures(gsea_table)

    result_tbl <- classify_single_signatures(
      single_table = result_tbl,
      nes_threshold = single_nes_threshold,
      p_value_threshold = classify_p_value_threshold,
      p_col = classify_p_col
    )
  }

  metadata <- list(
    input_type = input_type,
    library = library_used,
    data_package = data_package,
    requested_library_mode = library_mode,
    detected_library_mode = detected_mode,
    library_mode_used = mode_used,
    n_ranked_genes = length(ranks),
    n_term2gene_rows = nrow(term2gene),
    gsea_params = list(
      exponent = exponent,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size,
      eps = eps,
      pvalue_cutoff = pvalue_cutoff,
      p_adjust_method = p_adjust_method,
      seed = seed,
      by = by,
      nPerm = nPerm,
      verbose = verbose
    ),
    pairing_params = list(
      combine_p_method = combine_p_method,
      score_direction = if (mode_used == "paired") score_direction_used else NA_character_
    ),
    classify_params = list(
      activity_threshold = activity_threshold,
      classify_p_value_threshold = classify_p_value_threshold,
      classify_p_col = classify_p_col,
      single_nes_threshold = single_nes_threshold
    )
  )

  L1000PAIResult(
    input_rank = ranks,
    gsea_table = gsea_table,
    pair_table = result_tbl,
    metadata = metadata
  )
}
