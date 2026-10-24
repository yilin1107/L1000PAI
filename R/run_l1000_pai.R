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
#' @param library Optional L1000 library name to retrieve via
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
#' @param collect_annotation_tables Logical; if `TRUE`, package-level
#'   annotation tables are collected automatically for supported
#'   libraries (`chem`, `ligand`, `kinase_gpcr`). Default is `TRUE`.
#' @param rng_seed Optional explicit integer RNG seed passed to
#'   [run_l1000_gsea()].
#' @param keep_all_terms Logical; retain the complete size-eligible term
#'   universe when the fgsea backend cannot estimate a p-value for a pathway.
#'   Such terms receive conservative `pvalue = 1` and `NES = 0`. Default is
#'   `FALSE` for backward compatibility.
#' @param complete_prespecified_terms Logical; conservatively retain every
#'   prespecified term with ranked-gene overlap, even when its overlap is
#'   outside the fixed GSEA size limits. This does not relax those limits.
#' @param pvalue_semantics Interpretation of paired GSEA p-values passed to
#'   [collapse_pairs()]. Default is `"observed_tail"`.
#' @param pair_p_adjust_method Multiple-testing method for paired directional
#'   hypotheses. Default is `"BH"`.
#' @param joint_null Optional synchronized Up/Down null-draw data.frame passed
#'   to [score_joint_activity()]. When `NULL`, no GLS evidence is produced.
#' @param gls_ridge Non-negative relative covariance ridge used when
#'   `joint_null` is supplied. Default is zero.
#' @param gls_min_null_draws Minimum synchronized null draws per perturbagen.
#' @param classify_score_col Score column used by [classify_activity()]. The
#'   default preserves legacy NES-difference thresholds; use
#'   `"activity_gls_z"` to threshold standardized GLS activity.
#' @param classify_require_coherent Logical; require the adjusted directional
#'   conjunction whenever it is available, even if `classify_p_col` selects a
#'   global or legacy p-value. Default is `TRUE`.
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
    single_nes_threshold = 1,
    collect_annotation_tables = TRUE,
    rng_seed = NULL,
    keep_all_terms = FALSE,
    complete_prespecified_terms = FALSE,
    pvalue_semantics = c("observed_tail", "two_sided"),
    pair_p_adjust_method = "BH",
    joint_null = NULL,
    gls_ridge = 0,
    gls_min_null_draws = 50L,
    classify_score_col = "activity_score",
    classify_require_coherent = TRUE
) {
  library_mode <- match.arg(library_mode)
  score_direction <- match.arg(score_direction)
  pvalue_semantics <- match.arg(pvalue_semantics)

  # Standardize the input as a named ranked vector.
  if (is.numeric(x) && !is.null(names(x))) {
    ranks <- x
    input_type <- "rank_vector"
  } else {
    ranks <- rank_genes(x, gene_col = gene_col, stat_col = stat_col)
    input_type <- class(x)[1]
  }

  # Resolve the library source or validate the user-supplied TERM2GENE table.
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

  # Detect or honor the requested library mode.
  detected_mode <- detect_library_mode(term2gene)
  mode_used <- if (library_mode == "auto") detected_mode else library_mode
  if (mode_used != "paired" && !is.null(joint_null)) {
    stop("`joint_null` is only supported for paired signature libraries.")
  }

  # Resolve paired-library score direction.
  score_direction_used <- resolve_score_direction(
    library = if (is.null(library)) NULL else library,
    score_direction = score_direction
  )

  # Run GSEA against the selected perturbation library.
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
    verbose = verbose,
    rng_seed = rng_seed,
    keep_all_terms = keep_all_terms,
    complete_prespecified_terms = complete_prespecified_terms
  )

  if (nrow(gsea_table) == 0) {
    stop(
      "No enriched L1000 signatures were returned by GSEA. ",
      "Try adjusting `min_gs_size`, `max_gs_size`, or the input ranking."
    )
  }

  # Summarize paired or single-signature results.
  if (mode_used == "paired") {
    result_tbl <- collapse_pairs(
      gsea_table = gsea_table,
      combine_p_method = combine_p_method,
      score_direction = score_direction_used,
      pvalue_semantics = pvalue_semantics,
      pair_p_adjust_method = pair_p_adjust_method
    )

    if (!is.null(joint_null)) {
      result_tbl <- score_joint_activity(
        pair_table = result_tbl,
        joint_null = joint_null,
        score_direction = score_direction_used,
        ridge = gls_ridge,
        min_null_draws = gls_min_null_draws,
        p_adjust_method = pair_p_adjust_method
      )
    }

    result_tbl <- classify_activity(
      pair_table = result_tbl,
      activity_threshold = activity_threshold,
      p_value_threshold = classify_p_value_threshold,
      p_col = classify_p_col,
      score_col = classify_score_col,
      require_coherent = classify_require_coherent,
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
      rng_seed = rng_seed,
      keep_all_terms = keep_all_terms,
      complete_prespecified_terms = complete_prespecified_terms,
      verbose = verbose
    ),
    pairing_params = list(
      combine_p_method = combine_p_method,
      score_direction = if (mode_used == "paired") score_direction_used else NA_character_,
      pvalue_semantics = if (mode_used == "paired") pvalue_semantics else NA_character_,
      pair_p_adjust_method = pair_p_adjust_method
    ),
    classify_params = list(
      activity_threshold = activity_threshold,
      classify_p_value_threshold = classify_p_value_threshold,
      classify_p_col = classify_p_col,
      classify_score_col = classify_score_col,
      classify_require_coherent = classify_require_coherent,
      single_nes_threshold = single_nes_threshold
    ),
    joint_model_params = list(
      enabled = !is.null(joint_null),
      method = if (is.null(joint_null)) NA_character_ else "pair_specific_empirical_gls",
      covariance_source = if (is.null(joint_null)) "unavailable" else "synchronized_joint_null",
      ridge = gls_ridge,
      min_null_draws = gls_min_null_draws
    ),
    annotation_tables = NULL
  )

  if (isTRUE(collect_annotation_tables) &&
      library_used %in% available_l1000_annotation_catalogs()) {
    metadata$annotation_tables <- tryCatch(
      build_l1000_annotation_tables(result_tbl, library = library_used),
      error = function(e) {
        warning(
          "Failed to collect packaged annotation tables for `",
          library_used, "`: ", conditionMessage(e),
          call. = FALSE
        )
        NULL
      }
    )
  }

  L1000PAIResult(
    input_rank = ranks,
    gsea_table = gsea_table,
    pair_table = result_tbl,
    metadata = metadata
  )
}
