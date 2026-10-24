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
#' @param rng_seed Optional explicit integer RNG seed. When supplied, it is set
#'   immediately before GSEA and also enables the backend's seeded execution.
#'   This is preferred over relying on the session's pre-existing RNG state.
#' @param keep_all_terms Logical; for the `"fgsea"` backend, restore every
#'   size-eligible pathway that the backend evaluated, including pathways for
#'   which adaptive multilevel p-value estimation returned `NA`. Such rows are
#'   retained with conservative `pvalue = 1`, `p.adjust = 1`, and `qvalue = 1`.
#'   Default is `FALSE` for backward compatibility. Set `TRUE` when a complete
#'   prespecified testing universe is required.
#' @param complete_prespecified_terms Logical; for the `"fgsea"` backend,
#'   conservatively retain every prespecified TERM2GENE term with at least one
#'   ranked-gene overlap, including terms outside `min_gs_size` or
#'   `max_gs_size`. Size-ineligible rows receive `pvalue = 1`, `NES = 0`, and
#'   are not treated as estimated enrichment results. This does not relax the
#'   size thresholds. Default is `FALSE`.
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
    nPerm = NULL,
    rng_seed = NULL,
    keep_all_terms = FALSE,
    complete_prespecified_terms = FALSE
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

  if (!is.null(rng_seed)) {
    if (!is.numeric(rng_seed) || length(rng_seed) != 1L ||
        !is.finite(rng_seed) || rng_seed < 0 ||
        rng_seed > .Machine$integer.max || rng_seed != floor(rng_seed)) {
      stop("`rng_seed` must be NULL or a single integer between 0 and .Machine$integer.max.")
    }
    rng_seed <- as.integer(rng_seed)
    set.seed(rng_seed)
    backend_seed <- TRUE
  } else if (isTRUE(seed)) {
    set.seed(1L)
    backend_seed <- TRUE
  } else {
    backend_seed <- FALSE
  }

  if (!is.logical(keep_all_terms) || length(keep_all_terms) != 1L ||
      is.na(keep_all_terms)) {
    stop("`keep_all_terms` must be a single non-missing logical value.")
  }
  if (!is.logical(complete_prespecified_terms) ||
      length(complete_prespecified_terms) != 1L ||
      is.na(complete_prespecified_terms)) {
    stop("`complete_prespecified_terms` must be a single non-missing logical value.")
  }
  if ((isTRUE(keep_all_terms) || isTRUE(complete_prespecified_terms)) &&
      (!is.numeric(pvalue_cutoff) || length(pvalue_cutoff) != 1L ||
       !is.finite(pvalue_cutoff) || pvalue_cutoff != 1)) {
    stop("Complete-universe modes require `pvalue_cutoff = 1`.")
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
      seed = backend_seed,
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
      seed = backend_seed,
      by = by
    )
  }

  out <- as.data.frame(gsea_obj)

  if ((isTRUE(keep_all_terms) || isTRUE(complete_prespecified_terms)) &&
      by == "fgsea") {
    observed_terms <- if ("ID" %in% colnames(out)) {
      as.character(out$ID)
    } else if ("term" %in% colnames(out)) {
      as.character(out$term)
    } else {
      character()
    }
    eligible <- .eligible_l1000_terms(
      gene_ids = names(geneList),
      term2gene = term2gene,
      min_gs_size = min_gs_size,
      max_gs_size = max_gs_size
    )
    completion_universe <- if (isTRUE(complete_prespecified_terms)) {
      .overlapping_l1000_terms(names(geneList), term2gene)
    } else {
      eligible
    }
    missing_terms <- setdiff(completion_universe$term, observed_terms)
    if (length(missing_terms) > 0L) {
      fallback_rows <- .score_missing_fgsea_terms(
        gene_list = geneList,
        term2gene = term2gene,
        eligible = completion_universe,
        missing_terms = missing_terms,
        exponent = exponent
      )
      out <- rbind(out, fallback_rows[, colnames(out), drop = FALSE])
      # Recompute the within-arm adjustment over the complete prespecified
      # universe. Formal paired inference is recomputed downstream from these
      # raw p-values across both arms.
      out$p.adjust <- stats::p.adjust(out$pvalue, method = p_adjust_method)
      ineligible_n <- sum(!missing_terms %in% eligible$term)
      warning(
        "Retained ", length(missing_terms),
        " missing fgsea term(s) as conservative pvalue=1, NES=0 rows; ",
        ineligible_n, " were outside the fixed gene-set size limits.",
        call. = FALSE
      )
    }
  }

  # Standardize the returned enrichment term column across backends.
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

  # Move the standardized term column to the front for downstream code.
  out <- out[, c("term", setdiff(colnames(out), "term")), drop = FALSE]

  out
}

.eligible_l1000_terms <- function(gene_ids, term2gene, min_gs_size, max_gs_size) {
  sizes <- .overlapping_l1000_terms(gene_ids, term2gene)
  sizes <- sizes[
    sizes$setSize >= min_gs_size & sizes$setSize <= max_gs_size,
    ,
    drop = FALSE
  ]
  sizes$term <- as.character(sizes$term)
  sizes
}

.overlapping_l1000_terms <- function(gene_ids, term2gene) {
  in_rank <- term2gene$gene %in% gene_ids
  if (!any(in_rank)) {
    return(data.frame(term = character(), setSize = integer()))
  }
  sizes <- stats::aggregate(
    term2gene$gene[in_rank],
    by = list(term = term2gene$term[in_rank]),
    FUN = function(x) length(unique(x))
  )
  colnames(sizes)[2L] <- "setSize"
  sizes$term <- as.character(sizes$term)
  sizes
}

.score_missing_fgsea_terms <- function(
    gene_list,
    term2gene,
    eligible,
    missing_terms,
    exponent
) {
  # fgseaMultilevel can return NA for extremely unbalanced pathways. Preserve
  # their observed enrichment geometry while assigning conservative p=1 so a
  # prespecified paired library never shrinks in an outcome-dependent manner.
  terms <- as.character(missing_terms)
  set_sizes <- stats::setNames(eligible$setSize, eligible$term)
  result <- lapply(terms, function(term) {
    genes <- unique(term2gene$gene[term2gene$term == term])
    genes <- intersect(genes, names(gene_list))
    score_info <- fgsea::calcGseaStat(
      stats = gene_list,
      selectedStats = match(genes, names(gene_list)),
      gseaParam = exponent,
      returnAllExtremes = TRUE,
      returnLeadingEdge = TRUE
    )
    leading_edge <- names(gene_list)[score_info$leadingEdge]
    data.frame(
      ID = term,
      Description = term,
      setSize = as.integer(set_sizes[[term]]),
      enrichmentScore = as.numeric(score_info$res),
      NES = 0,
      pvalue = 1,
      p.adjust = 1,
      qvalue = 1,
      # The backend did not estimate this pathway; avoid fabricating a
      # clusterProfiler peak-rank diagnostic from the raw hit-index curve.
      rank = NA_integer_,
      leading_edge = "tags=NA%, list=NA%, signal=NA%",
      core_enrichment = paste(leading_edge, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, result)
}
