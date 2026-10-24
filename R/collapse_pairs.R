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
#' @param pvalue_semantics Interpretation of the input GSEA p-values.
#'   `"observed_tail"` (the default for clusterProfiler/fgsea and DOSE GSEA)
#'   treats each p-value as conditional on the observed enrichment direction
#'   and assigns p = 1 to the opposite direction. `"two_sided"` converts a
#'   conventional two-sided p-value to directional one-sided p-values using
#'   the NES sign.
#' @param pair_p_adjust_method Multiple-testing method used for paired and
#'   directional conjunction p-values. Default is `"BH"`.
#'
#' @return A data.frame with paired `Up`/`Down` statistics, activity score,
#'   legacy combined p-values, directional activation/inhibition conjunction
#'   p-values, and a unified `term` column. Activation and inhibition
#'   hypotheses are adjusted together. Pairs with a missing NES are retained
#'   with missing activity and directional state; they cannot support a
#'   directional activity call. Other estimable pairs are still processed.
#' @export
collapse_pairs <- function(
    gsea_table,
    combine_p_method = "max",
    score_direction = c("down_minus_up", "up_minus_down"),
    pvalue_semantics = c("observed_tail", "two_sided"),
    pair_p_adjust_method = "BH"
) {
  score_direction <- match.arg(score_direction)
  pvalue_semantics <- match.arg(pvalue_semantics)

  if (!is.data.frame(gsea_table)) {
    stop("`gsea_table` must be a data.frame.")
  }

  required_cols <- c("term", "NES", "pvalue")
  missing_cols <- setdiff(required_cols, colnames(gsea_table))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns in `gsea_table`: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  if (!combine_p_method %in% c("max", "fisher")) {
    stop("`combine_p_method` must be one of: 'max', 'fisher'.")
  }
  if (!is.character(pair_p_adjust_method) || length(pair_p_adjust_method) != 1L ||
      !pair_p_adjust_method %in% stats::p.adjust.methods) {
    stop(
      "`pair_p_adjust_method` must be one of: ",
      paste(stats::p.adjust.methods, collapse = ", "), "."
    )
  }

  if (!is.numeric(gsea_table$NES)) {
    stop("`gsea_table$NES` must be numeric.")
  }
  invalid_nes <- !is.na(gsea_table$NES) & !is.finite(gsea_table$NES)
  if (any(invalid_nes)) {
    stop("`gsea_table$NES` must contain finite values, or NA.")
  }
  if (!is.numeric(gsea_table$pvalue)) {
    stop("`gsea_table$pvalue` must be numeric.")
  }
  invalid_p <- !is.na(gsea_table$pvalue) &
    (!is.finite(gsea_table$pvalue) | gsea_table$pvalue < 0 | gsea_table$pvalue > 1)
  if (any(invalid_p)) {
    stop("`gsea_table$pvalue` must contain values between 0 and 1, or NA.")
  }

  term <- as.character(gsea_table$term)

  direction <- ifelse(
    grepl("([[:space:]_.-]+)(Up|up)$", term),
    "Up",
    ifelse(grepl("([[:space:]_.-]+)(Down|down|Dn|dn)$", term), "Down", NA_character_)
  )

  if (all(is.na(direction))) {
    stop("No paired Up/Down signatures were found in `gsea_table$term`.")
  }

  perturbagen <- term
  perturbagen <- sub("([[:space:]_.-]+)(Up|up)$", "", perturbagen)
  perturbagen <- sub("([[:space:]_.-]+)(Down|down|Dn|dn)$", "", perturbagen)
  perturbagen <- trimws(perturbagen)

  x <- gsea_table
  x$term <- term
  x$direction <- direction
  x$perturbagen <- perturbagen
  x <- x[!is.na(x$direction), , drop = FALSE]

  duplicate_signature <- duplicated(x[, c("perturbagen", "direction"), drop = FALSE])
  if (any(duplicate_signature)) {
    duplicate_ids <- unique(paste0(
      x$perturbagen[duplicate_signature], " (", x$direction[duplicate_signature], ")"
    ))
    stop(
      "Duplicate directional signatures would create ambiguous pairs: ",
      paste(utils::head(duplicate_ids, 10L), collapse = ", "),
      if (length(duplicate_ids) > 10L) ", ..." else "."
    )
  }

  up_tbl <- x[x$direction == "Up", , drop = FALSE]
  dn_tbl <- x[x$direction == "Down", , drop = FALSE]

  if (nrow(up_tbl) == 0 || nrow(dn_tbl) == 0) {
    stop("Both Up and Down signatures are required to form pairs.")
  }

  incomplete_ids <- union(
    setdiff(up_tbl$perturbagen, dn_tbl$perturbagen),
    setdiff(dn_tbl$perturbagen, up_tbl$perturbagen)
  )
  if (length(incomplete_ids) > 0L) {
    warning(
      "Incomplete Up/Down signatures were excluded for perturbagen(s): ",
      paste(utils::head(incomplete_ids, 10L), collapse = ", "),
      if (length(incomplete_ids) > 10L) ", ..." else ".",
      call. = FALSE
    )
  }

  keep_cols <- intersect(
    c(
      "perturbagen", "term", "Description", "setSize", "enrichmentScore",
      "NES", "pvalue", "p.adjust", "qvalues", "rank",
      "leading_edge", "core_enrichment"
    ),
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

  p_up_raw <- pair_tbl$pvalue_Up
  p_dn_raw <- pair_tbl$pvalue_Down
  p_up <- pmax(p_up_raw, .Machine$double.xmin)
  p_dn <- pmax(p_dn_raw, .Machine$double.xmin)

  pair_tbl$p_both_arms <- pmax(p_up_raw, p_dn_raw)
  pair_tbl$padj_both_arms <- .adjust_complete_p(
    pair_tbl$p_both_arms,
    method = pair_p_adjust_method
  )

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

  pair_tbl$activity_score_legacy <- pair_tbl$activity_score

  pair_tbl$padj_combined <- .adjust_complete_p(
    pair_tbl$p_combined,
    method = pair_p_adjust_method
  )

  pair_tbl <- .add_directional_pair_inference(
    pair_tbl = pair_tbl,
    score_direction = score_direction,
    pvalue_semantics = pvalue_semantics,
    p_adjust_method = pair_p_adjust_method
  )

  # Add a unified term column for downstream workflows shared with single-mode results.
  pair_tbl$term <- as.character(pair_tbl$perturbagen)

  # Move the unified term column to the front for easier downstream handling.
  pair_tbl <- pair_tbl[, c("term", setdiff(colnames(pair_tbl), "term")), drop = FALSE]

  pair_tbl <- pair_tbl[order(pair_tbl$activity_score, decreasing = TRUE), , drop = FALSE]
  rownames(pair_tbl) <- NULL

  pair_tbl
}

.add_directional_pair_inference <- function(
    pair_tbl,
    score_direction,
    pvalue_semantics,
    p_adjust_method
) {
  oriented <- .orient_pair_scores(
    up = pair_tbl$NES_Up,
    down = pair_tbl$NES_Down,
    score_direction = score_direction
  )

  up_positive <- .directional_arm_p(
    p = pair_tbl$pvalue_Up,
    oriented_score = oriented[, 1],
    alternative = "positive",
    pvalue_semantics = pvalue_semantics
  )
  down_positive <- .directional_arm_p(
    p = pair_tbl$pvalue_Down,
    oriented_score = oriented[, 2],
    alternative = "positive",
    pvalue_semantics = pvalue_semantics
  )
  up_negative <- .directional_arm_p(
    p = pair_tbl$pvalue_Up,
    oriented_score = oriented[, 1],
    alternative = "negative",
    pvalue_semantics = pvalue_semantics
  )
  down_negative <- .directional_arm_p(
    p = pair_tbl$pvalue_Down,
    oriented_score = oriented[, 2],
    alternative = "negative",
    pvalue_semantics = pvalue_semantics
  )

  p_activation <- pmax(up_positive, down_positive)
  p_inhibition <- pmax(up_negative, down_negative)
  adjusted <- .adjust_complete_p(
    c(p_activation, p_inhibition),
    method = p_adjust_method
  )
  n_pairs <- nrow(pair_tbl)
  padj_activation <- adjusted[seq_len(n_pairs)]
  padj_inhibition <- adjusted[n_pairs + seq_len(n_pairs)]

  missing_score <- is.na(oriented[, 1]) | is.na(oriented[, 2])
  # Explicit completeness masks keep NA out of vector assignment indices.
  # A missing arm is unknown, not an observed zero or a discordant pair.
  activated <- !missing_score & oriented[, 1] > 0 & oriented[, 2] > 0
  inhibited <- !missing_score & oriented[, 1] < 0 & oriented[, 2] < 0
  incomplete <- !missing_score & (oriented[, 1] == 0 | oriented[, 2] == 0)

  directional_state <- rep("Discordant", n_pairs)
  directional_state[activated] <- "Activated"
  directional_state[inhibited] <- "Inhibited"
  directional_state[incomplete] <- "Incomplete"
  directional_state[missing_score] <- NA_character_

  p_coherent <- rep(1, n_pairs)
  padj_coherent <- rep(1, n_pairs)
  p_coherent[activated] <- p_activation[activated]
  p_coherent[inhibited] <- p_inhibition[inhibited]
  padj_coherent[activated] <- padj_activation[activated]
  padj_coherent[inhibited] <- padj_inhibition[inhibited]
  p_coherent[missing_score] <- NA_real_
  padj_coherent[missing_score] <- NA_real_

  pair_tbl$p_activation <- p_activation
  pair_tbl$p_inhibition <- p_inhibition
  pair_tbl$padj_activation <- padj_activation
  pair_tbl$padj_inhibition <- padj_inhibition
  pair_tbl$p_coherent <- p_coherent
  pair_tbl$padj_coherent <- padj_coherent
  pair_tbl$directional_state <- directional_state
  pair_tbl$pvalue_semantics <- rep(pvalue_semantics, n_pairs)
  pair_tbl
}

.directional_arm_p <- function(
    p,
    oriented_score,
    alternative = c("positive", "negative"),
    pvalue_semantics = c("observed_tail", "two_sided")
) {
  alternative <- match.arg(alternative)
  pvalue_semantics <- match.arg(pvalue_semantics)
  out <- rep(NA_real_, length(p))
  complete <- !is.na(p) & !is.na(oriented_score)

  if (pvalue_semantics == "observed_tail") {
    supported <- if (alternative == "positive") {
      oriented_score > 0
    } else {
      oriented_score < 0
    }
    out[complete] <- 1
    out[complete & supported] <- p[complete & supported]
    return(out)
  }

  positive_p <- ifelse(
    oriented_score > 0,
    p / 2,
    ifelse(oriented_score < 0, 1 - p / 2, 0.5)
  )
  out[complete] <- if (alternative == "positive") {
    positive_p[complete]
  } else {
    1 - positive_p[complete]
  }
  out
}
