#' Score individual samples with L1000PAI
#'
#' Compute perturbation activity for each sample separately, so that the same
#' estimator applies to a cohort and to a single tumour. Each sample is scored
#' from its own within-sample gene ranking, and the zero point is set by a
#' size-matched random-set null computed inside that same sample. No other
#' sample is required. Scoring a sample alone or in a cohort requires the same
#' fixed gene universe, identifier mapping, preprocessing, library and settings
#' (including the seed for null-based results). This does not establish
#' cross-platform or cross-batch biological comparability.
#'
#' The paired Up/Down structure of an L1000 library is preserved: an Up and a
#' Down set score are computed separately and combined with the same direction
#' convention used by [collapse_pairs()], resolved through
#' [resolve_score_direction()]. For `chem`-style libraries the activity is
#' `Up - Down`; for `crispr` it is `Down - Up`.
#'
#' @section Choice of zero point:
#' `reference = "null"` (the default) expresses each score in units of the
#' size-matched random-set null within the same sample, so zero means "no
#' enrichment". `reference = "library_median"` centres instead on the median
#' score across all perturbagens in that sample, so zero means "the typical
#' perturbagen in this sample"; it is self-contained but inherits whatever the
#' library composition happens to be. `reference = "none"` returns the raw
#' statistic. Note that any purely within-sample centring removes only a
#' per-sample scalar and therefore cannot remove structured covariation between
#' perturbagens; that has to be handled downstream (for example by testing at
#' drug-class level, or by adjusting for confounders).
#'
#' @param expr A numeric matrix of expression values with genes in rows and
#'   samples in columns. `rownames(expr)` must be the gene identifiers used by
#'   the library. Values should already be on a continuous scale such as
#'   log2 intensity or log-CPM. Gene identifiers must be unique and non-empty;
#'   resolve duplicate probes upstream using a fixed mapping independent of
#'   the samples being scored. Constant sample columns are not informative and
#'   are rejected. A single-column matrix is allowed.
#' @param library Optional L1000 library name, one of
#'   `available_l1000_libraries()`. Ignored when `TERM2GENE` is supplied.
#' @param TERM2GENE Optional two-column data.frame with columns `term` and
#'   `gene`, as returned by [get_l1000_library()] or [read_l1000_gmt()].
#' @param method Per-sample set statistic. `"rank_mean"` (the default) uses the
#'   mean within-sample rank percentile of the set members, which is fast and
#'   robust. `"ssgsea"` uses the weighted single-sample enrichment statistic of
#'   Barbie et al., integrated over the walk, with weights
#'   `rank(expression, ties.method = "average")^tau`. Expression ties in the
#'   decreasing walk are ordered by gene identifier. This is the unnormalized
#'   integrated statistic; no across-sample range normalization is applied.
#'   Every scored arm must leave at least one gene outside the set.
#' @param reference Definition of the zero point. One of `"null"`,
#'   `"library_median"` or `"none"`. See Details.
#' @param score_direction One of `"auto"`, `"down_minus_up"` or
#'   `"up_minus_down"`, passed to [resolve_score_direction()].
#' @param tau Non-negative exponent applied to the within-sample ranks when
#'   `method = "ssgsea"`. Ignored otherwise.
#' @param n_null Number of size-matched random sets drawn per size bin, per
#'   sample. Set to `0` to skip the null, in which case no empirical p-values
#'   are returned and `reference = "null"` is not allowed.
#' @param null_size_bin Set sizes are rounded to a multiple of this value before
#'   a null is drawn, so that nearby set sizes share one null. L1000 signature
#'   sets are near-uniform in size, which makes this cheap and accurate.
#' @param min_size Minimum number of set members that must be present in
#'   `rownames(expr)` for a direction to be scored.
#' @param max_size Optional maximum number of present members per direction.
#' @param rng_seed Optional integer seed. When supplied it is set immediately
#'   before the null is drawn, making the result reproducible.
#' @param verbose Logical; emit progress messages.
#'
#' @return An object of class `l1000pai_sample_scores`: a list with matrices
#'   `score`, `z`, `p_activated` and `p_inhibited` (perturbagens in rows,
#'   samples in columns), a `sets` data.frame describing the scored
#'   perturbagens, a `null_summary` data.frame, and a `meta` list recording the
#'   call parameters. Use [sample_score_table()] to obtain a long data.frame.
#'
#' @seealso [sample_score_table()], [collapse_pairs()], [run_l1000_pai()]
#' @export
score_samples_l1000pai <- function(
    expr,
    library = NULL,
    TERM2GENE = NULL,
    method = c("rank_mean", "ssgsea"),
    reference = c("null", "library_median", "none"),
    score_direction = c("auto", "down_minus_up", "up_minus_down"),
    tau = 0.25,
    n_null = 1000L,
    null_size_bin = 5L,
    min_size = 10L,
    max_size = NULL,
    rng_seed = NULL,
    verbose = TRUE
) {
  method <- match.arg(method)
  reference <- match.arg(reference)
  score_direction <- match.arg(score_direction)

  # ---- validate expr -------------------------------------------------------
  if (!is.matrix(expr) || !is.numeric(expr)) {
    stop("`expr` must be a numeric matrix with genes in rows and samples in columns.")
  }
  if (nrow(expr) < 2L) {
    stop("`expr` must contain at least two genes.")
  }
  if (is.null(rownames(expr))) {
    stop("`expr` must have rownames giving gene identifiers.")
  }
  if (anyNA(rownames(expr)) || any(!nzchar(trimws(rownames(expr))))) {
    stop("`expr` gene identifiers must be non-missing and non-empty.")
  }
  if (ncol(expr) < 1L) {
    stop("`expr` must have at least one sample column.")
  }
  if (is.null(colnames(expr))) {
    colnames(expr) <- paste0("sample", seq_len(ncol(expr)))
  }
  if (anyNA(expr)) {
    stop("`expr` must not contain NA values; impute or filter them first.")
  }
  if (!all(is.finite(expr))) {
    stop("`expr` must contain finite values only.")
  }

  # A cohort-dependent probe choice would change an existing sample's score.
  if (anyDuplicated(rownames(expr))) {
    stop(
      "`expr` has duplicate gene identifiers. Resolve duplicate probes upstream ",
      "using a fixed, sample-independent mapping before scoring."
    )
  }
  # Canonical order also makes seeded random sets invariant to input row order.
  expr <- expr[order(rownames(expr), method = "radix"), , drop = FALSE]
  constant <- apply(expr, 2L, function(v) all(v == v[1L]))
  if (any(constant)) {
    stop(
      "`expr` contains a constant sample column with no informative gene ranking: ",
      paste(colnames(expr)[constant], collapse = ", "), "."
    )
  }

  if (!is.numeric(tau) || length(tau) != 1L || !is.finite(tau) || tau < 0) {
    stop("`tau` must be a single non-negative number.")
  }
  if (!is.numeric(n_null) || length(n_null) != 1L || !is.finite(n_null) ||
      n_null < 0 || n_null > .Machine$integer.max || n_null != floor(n_null)) {
    stop("`n_null` must be a single non-negative integer no greater than .Machine$integer.max.")
  }
  n_null <- as.integer(n_null)
  if (identical(reference, "null") && n_null == 0L) {
    stop("`reference = \"null\"` requires `n_null` greater than 0.")
  }
  if (!is.numeric(null_size_bin) || length(null_size_bin) != 1L ||
      !is.finite(null_size_bin) || null_size_bin < 1 ||
      null_size_bin > .Machine$integer.max || null_size_bin != floor(null_size_bin)) {
    stop("`null_size_bin` must be a single positive integer no greater than .Machine$integer.max.")
  }
  null_size_bin <- as.integer(null_size_bin)
  if (!is.numeric(min_size) || length(min_size) != 1L || !is.finite(min_size) ||
      min_size < 1 || min_size > .Machine$integer.max || min_size != floor(min_size)) {
    stop("`min_size` must be a single positive integer no greater than .Machine$integer.max.")
  }
  min_size <- as.integer(min_size)
  if (!is.null(max_size)) {
    if (!is.numeric(max_size) || length(max_size) != 1L || !is.finite(max_size) ||
        max_size < min_size || max_size > .Machine$integer.max ||
        max_size != floor(max_size)) {
      stop("`max_size` must be NULL or a single integer between `min_size` and .Machine$integer.max.")
    }
    max_size <- as.integer(max_size)
  }
  if (!is.null(rng_seed)) {
    if (!is.numeric(rng_seed) || length(rng_seed) != 1L ||
        !is.finite(rng_seed) || rng_seed < 0 ||
        rng_seed > .Machine$integer.max || rng_seed != floor(rng_seed)) {
      stop("`rng_seed` must be NULL or a single integer between 0 and .Machine$integer.max.")
    }
    rng_seed <- as.integer(rng_seed)
  }

  # ---- library -------------------------------------------------------------
  if (is.null(TERM2GENE)) {
    if (is.null(library)) {
      stop("Supply either `library` or `TERM2GENE`.")
    }
    TERM2GENE <- get_l1000_library(library)
  }
  if (!is.data.frame(TERM2GENE) ||
      !all(c("term", "gene") %in% colnames(TERM2GENE))) {
    stop("`TERM2GENE` must be a data.frame with columns `term` and `gene`.")
  }
  resolved_direction <- resolve_score_direction(library, score_direction)

  pairs <- .l1000pai_pair_sets(
    TERM2GENE = TERM2GENE,
    universe = rownames(expr),
    min_size = min_size,
    max_size = max_size
  )
  if (length(pairs$up) == 0L) {
    stop(
      "No Up/Down signature pair had at least `min_size` members present in ",
      "`rownames(expr)`. Check that the gene identifiers match the library."
    )
  }
  perturbagens <- names(pairs$up)
  n_pert <- length(perturbagens)
  n_samp <- ncol(expr)
  if (verbose) {
    message(
      "Scoring ", n_samp, " sample(s) against ", n_pert,
      " perturbagen pair(s) [method = ", method,
      ", direction = ", resolved_direction, "]."
    )
  }

  # ---- per-sample precomputation ------------------------------------------
  prep <- .l1000pai_sample_prep(expr, method = method, tau = tau)

  up_stat <- matrix(
    NA_real_, nrow = n_pert, ncol = n_samp,
    dimnames = list(perturbagens, colnames(expr))
  )
  dn_stat <- up_stat
  for (i in seq_len(n_pert)) {
    up_stat[i, ] <- .l1000pai_set_stat(prep, pairs$up[[i]], method = method)
    dn_stat[i, ] <- .l1000pai_set_stat(prep, pairs$dn[[i]], method = method)
  }

  score <- if (identical(resolved_direction, "up_minus_down")) {
    up_stat - dn_stat
  } else {
    dn_stat - up_stat
  }

  # ---- size-matched random-set null ---------------------------------------
  null_median <- matrix(0, nrow = n_pert, ncol = n_samp,
                        dimnames = dimnames(score))
  null_scale <- matrix(1, nrow = n_pert, ncol = n_samp,
                       dimnames = dimnames(score))
  p_act <- matrix(NA_real_, nrow = n_pert, ncol = n_samp, dimnames = dimnames(score))
  p_inh <- p_act
  null_summary <- NULL

  if (n_null > 0L) {
    if (!is.null(rng_seed)) {
      set.seed(rng_seed)
    }
    n_up <- vapply(pairs$up, length, integer(1L))
    n_dn <- vapply(pairs$dn, length, integer(1L))
    bin_up <- pmax(null_size_bin, round(n_up / null_size_bin) * null_size_bin)
    bin_dn <- pmax(null_size_bin, round(n_dn / null_size_bin) * null_size_bin)
    bin_key <- paste(bin_up, bin_dn, sep = "x")
    keys <- unique(bin_key)
    if (verbose) {
      message(
        "Drawing ", n_null, " size-matched random sets for ",
        length(keys), " size bin(s)."
      )
    }
    universe <- rownames(expr)
    n_genes <- length(universe)
    summary_rows <- vector("list", length(keys))

    for (k in seq_along(keys)) {
      key <- keys[[k]]
      idx <- which(bin_key == key)
      size_up <- bin_up[idx[1L]]
      size_dn <- bin_dn[idx[1L]]
      if (size_up + size_dn > n_genes) {
        stop(
          "Size bin ", key, " requires more genes than are present in `expr`."
        )
      }
      if (verbose) {
        message(
          "  [", k, "/", length(keys), "] size bin ", key,
          " covering ", length(idx), " perturbagen(s)."
        )
      }
      null_draws <- matrix(NA_real_, nrow = n_null, ncol = n_samp)
      for (b in seq_len(n_null)) {
        picked <- sample.int(n_genes, size_up + size_dn, replace = FALSE)
        u <- picked[seq_len(size_up)]
        d <- picked[(size_up + 1L):(size_up + size_dn)]
        su <- .l1000pai_set_stat(prep, u, method = method)
        sd_ <- .l1000pai_set_stat(prep, d, method = method)
        null_draws[b, ] <- if (identical(resolved_direction, "up_minus_down")) {
          su - sd_
        } else {
          sd_ - su
        }
      }
      med <- apply(null_draws, 2L, stats::median)
      scl <- apply(null_draws, 2L, stats::mad)
      # MAD 為 0 時退回標準差，再退回極小值，避免除以 0
      bad <- !is.finite(scl) | scl <= 0
      if (any(bad)) {
        alt <- apply(null_draws[, bad, drop = FALSE], 2L, stats::sd)
        alt[!is.finite(alt) | alt <= 0] <- .Machine$double.eps
        scl[bad] <- alt
      }
      for (i in idx) {
        null_median[i, ] <- med
        null_scale[i, ] <- scl
        obs <- score[i, ]
        ge <- colSums(null_draws >= matrix(obs, nrow = n_null, ncol = n_samp, byrow = TRUE))
        le <- colSums(null_draws <= matrix(obs, nrow = n_null, ncol = n_samp, byrow = TRUE))
        p_act[i, ] <- (1 + ge) / (n_null + 1)
        p_inh[i, ] <- (1 + le) / (n_null + 1)
      }
      summary_rows[[k]] <- data.frame(
        size_bin = key,
        size_up = size_up,
        size_dn = size_dn,
        n_perturbagens = length(idx),
        null_median_mean = mean(med),
        null_scale_mean = mean(scl),
        stringsAsFactors = FALSE
      )
    }
    null_summary <- do.call(rbind, summary_rows)
  }

  # ---- zero point ----------------------------------------------------------
  z <- switch(
    reference,
    null = (score - null_median) / null_scale,
    library_median = {
      med <- apply(score, 2L, stats::median)
      scl <- apply(score, 2L, stats::mad)
      bad <- !is.finite(scl) | scl <= 0
      if (any(bad)) {
        alt <- apply(score[, bad, drop = FALSE], 2L, stats::sd)
        alt[!is.finite(alt) | alt <= 0] <- .Machine$double.eps
        scl[bad] <- alt
      }
      sweep(sweep(score, 2L, med, "-"), 2L, scl, "/")
    },
    none = score
  )

  sets <- data.frame(
    perturbagen = perturbagens,
    n_up = vapply(pairs$up, length, integer(1L)),
    n_dn = vapply(pairs$dn, length, integer(1L)),
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  out <- list(
    score = score,
    z = z,
    p_activated = p_act,
    p_inhibited = p_inh,
    sets = sets,
    null_summary = null_summary,
    meta = list(
      method = method,
      reference = reference,
      score_direction = resolved_direction,
      library = library,
      tau = if (identical(method, "ssgsea")) tau else NA_real_,
      n_null = n_null,
      null_size_bin = null_size_bin,
      min_size = min_size,
      max_size = max_size,
      n_genes = nrow(expr),
      n_samples = n_samp,
      n_perturbagens = n_pert,
      rng_seed = rng_seed
    )
  )
  class(out) <- "l1000pai_sample_scores"
  out
}


#' Convert sample-level scores to a long table
#'
#' @param x An object of class `l1000pai_sample_scores`.
#' @param what Which matrix to melt. One of `"z"`, `"score"`, `"p_activated"`
#'   or `"p_inhibited"`. When `"all"`, every matrix is returned as one column
#'   each.
#' @param perturbagens Optional character vector restricting the output.
#'
#' @return A data.frame with columns `perturbagen`, `sample` and the requested
#'   value column(s).
#' @export
sample_score_table <- function(
    x,
    what = c("z", "score", "p_activated", "p_inhibited", "all"),
    perturbagens = NULL
) {
  what <- match.arg(what)
  if (!inherits(x, "l1000pai_sample_scores")) {
    stop("`x` must be an object of class `l1000pai_sample_scores`.")
  }
  keep <- rownames(x$score)
  if (!is.null(perturbagens)) {
    missing_pert <- setdiff(perturbagens, keep)
    if (length(missing_pert) > 0L) {
      stop(
        "Unknown perturbagen(s): ",
        paste(utils::head(missing_pert, 10L), collapse = ", "),
        if (length(missing_pert) > 10L) ", ..." else "."
      )
    }
    keep <- perturbagens
  }

  fields <- if (identical(what, "all")) {
    c("z", "score", "p_activated", "p_inhibited")
  } else {
    what
  }

  base <- x$score[keep, , drop = FALSE]
  out <- data.frame(
    perturbagen = rep(rownames(base), times = ncol(base)),
    sample = rep(colnames(base), each = nrow(base)),
    stringsAsFactors = FALSE
  )
  for (f in fields) {
    m <- x[[f]][keep, , drop = FALSE]
    out[[f]] <- as.vector(m)
  }
  out
}


#' @export
print.l1000pai_sample_scores <- function(x, ...) {
  m <- x$meta
  cat("<l1000pai_sample_scores>\n")
  cat("  perturbagens : ", m$n_perturbagens, "\n", sep = "")
  cat("  samples      : ", m$n_samples, "\n", sep = "")
  cat("  method       : ", m$method, "\n", sep = "")
  cat("  reference    : ", m$reference, "\n", sep = "")
  cat("  direction    : ", m$score_direction, "\n", sep = "")
  cat("  null draws   : ", m$n_null, "\n", sep = "")
  invisible(x)
}


# ---------------------------------------------------------------- internals

#' Build paired Up/Down member index lists restricted to a gene universe
#'
#' @param TERM2GENE Two-column data.frame with `term` and `gene`.
#' @param universe Character vector of measured gene identifiers.
#' @param min_size Minimum members present per direction.
#' @param max_size Optional maximum members present per direction.
#'
#' @return A list with `up` and `dn`, each a named list of integer positions
#'   into `universe`.
#' @keywords internal
#' @noRd
.l1000pai_pair_sets <- function(TERM2GENE, universe, min_size, max_size = NULL) {
  term <- as.character(TERM2GENE$term)
  gene <- as.character(TERM2GENE$gene)

  # 與 collapse_pairs() 使用相同的方向後綴規則
  is_up <- grepl("([[:space:]_.-]+)(Up|up)$", term)
  is_dn <- grepl("([[:space:]_.-]+)(Down|down|Dn|dn)$", term)
  if (!any(is_up) || !any(is_dn)) {
    stop("No paired Up/Down signatures were found in `TERM2GENE$term`.")
  }

  perturbagen <- term
  perturbagen <- sub("([[:space:]_.-]+)(Up|up)$", "", perturbagen)
  perturbagen <- sub("([[:space:]_.-]+)(Down|down|Dn|dn)$", "", perturbagen)
  perturbagen <- trimws(perturbagen)

  pos <- match(gene, universe)
  ok <- !is.na(pos) & (is_up | is_dn)

  direction <- ifelse(is_up, "up", ifelse(is_dn, "dn", NA_character_))
  df <- data.frame(
    perturbagen = perturbagen[ok],
    direction = direction[ok],
    pos = pos[ok],
    stringsAsFactors = FALSE
  )
  df <- df[!duplicated(df[, c("perturbagen", "direction", "pos")]), , drop = FALSE]

  split_up <- split(df$pos[df$direction == "up"], df$perturbagen[df$direction == "up"])
  split_dn <- split(df$pos[df$direction == "dn"], df$perturbagen[df$direction == "dn"])

  common <- intersect(names(split_up), names(split_dn))
  if (length(common) == 0L) {
    stop("No perturbagen had both an Up and a Down signature.")
  }
  common <- sort(common)
  split_up <- split_up[common]
  split_dn <- split_dn[common]

  size_ok <- vapply(common, function(p) {
    nu <- length(split_up[[p]])
    nd <- length(split_dn[[p]])
    big_enough <- nu >= min_size && nd >= min_size
    small_enough <- is.null(max_size) || (nu <= max_size && nd <= max_size)
    big_enough && small_enough
  }, logical(1L))

  list(up = split_up[size_ok], dn = split_dn[size_ok])
}


#' Per-sample precomputation shared by all sets
#'
#' @param expr Numeric expression matrix.
#' @param method `"rank_mean"` or `"ssgsea"`.
#' @param tau Weight exponent for `"ssgsea"`.
#'
#' @return A list of precomputed per-sample quantities.
#' @keywords internal
#' @noRd
.l1000pai_sample_prep <- function(expr, method, tau) {
  n_genes <- nrow(expr)
  if (identical(method, "rank_mean")) {
    rank_pct <- apply(expr, 2L, function(v) rank(v, ties.method = "average") / n_genes)
    rank_pct <- matrix(rank_pct, nrow = n_genes, dimnames = dimnames(expr))
    return(list(method = method, n_genes = n_genes, rank_pct = rank_pct))
  }

  # Integral coefficient: N - decreasing-walk position + 1. Resolve ties by
  # gene identifier, independently of the incoming expression row order.
  rev_rank <- apply(expr, 2L, function(v) {
    ord <- order(-v, rownames(expr), method = "radix")
    p <- integer(n_genes)
    p[ord] <- seq_len(n_genes)
    n_genes - p + 1L
  })
  rev_rank <- matrix(as.numeric(rev_rank), nrow = n_genes, dimnames = dimnames(expr))
  sample_ranks <- apply(expr, 2L, rank, ties.method = "average")
  sample_ranks <- matrix(sample_ranks, nrow = n_genes, dimnames = dimnames(expr))
  w <- sample_ranks^tau
  wr <- w * rev_rank
  if (any(!is.finite(wr))) {
    stop("`tau` produces non-finite ssGSEA rank weights; use a smaller exponent.")
  }
  list(
    method = method,
    n_genes = n_genes,
    rev_rank = rev_rank,
    w = w,
    wr = wr,
    total = n_genes * (n_genes + 1) / 2
  )
}


#' Statistic for one gene set across all samples
#'
#' For `"ssgsea"` this uses the closed form of the integrated weighted walk,
#' which is algebraically identical to accumulating `P_w - P_ns` over the
#' ranking but costs `O(set size)` instead of `O(genes)` per set.
#'
#' @param prep Output of `.l1000pai_sample_prep()`.
#' @param idx Integer positions of the set members.
#' @param method `"rank_mean"` or `"ssgsea"`.
#'
#' @return A numeric vector with one value per sample.
#' @keywords internal
#' @noRd
.l1000pai_set_stat <- function(prep, idx, method) {
  if (length(idx) == 0L) {
    stop("A scored gene set must contain at least one gene.")
  }
  if (identical(method, "rank_mean")) {
    return(colMeans(prep$rank_pct[idx, , drop = FALSE]))
  }
  ns <- length(idx)
  if (ns >= prep$n_genes) {
    stop("ssGSEA requires a non-empty background: each gene set must contain fewer genes than the expression universe.")
  }
  num <- colSums(prep$wr[idx, , drop = FALSE])
  den <- colSums(prep$w[idx, , drop = FALSE])
  in_sum <- colSums(prep$rev_rank[idx, , drop = FALSE])
  stat <- num / den - (prep$total - in_sum) / (prep$n_genes - ns)
  if (any(!is.finite(stat))) {
    stop("ssGSEA could not produce a finite statistic; check set sizes and `tau`.")
  }
  stat
}
