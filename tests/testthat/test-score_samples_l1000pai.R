test_that("paired sets are built from every direction-suffix spelling", {
  universe <- paste0("G", 1:60)
  t2g <- do.call(rbind, lapply(
    c("Fostamatinib up", "Fostamatinib dn",
      "AS-605240 up", "AS-605240 dn",
      "DrugX_Up", "DrugX_Down",
      "DrugY-Up", "DrugY-Dn",
      "Lonely up", "NoDirection"),
    function(tm) data.frame(term = tm, gene = paste0("G", 1:20), stringsAsFactors = FALSE)
  ))

  pairs <- L1000PAI:::.l1000pai_pair_sets(t2g, universe, min_size = 5L)

  expect_setequal(
    names(pairs$up),
    c("AS-605240", "DrugX", "DrugY", "Fostamatinib")
  )
  # 只有單一方向或沒有方向後綴的 term 不應成對
  expect_false("Lonely" %in% names(pairs$up))
  expect_false("NoDirection" %in% names(pairs$up))
  expect_equal(names(pairs$up), names(pairs$dn))
})


test_that("min_size and max_size filter directions by present members", {
  universe <- paste0("G", 1:100)
  t2g <- rbind(
    data.frame(term = "A up", gene = paste0("G", 1:30), stringsAsFactors = FALSE),
    data.frame(term = "A dn", gene = paste0("G", 31:60), stringsAsFactors = FALSE),
    data.frame(term = "B up", gene = paste0("G", 1:4), stringsAsFactors = FALSE),
    data.frame(term = "B dn", gene = paste0("G", 5:8), stringsAsFactors = FALSE)
  )

  keep <- L1000PAI:::.l1000pai_pair_sets(t2g, universe, min_size = 10L)
  expect_equal(names(keep$up), "A")

  capped <- L1000PAI:::.l1000pai_pair_sets(t2g, universe, min_size = 2L, max_size = 10L)
  expect_equal(names(capped$up), "B")
})


test_that("ssgsea equals an independent rank-weighted integrated walk with ties", {
  set.seed(42)
  n_genes <- 400L
  continuous <- stats::rnorm(n_genes)
  genes <- paste0("G", seq_len(n_genes))
  idx <- sort(sample.int(n_genes, 40L))

  # The reference accumulates the walk directly; it does not use the closed
  # form or any package preprocessing helper. Rounding introduces ties.
  for (x in list(continuous, round(continuous), rep(c(-3, 0, 2, 5), 100L))) {
    expr <- matrix(x, ncol = 1L, dimnames = list(genes, "s1"))
    for (tau in c(0, 0.25, 1)) {
      ranks <- rank(x, ties.method = "average")
      ord <- order(-x, genes, method = "radix")
      in_set <- ord %in% idx
      w <- ranks[ord]^tau
      p_w <- cumsum(ifelse(in_set, w, 0)) / sum(w[in_set])
      p_ns <- cumsum(!in_set) / (n_genes - sum(in_set))
      slow <- sum(p_w - p_ns)

      prep <- L1000PAI:::.l1000pai_sample_prep(expr, method = "ssgsea", tau = tau)
      fast <- L1000PAI:::.l1000pai_set_stat(prep, idx, method = "ssgsea")
      expect_equal(unname(fast), slow, tolerance = 1e-10)
    }
  }
})


test_that("rank_mean is the mean within-sample rank percentile", {
  expr <- matrix(
    c(5, 1, 3, 2, 4,
      1, 2, 3, 4, 5),
    nrow = 5L,
    dimnames = list(paste0("G", 1:5), c("a", "b"))
  )
  prep <- L1000PAI:::.l1000pai_sample_prep(expr, method = "rank_mean", tau = 0.25)
  got <- L1000PAI:::.l1000pai_set_stat(prep, c(1L, 5L), method = "rank_mean")

  expected_a <- mean(c(rank(expr[, "a"])[1], rank(expr[, "a"])[5]) / 5)
  expected_b <- mean(c(rank(expr[, "b"])[1], rank(expr[, "b"])[5]) / 5)
  expect_equal(unname(got), c(expected_a, expected_b))
})


make_toy <- function(n_genes = 800L, n_samples = 3L, planted = TRUE, seed = 7L) {
  set.seed(seed)
  expr <- matrix(
    stats::rnorm(n_genes * n_samples),
    nrow = n_genes,
    dimnames = list(paste0("G", seq_len(n_genes)),
                    paste0("s", seq_len(n_samples)))
  )
  up_genes <- paste0("G", 1:60)
  dn_genes <- paste0("G", 61:120)
  if (planted) {
    expr[up_genes, 1L] <- expr[up_genes, 1L] + 2
    expr[dn_genes, 1L] <- expr[dn_genes, 1L] - 2
  }
  t2g <- rbind(
    data.frame(term = "DrugA up", gene = up_genes, stringsAsFactors = FALSE),
    data.frame(term = "DrugA dn", gene = dn_genes, stringsAsFactors = FALSE),
    data.frame(term = "DrugB up", gene = paste0("G", 201:260), stringsAsFactors = FALSE),
    data.frame(term = "DrugB dn", gene = paste0("G", 261:320), stringsAsFactors = FALSE)
  )
  list(expr = expr, t2g = t2g)
}


test_that("a planted activation is recovered in the planted sample only", {
  toy <- make_toy()
  res <- score_samples_l1000pai(
    toy$expr,
    TERM2GENE = toy$t2g,
    score_direction = "up_minus_down",
    method = "rank_mean",
    n_null = 400L,
    rng_seed = 1L,
    verbose = FALSE
  )

  expect_s3_class(res, "l1000pai_sample_scores")
  expect_true(res$z["DrugA", "s1"] > 5)
  expect_lt(res$p_activated["DrugA", "s1"], 0.01)

  # 沒有植入訊號的樣本與對照藥不應被叫出來
  expect_true(all(abs(res$z["DrugA", c("s2", "s3")]) < 5))
  expect_gt(res$p_activated["DrugB", "s1"], 0.01)
})


test_that("a single sample can be scored on its own", {
  toy <- make_toy(n_samples = 1L)
  res <- score_samples_l1000pai(
    toy$expr,
    TERM2GENE = toy$t2g,
    score_direction = "up_minus_down",
    n_null = 200L,
    rng_seed = 1L,
    verbose = FALSE
  )

  expect_equal(ncol(res$score), 1L)
  expect_equal(ncol(res$z), 1L)
  expect_true(is.finite(res$z["DrugA", 1L]))
  expect_true(res$z["DrugA", 1L] > 5)
})


test_that("score_direction flips the sign of the paired statistic", {
  toy <- make_toy()
  up_first <- score_samples_l1000pai(
    toy$expr, TERM2GENE = toy$t2g, score_direction = "up_minus_down",
    n_null = 0L, reference = "none", verbose = FALSE
  )
  dn_first <- score_samples_l1000pai(
    toy$expr, TERM2GENE = toy$t2g, score_direction = "down_minus_up",
    n_null = 0L, reference = "none", verbose = FALSE
  )
  expect_equal(up_first$score, -dn_first$score)
})


test_that("library_median reference sets the per-sample median to zero", {
  toy <- make_toy()
  res <- score_samples_l1000pai(
    toy$expr, TERM2GENE = toy$t2g, reference = "library_median",
    n_null = 0L, verbose = FALSE
  )
  expect_equal(unname(apply(res$z, 2L, stats::median)), rep(0, ncol(res$z)))
})


test_that("results are reproducible for a fixed rng_seed", {
  toy <- make_toy()
  args <- list(
    toy$expr, TERM2GENE = toy$t2g, score_direction = "up_minus_down",
    n_null = 200L, rng_seed = 99L, verbose = FALSE
  )
  a <- do.call(score_samples_l1000pai, args)
  b <- do.call(score_samples_l1000pai, args)
  expect_equal(a$z, b$z)
  expect_equal(a$p_activated, b$p_activated)
})


test_that("sample_score_table melts the requested matrices", {
  toy <- make_toy()
  res <- score_samples_l1000pai(
    toy$expr, TERM2GENE = toy$t2g, n_null = 100L, rng_seed = 5L, verbose = FALSE
  )

  tbl <- sample_score_table(res, what = "z")
  expect_equal(nrow(tbl), nrow(res$z) * ncol(res$z))
  expect_named(tbl, c("perturbagen", "sample", "z"))

  all_tbl <- sample_score_table(res, what = "all")
  expect_true(all(c("z", "score", "p_activated", "p_inhibited") %in% names(all_tbl)))

  one <- sample_score_table(res, what = "score", perturbagens = "DrugA")
  expect_equal(unique(one$perturbagen), "DrugA")

  expect_error(sample_score_table(res, perturbagens = "NoSuchDrug"), "Unknown perturbagen")
})


test_that("invalid input is rejected", {
  toy <- make_toy()

  expect_error(score_samples_l1000pai(as.data.frame(toy$expr), TERM2GENE = toy$t2g),
               "must be a numeric matrix")

  no_names <- toy$expr
  rownames(no_names) <- NULL
  expect_error(score_samples_l1000pai(no_names, TERM2GENE = toy$t2g), "must have rownames")

  with_na <- toy$expr
  with_na[1L, 1L] <- NA_real_
  expect_error(score_samples_l1000pai(with_na, TERM2GENE = toy$t2g), "must not contain NA")

  expect_error(score_samples_l1000pai(toy$expr), "Supply either")

  expect_error(
    score_samples_l1000pai(toy$expr, TERM2GENE = toy$t2g, reference = "null", n_null = 0L),
    "requires `n_null` greater than 0"
  )

  expect_error(
    score_samples_l1000pai(toy$expr, TERM2GENE = toy$t2g, rng_seed = -1),
    "`rng_seed` must be NULL"
  )

  unpaired <- data.frame(
    term = "DrugA up", gene = paste0("G", 1:60), stringsAsFactors = FALSE
  )
  expect_error(
    score_samples_l1000pai(toy$expr, TERM2GENE = unpaired, verbose = FALSE),
    "No paired Up/Down signatures"
  )
})


test_that("duplicate gene identifiers cannot depend on the scoring cohort", {
  # The previous single/cohort duplicate rules gave +0.5 and -0.5 for s1.
  expr <- matrix(
    c(10, 0, 5, 5, 5, 5, 10, 100, 5, 5, 5, 5), nrow = 6L,
    dimnames = list(c("A", "A", "B", "C", "D", "E"), c("s1", "s2")))
  t2g <- data.frame(term = c("Drug up", "Drug dn"), gene = c("A", "B"))
  for (method in c("rank_mean", "ssgsea")) {
    for (x in list(expr[, "s1", drop = FALSE], expr)) {
      expect_error(
        score_samples_l1000pai(
          x, TERM2GENE = t2g, method = method, reference = "none", n_null = 0L,
          min_size = 1L, verbose = FALSE),
        "duplicate gene identifiers.*fixed, sample-independent mapping"
      )
    }
  }
})


test_that("ssgsea is unchanged by expression shifts or positive scaling", {
  # Previously adding 1000 changed the sign despite identical gene ranks.
  expr <- matrix(1:10, ncol = 1L, dimnames = list(paste0("G", 1:10), "s1"))
  t2g <- data.frame(
    term = rep(c("Drug up", "Drug dn"), each = 2L),
    gene = paste0("G", c(1, 9, 4, 7)))
  score <- function(x) score_samples_l1000pai(
    x, TERM2GENE = t2g, method = "ssgsea", reference = "null", n_null = 50L,
    null_size_bin = 1L, min_size = 2L, score_direction = "up_minus_down",
    rng_seed = 91L, verbose = FALSE)
  base <- score(expr)
  for (shifted in list(expr + 1000, expr - 500, expr * 3)) {
    got <- score(shifted)
    for (field in c("score", "z", "p_activated", "p_inhibited")) {
      expect_equal(got[[field]], base[[field]])
    }
  }
  expect_gt(unname(base$score[1, 1]), 0)
})


test_that("a unique-gene sample has identical seeded results alone and in a cohort", {
  toy <- make_toy()
  # Include ties and test a sample that is not the cohort's first column.
  toy$expr <- round(toy$expr, 1)
  for (method in c("rank_mean", "ssgsea")) {
    score <- function(x) score_samples_l1000pai(
      x, TERM2GENE = toy$t2g, method = method, n_null = 100L,
      rng_seed = 19L, verbose = FALSE)
    cohort <- score(toy$expr)
    alone <- score(toy$expr[, "s2", drop = FALSE])
    reordered <- score(toy$expr[rev(seq_len(nrow(toy$expr))), c("s3", "s2", "s1")])
    for (field in c("score", "z", "p_activated", "p_inhibited")) {
      expect_equal(alone[[field]], cohort[[field]][, "s2", drop = FALSE])
      expect_equal(reordered[[field]][, colnames(cohort[[field]])], cohort[[field]])
    }
  }
})


test_that("zero-expression arms remain finite but constant samples are rejected", {
  expr <- matrix(c(rep(0, 10), 1:30), ncol = 1L,
                 dimnames = list(paste0("G", 1:40), "s1"))
  t2g <- data.frame(term = rep(c("Drug up", "Drug dn"), each = 10L),
                    gene = paste0("G", 1:20))
  for (method in c("rank_mean", "ssgsea")) {
    score <- function(x) score_samples_l1000pai(
      x, TERM2GENE = t2g, method = method, n_null = 50L, null_size_bin = 1L,
      min_size = 10L, rng_seed = 1L, verbose = FALSE)
    got <- score(expr)
    for (field in c("score", "z", "p_activated", "p_inhibited")) {
      expect_true(all(is.finite(got[[field]])))
    }
    for (value in c(0, 7, -4)) {
      constant <- expr
      constant[,] <- value
      expect_error(score(constant), "constant sample column")
      expect_error(score(cbind(expr, constant = constant[, 1])), "constant sample column")
    }
  }
})


test_that("empty and full-universe ssGSEA sets fail clearly", {
  expr <- matrix(1:10, ncol = 1L, dimnames = list(paste0("G", 1:10), "s1"))
  prep <- L1000PAI:::.l1000pai_sample_prep(expr, "ssgsea", 0.25)
  expect_error(L1000PAI:::.l1000pai_set_stat(prep, integer(), "ssgsea"),
               "must contain at least one gene")
  expect_error(L1000PAI:::.l1000pai_set_stat(prep, 1:10, "ssgsea"),
               "non-empty background")
  t2g <- rbind(data.frame(term = "Drug up", gene = rownames(expr)),
               data.frame(term = "Drug dn", gene = rownames(expr)[1:2]))
  expect_error(score_samples_l1000pai(
    expr, TERM2GENE = t2g, method = "ssgsea", reference = "none", n_null = 0L,
    min_size = 1L, verbose = FALSE), "non-empty background")
})


test_that("ssgsea direction flips both scores and the empirical tails", {
  toy <- make_toy()
  score <- function(direction) score_samples_l1000pai(
    toy$expr, TERM2GENE = toy$t2g, method = "ssgsea", n_null = 80L,
    score_direction = direction, rng_seed = 27L, verbose = FALSE)
  up <- score("up_minus_down")
  dn <- score("down_minus_up")
  expect_equal(up$score, -dn$score)
  expect_equal(up$z, -dn$z)
  expect_equal(up$p_activated, dn$p_inhibited)
  expect_equal(up$p_inhibited, dn$p_activated)
})


test_that("integer controls are not silently truncated or overflowed", {
  toy <- make_toy()
  defaults <- list(expr = toy$expr, TERM2GENE = toy$t2g, min_size = 1L,
                   verbose = FALSE)
  for (param in c("n_null", "null_size_bin", "min_size", "max_size")) {
    for (bad in list(1.5, NA_real_, Inf, NaN, .Machine$integer.max + 1)) {
      args <- defaults
      args[[param]] <- bad
      expect_error(do.call(score_samples_l1000pai, args), paste0("`", param, "`"))
    }
  }
})


test_that("invalid gene universes and overflowing ssGSEA weights fail explicitly", {
  toy <- make_toy()
  for (identifier in c("", "  ", NA_character_)) {
    expr <- toy$expr
    rownames(expr)[1] <- identifier
    expect_error(score_samples_l1000pai(expr, TERM2GENE = toy$t2g),
                 "gene identifiers must be non-missing and non-empty")
  }
  for (rows in list(integer(), 1L)) {
    expect_error(score_samples_l1000pai(toy$expr[rows, , drop = FALSE],
                                       TERM2GENE = toy$t2g),
                 "at least two genes")
  }
  expect_error(score_samples_l1000pai(
    toy$expr, TERM2GENE = toy$t2g, method = "ssgsea", tau = 1000,
    reference = "none", n_null = 0L, verbose = FALSE),
    "non-finite ssGSEA rank weights")
})
