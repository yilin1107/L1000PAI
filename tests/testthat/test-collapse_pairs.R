test_that("collapse_pairs supports up_minus_down scoring", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(2.0, -1.5),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(
    gsea_table,
    combine_p_method = "max",
    score_direction = "up_minus_down"
  )

  expect_equal(res$activity_score, 2.0 - (-1.5))
  expect_equal(res$activity_score_legacy, res$activity_score)
  expect_equal(res$directional_state, "Activated")
  expect_equal(res$p_coherent, 0.02)
  expect_equal(res$padj_coherent, 0.04)
  expect_equal(res$pvalue_semantics, "observed_tail")
})

test_that("collapse_pairs supports up and dn suffixes", {
  gsea_table <- data.frame(
    term = c("TP53 up", "TP53 dn"),
    NES = c(-1.2, 1.8),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(
    gsea_table,
    combine_p_method = "max",
    score_direction = "down_minus_up"
  )

  expect_equal(res$term, "TP53")
  expect_equal(res$activity_score, 1.8 - (-1.2))
})

test_that("observed-tail directional IUT adjusts activation and inhibition over 2m hypotheses", {
  gsea_table <- data.frame(
    term = c(
      "Act_Up", "Act_Down",
      "Inh_Up", "Inh_Down",
      "Disc_Up", "Disc_Down"
    ),
    NES = c(-1.2, 1.8, 1.3, -1.7, 1.1, 1.4),
    pvalue = c(0.01, 0.02, 0.03, 0.04, 0.005, 0.006),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(gsea_table, score_direction = "down_minus_up")
  act <- res[res$term == "Act", ]
  inh <- res[res$term == "Inh", ]
  disc <- res[res$term == "Disc", ]

  expect_equal(act$p_activation, 0.02)
  expect_equal(act$p_inhibition, 1)
  expect_equal(act$p_coherent, 0.02)
  expect_equal(act$padj_activation, 0.12)
  expect_equal(act$padj_coherent, 0.12)
  expect_equal(act$directional_state, "Activated")

  expect_equal(inh$p_activation, 1)
  expect_equal(inh$p_inhibition, 0.04)
  expect_equal(inh$p_coherent, 0.04)
  expect_equal(inh$padj_inhibition, 0.12)
  expect_equal(inh$padj_coherent, 0.12)
  expect_equal(inh$directional_state, "Inhibited")

  expect_equal(disc$p_activation, 1)
  expect_equal(disc$p_inhibition, 1)
  expect_equal(disc$p_coherent, 1)
  expect_equal(disc$padj_coherent, 1)
  expect_equal(disc$directional_state, "Discordant")

  # The legacy two-arm p-values are adjusted over m pairs, whereas the two
  # directional conjunction families above are deliberately adjusted over 2m.
  expect_equal(act$p_both_arms, 0.02)
  expect_equal(act$padj_both_arms, 0.03)
  expect_equal(inh$p_both_arms, 0.04)
  expect_equal(inh$padj_both_arms, 0.04)
  expect_equal(disc$p_both_arms, 0.006)
  expect_equal(disc$padj_both_arms, 0.018)
  expect_true(all(res$pvalue_semantics == "observed_tail"))
})

test_that("two-sided p-values are converted to directional one-sided evidence", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(-1.2, 1.8),
    pvalue = c(0.02, 0.04),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(
    gsea_table,
    score_direction = "down_minus_up",
    pvalue_semantics = "two_sided"
  )

  expect_equal(res$p_activation, 0.02)
  expect_equal(res$p_inhibition, 0.99)
  expect_equal(res$p_coherent, 0.02)
  expect_equal(res$padj_activation, 0.04)
  expect_equal(res$padj_coherent, 0.04)
  expect_equal(res$directional_state, "Activated")
  expect_equal(res$pvalue_semantics, "two_sided")
})

test_that("Fisher combination has exact expected statistic and p-value", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(-1.2, 1.8),
    pvalue = c(0.01, 0.02),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(gsea_table, combine_p_method = "fisher")

  expect_equal(res$combined_stat, 17.03438638283247, tolerance = 1e-12)
  expect_equal(res$p_combined, 0.001903438638283249, tolerance = 1e-15)
})

test_that("a single extreme arm cannot pass the coherent conjunction gate", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(-4, 0.2),
    pvalue = c(1e-12, 0.8),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(gsea_table, combine_p_method = "fisher")

  expect_lt(res$p_combined, 1e-10)
  expect_equal(res$p_activation, 0.8)
  expect_equal(res$p_coherent, 0.8)
  expect_equal(res$padj_coherent, 1)
})

test_that("missing p-values propagate without changing the adjustment family", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(-1.2, 1.8),
    pvalue = c(NA_real_, 0.02),
    stringsAsFactors = FALSE
  )

  res <- collapse_pairs(gsea_table)

  expect_true(is.na(res$p_both_arms))
  expect_true(is.na(res$padj_both_arms))
  expect_true(is.na(res$p_activation))
  expect_true(is.na(res$padj_activation))
  expect_true(is.na(res$p_coherent))
  expect_true(is.na(res$padj_coherent))
})

test_that("collapse_pairs rejects ambiguous pairs and invalid inputs", {
  duplicate <- data.frame(
    term = c("DrugA_Up", "DrugA_Up", "DrugA_Down"),
    NES = c(-1, -1.1, 1),
    pvalue = c(0.01, 0.02, 0.03),
    stringsAsFactors = FALSE
  )
  expect_error(collapse_pairs(duplicate), "Duplicate directional signatures")

  invalid_p <- data.frame(
    term = c("DrugA_Up", "DrugA_Down"),
    NES = c(-1, 1),
    pvalue = c(-0.1, 0.02),
    stringsAsFactors = FALSE
  )
  expect_error(collapse_pairs(invalid_p), "between 0 and 1")
  invalid_p$pvalue[1] <- 1.1
  expect_error(collapse_pairs(invalid_p), "between 0 and 1")
  invalid_p$pvalue[1] <- Inf
  expect_error(collapse_pairs(invalid_p), "between 0 and 1")

  nonnumeric <- invalid_p
  nonnumeric$pvalue <- as.character(nonnumeric$pvalue)
  expect_error(collapse_pairs(nonnumeric), "pvalue.*numeric")
  nonnumeric <- invalid_p
  nonnumeric$NES <- as.character(nonnumeric$NES)
  expect_error(collapse_pairs(nonnumeric), "NES.*numeric")

  expect_error(collapse_pairs(invalid_p, combine_p_method = "bad"), "combine_p_method")
  expect_error(collapse_pairs(invalid_p, pair_p_adjust_method = "bad"), "pair_p_adjust_method")
})

test_that("incomplete pairs are reported and excluded", {
  gsea_table <- data.frame(
    term = c("DrugA_Up", "DrugA_Down", "DrugB_Up"),
    NES = c(-1.2, 1.8, -0.8),
    pvalue = c(0.01, 0.02, 0.2),
    stringsAsFactors = FALSE
  )

  expect_warning(
    res <- collapse_pairs(gsea_table),
    "Incomplete Up/Down signatures.*DrugB"
  )
  expect_equal(res$term, "DrugA")
})

test_that("strip_direction_suffix removes dn suffixes", {
  expect_equal(
    strip_direction_suffix(c("TP53 up", "TP53 dn", "RB1-dn")),
    c("TP53", "TP53", "RB1")
  )
})

test_that("mixed estimable and missing NES pairs retain unknown results safely", {
  input <- data.frame(
    term = paste0(rep(c("Valid", "MissingUp", "MissingDown", "MissingBoth", "Zero"),
                      each = 2), rep(c("_Up", "_Down"), 5)),
    NES = c(-2, 2, NA, 2, -2, NA, NA, NA, 0, 2),
    pvalue = rep(c(0.001, 0.002), 5)
  )
  for (direction in c("down_minus_up", "up_minus_down")) {
    for (semantics in c("observed_tail", "two_sided")) {
      for (combine in c("max", "fisher")) {
        got <- collapse_pairs(input, score_direction = direction,
                              pvalue_semantics = semantics,
                              combine_p_method = combine)
        expect_equal(nrow(got), 5L)
        good <- got[got$term == "Valid", ]
        expected_sign <- if (direction == "down_minus_up") 1 else -1
        expect_equal(good$activity_score, expected_sign * 4)
        expect_equal(good$directional_state,
                     if (expected_sign == 1) "Activated" else "Inhibited")
        missing <- got[grepl("^Missing", got$term), ]
        expect_true(all(is.na(missing$activity_score)))
        expect_true(all(is.na(missing$directional_state)))
        expect_true(all(is.na(missing$p_coherent)))
        expect_true(all(is.na(missing$padj_coherent)))
        expect_equal(got$directional_state[got$term == "Zero"], "Incomplete")
        classified <- classify_activity(got, score_direction = direction)
        expect_true(all(as.character(classified$activity_class[
          grepl("^Missing", classified$term)]) == "Not_significant"))
        expect_false(any(classified$is_significant[
          grepl("^Missing", classified$term)]))
      }
    }
  }
})

test_that("entirely missing NES and NaN inputs remain non-estimable", {
  input <- data.frame(term = c("A_Up", "A_Down", "B_Up", "B_Down"),
                      NES = c(NA_real_, NA_real_, NaN, 0),
                      pvalue = c(.01, .02, NA_real_, .02))
  got <- collapse_pairs(input)
  expect_equal(nrow(got), 2L)
  expect_true(all(is.na(got$directional_state)))
  expect_true(all(is.na(got$activity_score)))
  expect_true(all(is.na(got$p_coherent)))
  expect_true(all(is.na(got$padj_coherent)))
  classified <- classify_activity(got)
  expect_true(all(as.character(classified$activity_class) == "Not_significant"))
  expect_false(any(classified$is_significant))
})
