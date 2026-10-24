test_that("the synthetic fixture gives calls in both directions", {
  fx <- toy_plot_fixture()
  cls <- table(result_table(fx$x)$activity_class)
  expect_gt(cls[["Activated"]], 10)
  expect_gt(cls[["Inhibited"]], 10)
})

test_that("target plots return a plot with values and a recommended size", {
  fx <- toy_plot_fixture(); x <- fx$x
  top <- result_table(x)$term[which.max(result_table(x)$activity_score)]
  plots <- list(
    pai_dotplot(x, n = 5, library = "crispr"),
    pai_armplot(x, n_label = 4, library = "crispr"),
    pai_gseaplot(x, top, term2gene = fx$t2g, library = "crispr"),
    pai_evidence_card(x, n = 8, library = "crispr"),
    pai_compareplot(list(first = x, second = x), n = 4, library = "crispr"),
    pai_cnetplot(x, n = 4, n_genes = 4, library = "crispr"))
  for (p in plots) {
    expect_s3_class(p, "gg")
    expect_true(!is.null(attr(p, "plot_data")))
    expect_length(attr(p, "size_mm"), 2)
    expect_true(draw_ok(p))
  }
  expect_error(pai_familyplot(x), "Chem")
})

test_that("the evidence card shares n rows between the directions", {
  sr <- L1000PAI:::.split_rows
  expect_equal(unname(sr(20, c(Activated = 1093, Inhibited = 743))), c(12, 8))
  expect_equal(unname(sr(7, c(Activated = 1093, Inhibited = 743))), c(4, 3))
  expect_equal(unname(sr(20, c(Activated = 3, Inhibited = 743))), c(1, 19))
  expect_equal(unname(sr(10, c(Activated = 3, Inhibited = 4))), c(3, 4))
  expect_equal(unname(sr(1, c(Activated = 3, Inhibited = 743))), c(0, 1))
  fx <- toy_plot_fixture(); x <- fx$x
  card <- attr(pai_evidence_card(x, n = 9, library = "crispr"), "plot_data")
  n_calls <- table(result_table(x)$activity_class)[c("Activated", "Inhibited")]
  want <- sr(9, c(Activated = n_calls[["Activated"]], Inhibited = n_calls[["Inhibited"]]))
  expect_equal(nrow(card), 9)
  expect_equal(sum(card$activity_class == "Activated"), unname(want[["Activated"]]))
  inh <- attr(pai_evidence_card(x, n = 5, direction = "Inhibited", library = "crispr"), "plot_data")
  expect_true(all(inh$activity_class == "Inhibited"))
  expect_equal(nrow(inh), 5)
})

test_that("pai_enrich() takes a GMT file and feeds every pathway plot", {
  fx <- toy_plot_fixture(); x <- fx$x
  expect_error(pai_enrich(x, collections = fx$gmt), "gene-level library")
  e <- pai_enrich(x, collections = c(Toy = fx$gmt), min_size = 5, library = "crispr")
  expect_s3_class(e, "pai_enrichment")
  expect_output(print(e), "pai_enrichment")
  expect_true(all(e$table$collection == "Toy"))
  tt <- pai_term_table(e)
  expect_s3_class(tt, "pai_term_table")
  expect_identical(pai_term_table(tt), tt)
  # balanced z of one term, from its counts and the overall share of Activated calls
  r <- tt$table[which.max(tt$table$count), ]
  p <- tt$p_activated; n <- r$n_activated + r$n_inhibited
  expect_equal(r$zscore_balanced, (r$n_activated - n * p) / sqrt(n * p * (1 - p)))
  expect_equal(r$zscore, (r$n_activated - r$n_inhibited) / sqrt(r$count))
  sig <- sum(tt$table$p.adjust < 0.05)
  expect_gte(sig, 3)
  plots <- suppressMessages(list(
    pai_enrichplot(e, n = 3), pai_gobubble(e), pai_gobubble(e, width = 92), pai_gobar(e, n = 4),
    pai_gobar(e, n = 4, z = "goplot", order = "padj"), pai_gocircle(e, n = 4), pai_gocircle(e, n = 4, members = "called", label_genes = 0),
    pai_gochord(e, n = 3, n_genes = 10), pai_goheat(e, n = 4, n_genes = 10), pai_goheat(e, n = 4, n_genes = 10, order = "cluster"),
    pai_govenn(e, n = 3), pai_gocluster(e, n = 3, n_genes = 15), pai_armmap(e, n = 3), pai_ladder(e, n = 4),
    pai_ladder(e, n = 4, order = "median"), pai_term_cnetplot(e, n = 3, n_genes = 5),
    pai_term_cnetplot(e, n = 3, n_genes = 5, layout = "network")))
  for (p in plots) {
    expect_true(!is.null(attr(p, "plot_data")))
    expect_true(draw_ok(p))
  }
  # GOCircle draws every scored member and names one extreme per direction by default
  d0 <- attr(pai_gocircle(e, n = 4), "plot_data")
  d1 <- attr(pai_gocircle(e, n = 4, members = "all", label_genes = 1), "plot_data")
  expect_equal(d0$genes, d1$genes)
  dc <- attr(pai_gocircle(e, n = 4, members = "called", label_genes = 0), "plot_data")
  expect_true(nrow(dc$genes) <= nrow(d0$genes))
  # z defaults to the balanced z-score
  b <- attr(pai_gobar(e, n = 4), "plot_data")
  expect_equal(b$zval, b$zscore_balanced)
  # terms by name, and a clear stop when a term is not there
  nm <- tt$table$Description[order(tt$table$p.adjust)][1:2]
  lad <- attr(pai_ladder(e, terms = nm), "plot_data")
  expect_equal(lad$terms$Description, nm)
  expect_error(pai_ladder(e, terms = "no such pathway"), "Not found")
  expect_error(pai_gobubble(x), "pai_enrich")
})

test_that("hub analysis and ranked pathway GSEA work on a GMT of one's own", {
  fx <- toy_plot_fixture(); x <- fx$x
  e <- pai_enrich(x, collections = c(Toy = fx$gmt), min_size = 5, library = "crispr")
  h <- pai_enrich_hubs(e, direction = "Inhibited", similarity = 0.2)
  expect_s3_class(h, "pai_hubs")
  expect_true(all(!is.na(h$genes$activity_score)))
  for (p in list(pai_hubplot(h, n_genes = 8, n_modules = 3), pai_chordplot(h, n_genes = 8, n_modules = 3),
                 pai_enrich_cnetplot(h, n_modules = 3, n_genes = 8), pai_emapplot(h, n = 10, n_modules = 3)))
    expect_true(draw_ok(p))
  g <- pai_gsea_pathways(x, collections = c(Toy = fx$gmt), min_size = 5, library = "crispr")
  expect_s3_class(g, "pai_pathway_gsea")
  expect_true(draw_ok(pai_pathwayplot(g, ora = e, n = 2)))
})

test_that("palette, style sheet, report and pai_save()", {
  fx <- toy_plot_fixture(); x <- fx$x
  old <- pai_use_palette("classic")
  expect_equal(old, "house")
  expect_equal(L1000PAI:::.pal$marks[["Discordant"]], "#6E5A8A")
  expect_equal(pai_use_palette(old), "classic")
  expect_equal(L1000PAI:::.pal$marks[["Discordant"]], "#A97FD0")
  expect_true(draw_ok(pai_style_guide()))
  e <- pai_enrich(x, collections = c(Toy = fx$gmt), min_size = 5, library = "crispr")
  rep <- pai_report(x, e, n_terms = 3)
  expect_true(draw_ok(rep))
  f <- tempfile(fileext = ".png")
  expect_identical(pai_save(pai_dotplot(x, n = 3, library = "crispr"), f, dpi = 72), f)
  expect_true(file.exists(f))
  expect_error(pai_save(1, f), "plotting function")
})
