# A small synthetic paired knockout library and pathway collection, so that the plotting functions
# can be tested without the L1000 libraries or a network connection. Built once per session.
toy_plot_fixture <- local({
  cache <- NULL
  function() {
    if (!is.null(cache)) return(cache)
    set.seed(11)
    genes <- sprintf("G%04d", seq_len(1500))
    stat <- sort(stats::rnorm(1500, sd = 2), decreasing = TRUE)
    ranks <- data.frame(gene = genes, stat = stat, stringsAsFactors = FALSE)
    top <- genes[1:300]; bottom <- genes[1201:1500]; mid <- genes[301:1200]
    targets <- genes[seq(5, 1500, by = 10)]                       # 150 knocked-out genes
    role <- rep(c("act", "inh", "none"), c(45, 35, 70))
    t2g <- do.call(rbind, lapply(seq_along(targets), function(i) {
      # down_minus_up: a target is Activated when its down arm sits at the top of the ranking and
      # its up arm at the bottom
      sets <- switch(role[i],
        act = list(up = sample(bottom, 40), dn = sample(top, 40)),
        inh = list(up = sample(top, 40), dn = sample(bottom, 40)),
        none = list(up = sample(mid, 40), dn = sample(mid, 40)))
      data.frame(term = rep(paste0(targets[i], c("_Up", "_Down")), each = 40), gene = c(sets$up, sets$dn),
        stringsAsFactors = FALSE)
    }))
    x <- suppressWarnings(run_l1000_pai(ranks, term2gene = t2g, score_direction = "down_minus_up", seed = TRUE, rng_seed = 1,
      collect_annotation_tables = FALSE))
    act <- targets[role == "act"]; inh <- targets[role == "inh"]; none <- targets[role == "none"]
    paths <- c(
      list(PATH_INH_A = c(inh[1:14], none[1:3]), PATH_INH_B = c(inh[10:24], none[4:6]), PATH_INH_C = c(inh[20:32], none[7:9]),
           PATH_ACT_A = c(act[1:15], none[10:12]), PATH_ACT_B = c(act[12:28], none[13:14])),
      stats::setNames(lapply(1:7, function(j) sample(targets, 15)), sprintf("PATH_MIXED_%d", 1:7)))
    gmt <- tempfile(fileext = ".gmt")
    writeLines(vapply(names(paths), function(p) paste(c(p, tolower(gsub("_", " ", p)), paths[[p]]), collapse = "\t"),
      character(1)), gmt)
    cache <<- list(x = x, t2g = t2g, gmt = gmt, ranks = ranks)
    cache
  }
})

# Draw a plot to a file, so that errors raised at drawing time are caught too. The cairo device is
# the one pai_save() uses for PDF: it draws the arrows, dagger and >= sign of the labels, which the
# base pdf() device replaces with dots (an error under R CMD check --as-cran).
draw_ok <- function(p) {
  f <- tempfile(fileext = ".pdf")
  if (isTRUE(capabilities("cairo"))) grDevices::cairo_pdf(f, width = 7, height = 7) else grDevices::pdf(NULL)
  on.exit({ grDevices::dev.off(); unlink(f) })
  print(p)
  invisible(TRUE)
}
