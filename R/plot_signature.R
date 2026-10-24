# Displays that rest on the two arms of every L1000PAI score.
#
# These plots rest on what only L1000PAI has: every target is scored on two arms (the genes the
# perturbation raises and the genes it lowers), and a call needs both arms to agree.
#   pai_evidence_card()  the result table of chosen or top-n targets with a two-arm glyph per row
#   pai_armmap()         where the top-n enriched pathways' members sit in the oriented arm plane
#   pai_ladder()         every member's activity score, pathway by pathway, against all targets
#   pai_term_cnetplot()  term-gene network with direction: Inhibited targets left, Activated right
#   pai_report()         one page combining the above
#   pai_save()           save at the recommended size (it grows with n)
# n works the same way everywhere: n = 10 or n = 20 draws the top 10 or 20 terms (rows, in total, for
# the evidence card); `terms` / `targets` choose them by name instead.

# ---- shared helpers ---------------------------------------------------------------------------
# Members of a term among the scored targets (same background as the ORA).
.term_members <- function(tt, keys) {
  stopifnot(inherits(tt, "pai_term_table"), !is.null(tt$sets))
  bg <- tt$background
  cols <- unique(sub("\\|.*$", "", keys)); out <- list()
  for (cl in cols) {
    gs <- .read_gmt2(tt$sets[[cl]]$file)$TERM2GENE
    for (k in keys[startsWith(keys, paste0(cl, "|"))]) {
      id <- sub("^[^|]*\\|", "", k)
      out[[k]] <- sort(intersect(unique(gs$gene[gs$term == id]), bg))
    }
  }
  out[keys]
}

# Which terms to draw: chosen by key ("Reactome|DNA Replication", "GO_BP|GO:0007049") or by name,
# or else the top n non-redundant significant terms. Fewer than n available -> a message, not an error.
.select_terms <- function(tt, terms = NULL, n = 10, cutoff = 0.05, max_overlap = 0.5, collections = NULL,
                          direction = c("both", "Activated", "Inhibited")) {
  direction <- match.arg(direction)
  if (!is.null(terms)) {
    tab <- tt$table
    # key first; then the exact name; then the name ignoring case. A name shared by several
    # databases (e.g. "DNA replication") takes the most significant one and says so.
    i <- vapply(terms, function(z) {
      hit <- which(tab$key == z)
      if (!length(hit)) hit <- which(tab$Description == z)
      if (!length(hit)) hit <- which(tolower(tab$Description) == tolower(z))
      if (!length(hit)) return(NA_integer_)
      if (length(hit) > 1) {
        hit <- hit[base::order(tab$p.adjust[hit], tab$key[hit])]
        message(sprintf("\"%s\" matches %d terms; using %s (give a key such as %s to choose another).",
          z, length(hit), tab$key[hit[1]], tab$key[hit[2]]))
      }
      hit[1]
    }, integer(1))
    if (anyNA(i)) stop("Not found among the tested terms: ", paste(terms[is.na(i)], collapse = ", "),
      ". Use a key such as \"Reactome|DNA Replication\" or the term name.", call. = FALSE)
    s <- tt$table[i, ]; s$tid <- sprintf("T%d", seq_len(nrow(s))); rownames(s) <- NULL
    return(s)
  }
  stopifnot(is.numeric(n), length(n) == 1, n >= 1)
  s <- .pick_terms(tt, n, cutoff, max_overlap, collections, direction)
  if (!nrow(s)) {
    pool <- tt$table
    if (!is.null(collections)) pool <- pool[pool$collection %in% collections, ]
    if (direction != "both") pool <- pool[pool$direction == direction, ]
    stop(sprintf(paste0("No term has adj. P < %s%s (smallest: %s). To draw the top terms descriptively, use a larger ",
      "cutoff (e.g. cutoff = 1) or name them with terms = ; they are then not significant."), cutoff,
      if (direction == "both") "" else sprintf(" in the %s direction", direction),
      if (nrow(pool)) signif(min(pool$p.adjust), 2) else "no term tested"), call. = FALSE)
  }
  if (nrow(s) < n) message(sprintf("Only %d non-redundant terms pass adj. P < %s; drawing %d.", nrow(s), cutoff, nrow(s)))
  s
}

.sel_note <- function(sel, n_req, terms, max_overlap, direction) {
  if (!is.null(terms)) return(sprintf("%d chosen terms", nrow(sel)))
  sprintf("Top %d %sterms%s", nrow(sel), if (is.finite(max_overlap)) "non-redundant " else "",
    if (identical(direction, "both")) "" else sprintf(" (%s direction)", tolower(direction)))
}

# Recommended size travels with the plot; pai_save() uses it.
.with_size <- function(p, width, height) { attr(p, "size_mm") <- c(width = width, height = height); p }

# Save a plot at its recommended size
# 
# Functions whose size depends on n (top 10, top 20, ...) attach a recommended width and height
# in millimetres; pai_save() uses it unless width or height is given. PDF uses cairo (fonts
# embedded); PNG is written at `dpi`.
#' Save a plot at its recommended size
#'
#' Every plotting function attaches a recommended width and height in millimetres, which grow
#' with `n`; `pai_save()` uses them unless `width` or `height` is given. PDF is written with
#' cairo (fonts embedded); other formats at `dpi`.
#'
#' Save PDF files with `pai_save()` or with `ggplot2::ggsave(..., device = grDevices::cairo_pdf)`.
#' The labels use arrows, a dagger and the greater-or-equal sign, which the base `pdf()` device
#' (the default of `ggsave()` for PDF files and of `Rplots.pdf`) cannot draw: it replaces them
#' with dots, so that, for example, the counts "0 down, 17 up" of [pai_gobar()] lose their
#' direction. Screen devices, PNG and SVG are not affected.
#'
#' @param p A plot made by an L1000PAI plotting function.
#' @param file File name; the extension sets the format.
#' @param width,height Size in millimetres; `NULL` uses the recommended size.
#' @param dpi Resolution of raster formats.
#' @return `file`, invisibly.
#' @family plot style
#' @export
pai_save <- function(p, file, width = NULL, height = NULL, dpi = 300) {
  if (!inherits(p, c("gg", "patchwork"))) stop("`p` must be a plot made by an L1000PAI plotting function.", call. = FALSE)
  sz <- attr(p, "size_mm") %||% c(width = 180, height = 120)
  w <- width %||% sz[["width"]]; h <- height %||% sz[["height"]]
  dev <- if (grepl("\\.pdf$", file, ignore.case = TRUE)) grDevices::cairo_pdf else NULL
  if (is.null(dev)) ggplot2::ggsave(file, p, width = w, height = h, units = "mm", dpi = dpi, bg = "white")
  else ggplot2::ggsave(file, p, width = w, height = h, units = "mm", device = dev)
  invisible(file)
}

# ---- 1. Arm map: pathways in the oriented arm plane ----------------------------------------------
# Every scored target of a term is drawn at its two oriented arm NES. The arm plane of a CRISPR
# library is strongly clustered, so the summary is the share of the term's targets in each
# quadrant. The first panel shows all scored targets, as the reference for every other panel.
.quadrant <- function(x, y) ifelse(x >= 0 & y >= 0, "act", ifelse(x < 0 & y < 0, "inh", ifelse(x < 0, "dis_tl", "dis_br")))
.auto_ncol <- function(k) if (k <= 3) k else if (k <= 9) 3 else if (k <= 16) 4 else if (k <= 25) 5 else 6

# Where each enriched pathway sits on the two arms
# 
# @param x L1000PAIResult of a gene-level library.
# @param tt pai_term_table(x, pai_enrich(x)).
# @param n Number of terms (top n non-redundant terms by adjusted P), e.g. 10 or 20.
# @param terms Instead of n: term keys ("Reactome|DNA Replication") or term names.
# @param direction "both", "Activated" or "Inhibited": keep terms whose better ORA direction is this one.
# @param max_overlap Jaccard limit of the non-redundant pick; Inf keeps every term.
# @param ncol Panels per row; NULL chooses from n.
# @param width Page width in mm the layout is planned for (text positions and the recommended height).
#' Where each enriched pathway sits on the two arms
#'
#' Every scored target of each top term on the oriented arm plane, with the share of the term's
#' targets in each quadrant. The framed first panel shows all scored targets as the reference.
#'
#' @inheritParams pai_plot_params
#' @param ncol Panels per row; `NULL` chooses from `n`.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_armmap <- function(e, n = 8, terms = NULL, direction = c("both", "Activated", "Inhibited"), cutoff = 0.05,
                       max_overlap = 0.5, collections = NULL, ncol = NULL, width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction)
  sd <- tt$score_direction; sem <- tt$semantics
  df <- tt$targets
  if (!all(c("NES_Up", "NES_Down") %in% names(df))) stop("pai_armmap() needs a paired library (Up and Down arms).", call. = FALSE)
  df <- df[!is.na(df$NES_Up) & !is.na(df$NES_Down), ]
  o <- .pai_orient(df, sd); df$x_oriented <- o$x; df$y_oriented <- o$y
  df$cls <- df$activity_class
  df$quadrant <- .quadrant(df$x_oriented, df$y_oriented)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  mem <- .term_members(tt, sel$key)
  k <- nrow(sel); ids <- c("T0", sel$tid)
  pts <- do.call(rbind, c(list(transform(df[, c("gene", "x_oriented", "y_oriented", "activity_score", "cls", "quadrant")], tid = "T0")),
    lapply(seq_len(k), function(i) {
      d <- df[df$gene %in% mem[[i]], c("gene", "x_oriented", "y_oriented", "activity_score", "cls", "quadrant")]
      d$tid <- sel$tid[i]; d })))
  qs <- c("act", "inh", "dis_tl", "dis_br")
  bg_share <- as.numeric(table(factor(df$quadrant, qs))) / nrow(df)
  lim <- ceiling(max(abs(c(df$x_oriented, df$y_oriented))) * 1.02 * 2) / 2
  # layout planned for `width` mm: panel size, text sizes and the offsets of the two title lines
  nc <- ncol %||% .auto_ncol(k + 1); nr <- ceiling((k + 1) / nc)
  fw <- (width - 18 - (nc - 1) * 4.2) / nc                       # panel width in mm
  f <- if (fw >= 52) 1 else max(0.8, fw / 52 + 0.15)             # text scale for small panels
  shares <- do.call(rbind, lapply(ids, function(t) {
    d <- pts[pts$tid == t, ]; kk <- as.numeric(table(factor(d$quadrant, qs)))
    data.frame(tid = t, quadrant = qs, k = kk, n = nrow(d), share = kk / nrow(d), all_share = bg_share,
      x = c(0.95, -0.95, -0.95, 0.95) * lim, y = c(0.95, -0.95, 0.95, -0.95) * lim, h = c(1, 0, 0, 1), v = c(1, 0, 1, 0),
      cls = c("Activated", "Inhibited", "Discordant", "Discordant"))
  }))
  pct <- function(v) ifelse(v > 0 & v < 0.005, "<1%", sprintf("%d%%", as.integer(round(100 * v))))
  shares$lab <- pct(shares$share)
  info <- data.frame(tid = ids, n = vapply(ids, function(t) sum(pts$tid == t), numeric(1)))
  chars <- max(16, floor(fw / (if (fw >= 45) 1.42 else 1.22)) - 3)
  info$title <- c("All scored targets", sprintf("%s  %s", sel$tid, .short_label(sel$Description, chars)))
  info$line2 <- c(sprintf('"%s targets  \u00b7  reference"', .fmt_n(nrow(df))),
    .fmt_p10_expr(sel$p.adjust, prefix = sprintf(if (fw >= 45) "%s  \u00b7  %s targets  \u00b7  q =" else "%s \u00b7 n = %s \u00b7 q =",
      .db_tag(sel$collection), .fmt_n(info$n[-1]))))
  info$ref <- info$tid == "T0"
  bgd <- df[, c("x_oriented", "y_oriented")]
  pts <- pts[order(factor(pts$cls, c("Not_significant", "Discordant", "Inhibited", "Activated"))), ]
  bw_ <- 2
  brk <- function(l) { b <- seq(-floor(lim / bw_) * bw_, floor(lim / bw_) * bw_, by = bw_); b[abs(b) < lim - 1e-9] }  # no labels at the edges
  tsize <- if (fw >= 45) 7.2 else 6.4                             # panel title in pt; subtitle is 0.78 of it
  # Each panel is its own small plot with a header in millimetres, so titles never collide
  # with the panel whatever size the page is drawn at.
  cell <- function(j) {
    id <- ids[j]; d <- pts[pts$tid == id, ]; sh <- shares[shares$tid == id, ]
    left <- (j - 1) %% nc == 0; bottom <- j + nc > length(ids)
    g <- ggplot2::ggplot(d, ggplot2::aes(x_oriented, y_oriented)) +
      .arm_canvas(lim, sem, tt$threshold, iso = FALSE, labels = FALSE) +
      (if (id != "T0") ggplot2::geom_point(data = bgd, colour = "#ECEAE6", size = 0.22 * f, stroke = 0)) +
      ggplot2::geom_point(ggplot2::aes(colour = cls, shape = cls), size = (if (id == "T0") 0.35 else 0.8) * f, stroke = 0, alpha = 0.95) +
      ggplot2::geom_text(data = sh, inherit.aes = FALSE, ggplot2::aes(x = x, y = y, label = lab, hjust = h, vjust = v),
        colour = unname(.pal$text[as.character(sh$cls)]), size = 2.3 * f, fontface = "bold", family = .pai_font()) +
      (if (id == "T0") ggplot2::annotate("rect", xmin = -lim, xmax = lim, ymin = -lim, ymax = lim, fill = NA,
        colour = .tok$ink2, linewidth = 0.45)) +
      ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
      ggplot2::scale_shape_manual(values = c(Not_significant = 16, Discordant = 18, Inhibited = 16, Activated = 16), guide = "none") +
      ggplot2::scale_x_continuous(breaks = brk, labels = .lab_minus) + ggplot2::scale_y_continuous(breaks = brk, labels = .lab_minus) +
      ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
      ggplot2::labs(x = NULL, y = NULL) + theme_pai(axes = "xy") +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 1.5, 0.5, 0, "mm"))
    if (!bottom) g <- g + ggplot2::theme(axis.text.x = ggplot2::element_blank())
    if (!left) g <- g + ggplot2::theme(axis.text.y = ggplot2::element_blank())
    .stack(list(.pai_header(info$title[j], info$line2[j], size = tsize, subtitle_parse = TRUE), g), c(7.4, NA))
  }
  grid_p <- patchwork::wrap_elements(full = patchwork::wrap_plots(lapply(seq_along(ids), cell), ncol = nc)) +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  ylab <- .wrap_grob(grid::textGrob(o$ylab, rot = 90, gp = .gp(6.5, .tok$ink)))
  xlab <- .wrap_grob(grid::textGrob(o$xlab, gp = .gp(6.5, .tok$ink)))
  p <- patchwork::wrap_plots(list(ylab, grid_p, patchwork::plot_spacer(), xlab), ncol = 2,
    widths = grid::unit(c(4, 1), c("mm", "null")), heights = grid::unit(c(1, 5), c("null", "mm")))
  chips <- .pai_chips(c(stats::setNames(.pal$marks[c("Activated", "Inhibited", "Discordant", "Not_significant")],
      sprintf("%s target", sem$short[c("Activated", "Inhibited", "Discordant", "Not_significant")])),
      `targets outside the term` = "#E4E2DD"), shapes = c("circle", "circle", "diamond", "circle", "circle"), per_row = 5, size = 5.9)
  key <- .pai_badge("Percent", "share of the panel\u2019s targets in each quadrant; compare each term with the framed first panel")
  out <- .stack(list(.pai_header(title %||% "Where each pathway sits on the two arms",
      sprintf("%s; every scored target in the term is drawn on the oriented arm plane", .sel_note(sel, n, terms, max_overlap, direction))),
    chips, p, key), c(.header_mm(1), 4.6, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 0.5, 6.2, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- list(terms = sel, members = pts[pts$tid != "T0", ],
    shares = shares[shares$tid != "T0", c("tid", "quadrant", "k", "n", "share", "all_share")],
    reference = shares[shares$tid == "T0", c("quadrant", "k", "n", "share")])
  .with_size(out, width, round(30 + nr * (fw + 6.2), 1))
}

# ---- 2. Activity ladder: every member's score, term by term --------------------------------------
# Activity ladder of enriched pathways
# 
# One row per term; every scored target of the term at its activity score, the members' median
# (black bar), and all scored targets as reference (grey band = middle 50%, grey line = median).
# @param n Number of terms, e.g. 10 or 20; or give `terms` (keys or names).
# @param order "padj" (selection order, most significant on top) or "median" (highest member median on top).
#' Activity ladder of the enriched pathways
#'
#' One row per term: every scored target of the term at its activity score, the members' median
#' (black bar) and all scored targets as reference (grey band = middle 50%, grey line = median).
#'
#' @inheritParams pai_plot_params
#' @param order `"padj"` (most significant on top) or `"median"` (highest member median on top).
#' @param row_mm Height of one row in millimetres.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_ladder <- function(e, n = 10, terms = NULL, direction = c("both", "Activated", "Inhibited"), order = c("padj", "median"),
                       cutoff = 0.05, max_overlap = 0.5, collections = NULL, seed = 1, row_mm = 5.6, width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction); order <- match.arg(order)
  sem <- tt$semantics
  df <- tt$targets; df <- df[!is.na(df$activity_score), ]
  df$cls <- df$activity_class
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  mem <- .term_members(tt, sel$key)
  k <- nrow(sel)
  med <- vapply(mem, function(g) stats::median(df$activity_score[df$gene %in% g]), numeric(1))
  rank_top <- if (order == "padj") seq_len(k) else base::order(-med, sel$p.adjust)   # position from the top
  rowpos <- integer(k); rowpos[rank_top] <- rev(seq_len(k))                          # y = k for the top row
  pts <- do.call(rbind, lapply(seq_len(k), function(i) {
    d <- df[df$gene %in% mem[[i]], c("gene", "activity_score", "cls")]; d$tid <- sel$tid[i]; d$row <- rowpos[i]; d }))
  set.seed(seed); pts$y <- pts$row + stats::runif(nrow(pts), -0.27, 0.27)
  sumr <- do.call(rbind, lapply(split(pts, pts$tid), function(d) data.frame(tid = d$tid[1], row = d$row[1], n = nrow(d),
    median = stats::median(d$activity_score), n_act = sum(d$cls == "Activated"), n_inh = sum(d$cls == "Inhibited"))))
  sumr <- sumr[match(sel$tid, sumr$tid), ]
  gq <- stats::quantile(df$activity_score, c(0.25, 0.5, 0.75), names = FALSE, type = 7)
  thr <- tt$threshold
  lim <- max(abs(pts$activity_score)) * 1.04
  pts <- pts[base::order(factor(pts$cls, c("Not_significant", "Discordant", "Inhibited", "Activated"))), ]
  called <- pts$cls %in% c("Activated", "Inhibited")
  ylim <- c(0.45, k + 0.55)
  ys <- ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0))
  lad <- ggplot2::ggplot(pts, ggplot2::aes(activity_score, y)) +
    ggplot2::annotate("rect", xmin = gq[1], xmax = gq[3], ymin = -Inf, ymax = Inf, fill = .tok$pill) +
    ggplot2::geom_vline(xintercept = gq[2], colour = .tok$muted, linewidth = 0.25) +
    ggplot2::geom_vline(xintercept = c(-thr, thr), colour = .tok$muted, linewidth = 0.22, linetype = "22") +
    ggplot2::geom_hline(yintercept = seq(1.5, k - 0.5, 1), colour = .tok$rule, linewidth = 0.2) +
    ggplot2::geom_point(data = pts[!called, ], ggplot2::aes(colour = cls, shape = cls), size = 0.55, stroke = 0, alpha = 0.9) +
    ggplot2::geom_point(data = pts[called, ], ggplot2::aes(fill = cls), shape = 21, colour = "white", size = 1.2, stroke = 0.25) +
    ggplot2::geom_segment(data = sumr, inherit.aes = FALSE, ggplot2::aes(x = median, xend = median, y = row - 0.36, yend = row + 0.36),
      colour = .tok$ink, linewidth = 0.7, lineend = "butt") +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_shape_manual(values = c(Not_significant = 16, Discordant = 18), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), breaks = scales::breaks_width(2), labels = .lab_minus) + ys +
    ggplot2::labs(x = sem$score_label, y = NULL) +
    theme_pai(axes = "x") + ggplot2::theme(axis.text.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(0, 2, 0, 2))
  key <- sel; key$row <- rowpos; key$db <- .db_tag(key$collection)
  key$lab <- sprintf("%s  %s", key$tid, .short_label(key$Description, max(16, floor((width * 1.45 / 4.45 - 18) / 1.15))))
  left <- ggplot2::ggplot(key, ggplot2::aes(y = row)) +
    geom_pai_pill(ggplot2::aes(x = 0, y = row, label = db), data = key, hjust = 0, size = 1.8) +
    ggplot2::geom_text(ggplot2::aes(x = 0.2, label = lab), hjust = 0, size = 2.1, colour = .tok$ink, family = .pai_font()) +
    ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + ys + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  sumr$lab_q <- .fmt_p10_expr(sel$p.adjust); sumr$dir <- sel$direction
  bw <- max(sumr$n_act, sumr$n_inh, 1)
  right <- ggplot2::ggplot(sumr, ggplot2::aes(y = row)) +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$axis, linewidth = 0.25) +
    ggplot2::geom_rect(ggplot2::aes(xmin = -n_inh / bw, xmax = 0, ymin = row - 0.22, ymax = row + 0.22), fill = .pal$marks[["Inhibited"]]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = 0, xmax = n_act / bw, ymin = row - 0.22, ymax = row + 0.22), fill = .pal$marks[["Activated"]]) +
    ggplot2::geom_text(ggplot2::aes(x = -1.08, label = n_inh), hjust = 1, size = 2, colour = .pal$text[["Inhibited"]], family = .pai_font(), fontface = "bold") +
    ggplot2::geom_text(ggplot2::aes(x = 1.08, label = n_act), hjust = 0, size = 2, colour = .pal$text[["Activated"]], family = .pai_font(), fontface = "bold") +
    ggplot2::geom_text(ggplot2::aes(x = 1.75, label = lab_q, colour = dir), parse = TRUE, hjust = 0, size = 2, family = .pai_font()) +
    ggplot2::scale_colour_manual(values = .pal$text, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-1.6, 3.4), expand = c(0, 0)) + ys + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  colhead <- function(txt, hjust = 0, x = 0) .wrap_grob(grid::textGrob(txt, x = grid::unit(x, "npc"), hjust = hjust,
    gp = .gp(5.8, .tok$ink2, "bold")))
  rh <- ggplot2::ggplot() +
    ggplot2::annotate("text", x = c(-0.06, 0.06), y = 0, label = c(sem$short[["Inhibited"]], sem$short[["Activated"]]), hjust = c(1, 0),
      colour = c(.pal$text[["Inhibited"]], .pal$text[["Activated"]]), size = 2.05, fontface = "bold", family = .pai_font()) +
    ggplot2::annotate("text", x = 1.75, y = 0, label = "adj. P", hjust = 0, colour = .tok$ink2, size = 2.05, fontface = "bold", family = .pai_font()) +
    ggplot2::scale_x_continuous(limits = c(-1.6, 3.4), expand = c(0, 0)) + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  heads <- (colhead(if (order == "median") "Term  (highest member median on top)" else if (!is.null(terms)) "Term  (in the order given)"
      else "Term  (most significant on top)") |
    colhead("Every scored target in the term, at its activity score", 0.5, 0.5) | rh) + patchwork::plot_layout(widths = c(1.45, 2, 1))
  body <- (left | lad | right) + patchwork::plot_layout(widths = c(1.45, 2, 1))
  chips <- .pai_chips(c(stats::setNames(.pal$marks[c("Activated", "Inhibited", "Discordant", "Not_significant")],
      sem$short[c("Activated", "Inhibited", "Discordant", "Not_significant")]), `member median` = .tok$ink, `all targets, middle 50%` = "#E2E0DA"),
    shapes = c("circle", "circle", "diamond", "circle", "line", "square"), per_row = 6, size = 5.9)
  out <- .stack(list(.pai_header(title %||% "Activity ladder of the enriched pathways",
      sprintf("%s; the shaded band and grey line show all %s scored targets", .sel_note(sel, n, terms, max_overlap, direction), .fmt_n(nrow(df)))),
    chips, heads, body, .pai_badge("Read", sprintf("black bar outside the grey band: most of the term's targets lean to that side  \u00b7  dashed: |score| = %s", thr))),
    c(.header_mm(1), 4.6, 4, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 0.5, 4.6, 4, 1, 6.5), c("mm", "mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- list(terms = sel, members = pts, summary = sumr, global = gq)
  .with_size(out, width, round(.header_mm(1) + 0.5 + 4.6 + 4 + 6.5 + 11 + k * row_mm, 1))
}

# ---- 3. Evidence card: the result table with a two-arm glyph -----------------------------------
# Evidence card
# 
# The result-table columns with a two-arm glyph per row (upper bar = up arm,
# lower bar = down arm, both oriented so that right = toward activation).
# @param targets Targets to show, in this order; or leave NULL and give n.
# @param n Top n calls per direction by |activity score| (ties: adjusted P, then name), the rule of
#   pai_dotplot(); n = 10 with direction = "both" gives 10 Activated and 10 Inhibited rows.
# @param direction "both", "Activated" or "Inhibited".
#' Evidence card: the result table with a two-arm glyph
#'
#' One row per target with both NES, the activity score, the directional P values, the coherent
#' adjusted P, the directional state and the class; the glyph draws the two oriented arms (upper
#' bar = up arm, lower bar = down arm, right = toward activation).
#'
#' @inheritParams pai_plot_params
#' @param targets Targets to show, in this order. `NULL` takes the top `n` calls.
#' @param n Number of rows in total. With `direction = "both"` the rows are shared between the two
#'   directions in proportion to their numbers of calls (each direction with calls gets at least
#'   one row); within a direction targets are ranked by |activity score| (ties: coherent adjusted
#'   P, then name).
#' @param direction `"both"`, `"Activated"` or `"Inhibited"`.
#' @param row_mm Height of one row in millimetres.
#' @param subtitle Subtitle; `NULL` for the default.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_evidence_card <- function(x, targets = NULL, n = 20, direction = c("both", "Activated", "Inhibited"), library = NULL,
                              row_mm = 4.3, width = 130, title = NULL, subtitle = NULL) {
  direction <- match.arg(direction)
  lib <- .pai_library(x, library); sd <- .pai_score_direction(x, lib); sem <- .pai_semantics(lib, sd)
  df <- .pai_df(x); gene <- df$perturbagen %||% df$term
  auto <- is.null(targets)
  if (auto) {
    stopifnot(is.numeric(n), length(n) == 1, n >= 1)
    cls_all <- as.character(df$activity_class)
    dirs <- if (direction == "both") c("Activated", "Inhibited") else direction
    quota <- .split_rows(n, vapply(dirs, function(cl) sum(cls_all == cl), numeric(1)))
    idx <- unlist(lapply(dirs, function(cl) {
      w <- which(cls_all == cl)
      utils::head(w[base::order(-abs(df$activity_score[w]), df$padj_coherent[w], df$term[w])], quota[[cl]])
    }))
    if (!length(idx)) stop("No ", paste(dirs, collapse = "/"), " call in this result.", call. = FALSE)
    targets <- gene[idx]
  } else {
    idx <- match(targets, gene); idx[is.na(idx)] <- match(targets[is.na(idx)], df$term)
    if (anyNA(idx)) stop("Not in the result: ", paste(targets[is.na(idx)], collapse = ", "), call. = FALSE)
  }
  d <- df[idx, ]; d$label <- targets
  o <- .pai_orient(d, sd); d$x_oriented <- o$x; d$y_oriented <- o$y
  k <- nrow(d); d$row <- rev(seq_len(k)); d$cls <- as.character(d$activity_class)
  d$state <- as.character(d$directional_state)
  d$state_lab <- ifelse(d$state %in% names(sem$short), unname(sem$short[d$state]), d$state)
  drug_like <- identical(sem$short[["Activated"]], "Concordant")
  tw <- if (drug_like) c("concordance", "opposition") else c("activation", "inhibition")
  amax <- max(abs(c(df$NES_Up, df$NES_Down)), na.rm = TRUE)
  X <- c(target = 1, arms = 21, nes_up = 39.5, nes_dn = 51, score = 60, p_act = 73.5, p_inh = 85, q = 96.5, state = 99, class = 111)
  # names longer than about 8 characters (compound names) widen the name column; the Concordant and
  # Opposing pills need a little more room than Activated / Inhibited. Gene cards keep their size.
  shift <- max(0, ceiling(max(nchar(d$label)) * 1.3 + 2 - 13.5))
  X[-1] <- X[-1] + shift
  xr <- 122 + shift + if (drug_like) 3 else 0
  num <- data.frame(row = rep(d$row, 3), x = rep(X[c("nes_up", "nes_dn", "score")], each = k),
    lab = c(.fmt_num(d$NES_Up), .fmt_num(d$NES_Down), .fmt_num(d$activity_score)))
  pnum <- data.frame(row = rep(d$row, 3), x = rep(X[c("p_act", "p_inh", "q")], each = k),
    lab = c(.fmt_p10_expr(d$p_activation), .fmt_p10_expr(d$p_inhibition), .fmt_p10_expr(d$padj_coherent)))
  hd <- data.frame(x = X[c("target", "arms", "nes_up", "nes_dn", "score", "p_act", "p_inh", "q", "state", "class")],
    lab = c("target", "arms (oriented)", "NES Up", "NES Down", "score", "p act", "p inh", "q coh", "state", "class"),
    h = c(0, 0.5, 1, 1, 1, 1, 1, 1, 0, 0))
  aw <- 7.5
  arms <- rbind(data.frame(row = d$row, arm = "up", v = d$x_oriented, yo = 0.14),
                data.frame(row = d$row, arm = "dn", v = d$y_oriented, yo = -0.14))
  arms$x0 <- X[["arms"]]; arms$x1 <- X[["arms"]] + aw * arms$v / amax
  arms$side <- ifelse(arms$v >= 0, "Activated", "Inhibited")
  pill_df <- data.frame(row = d$row, x = X[["class"]], lab = sem$short[d$cls], cls = d$cls)
  # a rule between the class blocks when the card lists top calls of both directions
  brk <- which(utils::head(d$cls, -1) != utils::tail(d$cls, -1))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = d[d$row %% 2 == 0, ], ggplot2::aes(xmin = -0.5, xmax = xr, ymin = row - 0.5, ymax = row + 0.5), fill = .tok$pill) +
    ggplot2::geom_segment(ggplot2::aes(x = -0.5, xend = xr, y = k + 0.62, yend = k + 0.62), colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::geom_text(data = hd, ggplot2::aes(x = x, y = k + 1.05, label = lab, hjust = h), size = 1.95, colour = .tok$ink2,
      family = .pai_font(), fontface = "bold") +
    ggplot2::geom_text(data = d, ggplot2::aes(x = X[["target"]], y = row, label = label), hjust = 0, size = 2.2, fontface = "bold",
      colour = .tok$ink, family = .pai_font()) +
    ggplot2::geom_segment(data = d, ggplot2::aes(x = X[["arms"]] - aw, xend = X[["arms"]] + aw, y = row, yend = row), colour = .tok$rule, linewidth = 0.25) +
    ggplot2::geom_segment(data = d, ggplot2::aes(x = X[["arms"]], xend = X[["arms"]], y = row - 0.32, yend = row + 0.32), colour = .tok$ink2, linewidth = 0.25) +
    ggplot2::geom_segment(data = arms, ggplot2::aes(x = x0, xend = x1, y = row + yo, yend = row + yo, colour = side), linewidth = 1.15, lineend = "butt") +
    ggplot2::geom_text(data = num, ggplot2::aes(x = x, y = row, label = lab), hjust = 1, size = 2.05, colour = .tok$ink, family = .pai_font()) +
    ggplot2::geom_text(data = pnum, ggplot2::aes(x = x, y = row, label = lab), parse = TRUE, hjust = 1, size = 2.05, colour = .tok$ink, family = .pai_font()) +
    ggplot2::geom_text(data = d, ggplot2::aes(x = X[["state"]], y = row, label = state_lab, colour = state), hjust = 0, size = 2.0, family = .pai_font()) +
    ggplot2::geom_label(data = pill_df, ggplot2::aes(x = x, y = row, label = lab, fill = cls, colour = cls), hjust = 0, size = 1.95,
      family = .pai_font(), fontface = "bold", linewidth = 0, label.r = grid::unit(0.6, "mm"), label.padding = grid::unit(0.45, "mm")) +
    ggplot2::scale_colour_manual(values = .pal$text, guide = "none") +
    ggplot2::scale_fill_manual(values = .pal$tint, guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-0.5, xr), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.45, k + 1.4), expand = c(0, 0)) +
    ggplot2::theme_void() + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  if (length(brk)) p <- p + ggplot2::annotate("segment", x = -0.5, xend = xr, y = d$row[brk] - 0.5, yend = d$row[brk] - 0.5,
    colour = .tok$ink2, linewidth = 0.3)
  pills <- .pai_pills(list(list("L1000PAIResult", "object"), list("result_table()", "fn")),
    note = sprintf("L1000PAI %s  \u00b7  %s", as.character(utils::packageVersion("L1000PAI")), .pai_libname(lib)))
  arm_key <- .pai_chips(stats::setNames(c(.pal$marks[["Activated"]], .pal$marks[["Inhibited"]]), paste("arm toward", tw)),
    line = TRUE, size = 5.9)
  default_title <- if (auto) {
    if (direction == "both") sprintf("Top %d calls: %d %s, %d %s", nrow(d), sum(d$cls == "Activated"), sem$short[["Activated"]],
      sum(d$cls == "Inhibited"), sem$short[["Inhibited"]]) else sprintf("Top %d %s calls", nrow(d), sem$short[[direction]])
  } else "Evidence for selected targets"
  out <- .stack(list(.pai_header(title %||% default_title,
      subtitle %||% paste0("Arms glyph: upper bar = up arm, lower bar = down arm, both oriented so that right = toward ", tw[1])),
    pills, arm_key, p), c(.header_mm(1), 5.2, 4.4, NA)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 0.5, 5.2, 4.4, 1), c("mm", "mm", "mm", "null")))
  attr(out, "plot_data") <- d[, c("label", "NES_Up", "NES_Down", "x_oriented", "y_oriented", "activity_score", "p_activation",
    "p_inhibition", "padj_coherent", "directional_state", "activity_class")]
  .with_size(out, round(width * (xr + 0.5) / 122.5, 1), round(.header_mm(1) + 0.5 + 5.2 + 4.4 + (k + 1.2) * row_mm + 1, 1))
}

# n rows shared between the directions in proportion to their numbers of calls (largest remainder;
# every direction that has calls gets at least one row when n allows); a direction with fewer calls
# than its share passes the rest on.
.split_rows <- function(n, avail) {
  avail[] <- pmax(0, avail); q <- stats::setNames(numeric(length(avail)), names(avail))
  if (!sum(avail)) return(q)
  n <- min(n, sum(avail))
  raw <- n * avail / sum(avail); q[] <- floor(raw)
  need <- avail > 0 & q == 0
  if (any(need) && n >= sum(avail > 0)) q[need] <- 1
  while (sum(q) > n) { i <- which.max(q - raw); q[i] <- q[i] - 1 }
  while (sum(q) < n) {
    room <- avail - q; r <- raw - q; r[room <= 0] <- -Inf
    i <- which.max(r); q[i] <- q[i] + 1
  }
  pmin(q, avail)
}

# ---- 4. Direction-aware term-gene network (cnet) ------------------------------------------------
# Spread positions in [0, 1] near their barycentres, at least `gap` apart, in the given order.
.spread <- function(b, gap) {
  m <- length(b); if (!m) return(numeric()); y <- numeric(m)
  y[1] <- min(1, b[1]); if (m > 1) for (i in 2:m) y[i] <- min(b[i], y[i - 1] - gap)
  if (y[m] < 0) { y[m] <- 0; if (m > 1) for (i in (m - 1):1) y[i] <- max(y[i], y[i + 1] + gap) }
  pmin(y, 1)
}

# Term-gene network with direction
# 
# Terms: the top n non-redundant enriched terms (or `terms`). Genes: called targets of these
# terms, `n_genes` per direction, the most shared first (ties: |activity score|, then name).
# layout = "signed": Inhibited targets on the left and Activated targets on the right, each at its
# activity score; the terms sit in the middle column, each drawn as a bar split by the share of its
# called targets that are Inhibited (cool, left part) and Activated (warm, right part).
# layout = "network": the same nodes and colours in a force-directed layout.
#' Enriched terms and their called targets, by direction
#'
#' `layout = "signed"` puts the Inhibited targets on the left and the Activated targets on the
#' right, each at its activity score, with the terms in the middle as bars split by the share of
#' their called targets in each direction; `layout = "network"` draws the same nodes and colours
#' in a force-directed layout.
#'
#' @inheritParams pai_plot_params
#' @param n_genes Called targets per direction, those in most terms first.
#' @param layout `"signed"` or `"network"`.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_term_cnetplot <- function(e, n = 5, n_genes = 15, terms = NULL, direction = c("both", "Activated", "Inhibited"),
                              layout = c("signed", "network"), z = c("balanced", "goplot"), cutoff = 0.05, max_overlap = 0.5,
                              collections = NULL, seed = 1, width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction); layout <- match.arg(layout); z <- match.arg(z); zc <- .zcol(z)
  sem <- tt$semantics
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  k <- nrow(sel)
  memb <- do.call(rbind, lapply(seq_len(k), function(i) data.frame(tid = sel$tid[i],
    gene = strsplit(sel$geneID[i], "/", fixed = TRUE)[[1]], stringsAsFactors = FALSE)))
  memb$cls <- unname(tt$classes[memb$gene]); memb <- memb[memb$cls %in% c("Activated", "Inhibited"), ]
  deg <- table(memb$gene)
  cand <- data.frame(gene = names(deg), degree = as.integer(deg), stringsAsFactors = FALSE)
  cand$score <- unname(tt$scores[cand$gene]); cand$cls <- unname(tt$classes[cand$gene])
  side <- function(cl) { d <- cand[cand$cls == cl, ]; utils::head(d[base::order(-d$degree, -abs(d$score), d$gene), ], n_genes) }
  genes <- rbind(side("Inhibited"), side("Activated")); rownames(genes) <- NULL
  e <- memb[memb$gene %in% genes$gene, c("tid", "gene")]; rownames(e) <- NULL
  tn <- data.frame(tid = sel$tid, key = sel$key, label = sprintf("%s  %s", sel$tid, .short_label(sel$Description, 34)),
    db = .db_tag(sel$collection), n_act = sel$n_activated, n_inh = sel$n_inhibited, z = sel[[zc$col]], p = sel$p.adjust,
    stringsAsFactors = FALSE)
  tn$f_inh <- tn$n_inh / (tn$n_act + tn$n_inh)
  tn$line2 <- .fmt_p10_expr(tn$p, prefix = sprintf("%s  \u00b7  %d\u2193 %d\u2191  \u00b7  z %s  \u00b7  q =", tn$db, tn$n_inh, tn$n_act,
    .fmt_num(tn$z, 1)))
  lim <- max(abs(genes$score))
  fill_scale <- ggplot2::scale_fill_gradient2(low = .pal$marks[["Inhibited"]], mid = .pai_mid, high = .pal$marks[["Activated"]],
    limits = c(-lim, lim), labels = .lab_minus, name = "Activity score",
    guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(24, "mm"), barheight = ggplot2::unit(2, "mm"), title.position = "left"))
  size_scale <- ggplot2::scale_size_continuous(range = c(1.6, 3.6), limits = c(1, max(2, k)), breaks = unique(pmax(1, round(seq(1, max(2, k), length.out = 3)))),
    name = "Terms containing the target", guide = ggplot2::guide_legend(title.position = "left", nrow = 1,
      override.aes = list(shape = 16, colour = .tok$muted)))
  m_max <- max(table(factor(genes$cls, c("Inhibited", "Activated"))))
  panel_mm <- max(45, (m_max - 1) * 4.6, (k - 1) * 13)            # height of the [0, 1] range in mm
  if (layout == "signed") {
    tn$y <- if (k == 1) 0.5 else rev(seq(0, 1, length.out = k))
    genes$bary <- vapply(genes$gene, function(g) mean(tn$y[match(e$tid[e$gene == g], tn$tid)]), numeric(1))
    gap <- 1 / max(1, m_max - 1)
    place <- function(d) { d <- d[base::order(-round(d$bary, 9), -abs(d$score), d$gene), ]; d$y <- .spread(d$bary, gap); d }   # rounding: ties in the mean height
    genes <- rbind(place(genes[genes$cls == "Inhibited", ]), place(genes[genes$cls == "Activated", ]))
    genes$x <- genes$score
    bw <- 0.85; hb <- 1.25 / panel_mm                             # bar half-width (score units), half-height (y units)
    ee <- merge(e, genes[, c("gene", "x", "y", "cls")], by = "gene"); ee <- merge(ee, tn[, c("tid", "y")], by = "tid", suffixes = c("", "_t"))
    ee$x0 <- ifelse(ee$cls == "Inhibited", -bw, bw)
    bars <- rbind(data.frame(tid = tn$tid, xmin = -bw, xmax = -bw + 2 * bw * tn$f_inh, y = tn$y, cls = "Inhibited"),
                  data.frame(tid = tn$tid, xmin = -bw + 2 * bw * tn$f_inh, xmax = bw, y = tn$y, cls = "Activated"))
    xl <- c(min(-2, min(genes$x)) - 2.6, max(2, max(genes$x)) + 2.6)
    ylab_mm <- 3.2 / panel_mm
    ylo <- -hb - 7 / panel_mm; yhi <- 1 + 9 / panel_mm
    p <- ggplot2::ggplot() +
      ggforce::geom_diagonal(data = ee, ggplot2::aes(x = x0, y = y_t, xend = x, yend = y, colour = cls), linewidth = 0.3, alpha = 0.45,
        strength = 0.6) +
      ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - hb, ymax = y + hb, fill = NULL), fill = .pal$marks[bars$cls],
        colour = "white", linewidth = 0.25) +
      ggplot2::geom_label(data = tn, ggplot2::aes(x = 0, y = y + hb + 0.5 * ylab_mm, label = label), vjust = 0, size = 2.15, fontface = "bold",
        colour = .tok$ink, family = .pai_font(), fill = grDevices::adjustcolor("white", 0.88), linewidth = 0, label.padding = grid::unit(0.3, "mm")) +
      ggplot2::geom_label(data = tn, ggplot2::aes(x = 0, y = y - hb - 0.5 * ylab_mm, label = line2), parse = TRUE, vjust = 1, size = 1.85,
        colour = .tok$ink2, family = .pai_font(), fill = grDevices::adjustcolor("white", 0.88), linewidth = 0, label.padding = grid::unit(0.3, "mm")) +
      ggplot2::geom_point(data = genes, ggplot2::aes(x = x, y = y, fill = score, size = degree), shape = 21, colour = "white", stroke = 0.5) +
      ggplot2::geom_text(data = genes, ggplot2::aes(x = x + ifelse(cls == "Inhibited", -0.28, 0.28), y = y, label = gene,
        hjust = ifelse(cls == "Inhibited", 1, 0), fontface = ifelse(degree > 1, "bold", "plain")), size = 2.1, colour = .tok$ink, family = .pai_font()) +
      ggplot2::annotate("text", x = c(xl[1] + 0.2, xl[2] - 0.2), y = 1 + 5.5 / panel_mm, hjust = c(0, 1), vjust = 0, size = 2.3, fontface = "bold",
        label = c(paste("\u2190", sem$short[["Inhibited"]], "targets"), paste(sem$short[["Activated"]], "targets \u2192")),
        colour = c(.pal$text[["Inhibited"]], .pal$text[["Activated"]]), family = .pai_font()) +
      ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") + fill_scale + size_scale +
      ggplot2::scale_x_continuous(limits = xl, expand = c(0, 0), breaks = function(l) { b <- seq(-8, 8, 2); b[b != 0 & b > l[1] + 1.5 & b < l[2] - 1.5] },
        labels = .lab_minus) +
      ggplot2::scale_y_continuous(limits = c(ylo, yhi), expand = c(0, 0)) +
      ggplot2::labs(x = sem$score_label, y = NULL) +
      theme_pai(axes = "x") + ggplot2::theme(axis.text.y = ggplot2::element_blank(), legend.position = "bottom",
        legend.justification = "left", legend.box = "horizontal")
  } else {
    vert <- data.frame(name = c(tn$tid, genes$gene), type = c(rep("term", k), rep("gene", nrow(genes))))
    gr <- igraph::graph_from_data_frame(e[, c("tid", "gene")], directed = FALSE, vertices = vert)
    set.seed(seed); lay <- igraph::layout_with_fr(gr, niter = 2000)
    lay <- (lay - min(lay)) / diff(range(lay))
    nodes <- cbind(vert, x = lay[, 1], y = lay[, 2])
    tn$x <- nodes$x[match(tn$tid, nodes$name)]; tn$y <- nodes$y[match(tn$tid, nodes$name)]
    genes$x <- nodes$x[match(genes$gene, nodes$name)]; genes$y <- nodes$y[match(genes$gene, nodes$name)]
    ee <- merge(merge(e, genes[, c("gene", "x", "y", "cls")], by = "gene"), tn[, c("tid", "x", "y")], by = "tid", suffixes = c("", "_t"))
    bw <- 0.05; hb <- 0.012
    bars <- rbind(data.frame(tid = tn$tid, xmin = tn$x - bw, xmax = tn$x - bw + 2 * bw * tn$f_inh, y = tn$y, cls = "Inhibited"),
                  data.frame(tid = tn$tid, xmin = tn$x - bw + 2 * bw * tn$f_inh, xmax = tn$x + bw, y = tn$y, cls = "Activated"))
    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(data = ee, ggplot2::aes(x = x_t, y = y_t, xend = x, yend = y, colour = cls), linewidth = 0.3, alpha = 0.45) +
      ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = y - hb, ymax = y + hb), fill = .pal$marks[bars$cls],
        colour = "white", linewidth = 0.25) +
      ggplot2::geom_point(data = genes, ggplot2::aes(x = x, y = y, fill = score, size = degree), shape = 21, colour = "white", stroke = 0.5) +
      do.call(ggrepel::geom_text_repel, c(list(data = genes, mapping = ggplot2::aes(x, y, label = gene, fontface = ifelse(degree > 1, "bold", "plain")),
        size = 2.1, family = .pai_font(), colour = .tok$ink, segment.size = 0.2, segment.colour = .tok$axis, box.padding = 0.2,
        max.overlaps = Inf, seed = seed), .halo)) +
      ggplot2::geom_label(data = tn, ggplot2::aes(x = x, y = y + hb + 0.012, label = tid), vjust = 0, size = 2.3, fontface = "bold",
        colour = .tok$ink, family = .pai_font(), fill = grDevices::adjustcolor("white", 0.88), linewidth = 0, label.padding = grid::unit(0.3, "mm")) +
      ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") + fill_scale + size_scale +
      ggplot2::coord_equal(clip = "off") + theme_pai(axes = "none") +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.title = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
        legend.position = "bottom", legend.justification = "left", legend.box = "horizontal")
    ylo <- 0; yhi <- 1; panel_mm <- 110
    key <- tn; key$row <- max(k, 8) - seq_len(k) + 1                 # listed from the top
    key$txt <- sprintf("%s\n%d\u2193 %d\u2191  \u00b7  z %s", .short_label(sel$Description, 30), tn$n_inh, tn$n_act, .fmt_num(tn$z, 1))
    tab <- ggplot2::ggplot(key, ggplot2::aes(y = row)) +
      ggplot2::geom_text(ggplot2::aes(x = 0, label = tid), hjust = 0, size = 2.3, fontface = "bold", colour = .tok$ink, family = .pai_font()) +
      geom_pai_pill(ggplot2::aes(x = 0.55, y = row, label = db), data = key, hjust = 0, size = 1.75) +
      ggplot2::geom_text(ggplot2::aes(x = 1.75, label = txt), hjust = 0, size = 2.05, lineheight = 1.05, colour = .tok$ink, family = .pai_font()) +
      ggplot2::scale_x_continuous(limits = c(0, 7)) + ggplot2::scale_y_continuous(limits = c(0.4, max(k, 8) + 0.6)) + ggplot2::theme_void()
    p <- (p | tab) + patchwork::plot_layout(widths = c(1.7, 1))
  }
  chips <- .pai_chips(c(stats::setNames(.pal$marks[c("Inhibited", "Activated")], sprintf("%s target and its links", sem$short[c("Inhibited", "Activated")]))),
    shapes = "circle", size = 6)
  key <- .pai_badge("Term bar", sprintf("share of the term's called targets: Inhibited (cool, left) and Activated (warm, right)  \u00b7  z: %s",
    if (z == "balanced") sprintf("balanced, p = %.2f", tt$p_activated) else "(\u2191 \u2212 \u2193) / \u221an"))
  out <- .stack(list(.pai_header(title %||% "Enriched terms and their called targets, by direction",
      sprintf("%s; up to %d called targets per direction, most shared first; bold: in more than one term", .sel_note(sel, n, terms, max_overlap, direction), n_genes)),
    chips, p, key), c(.header_mm(1), 4.6, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 0.5, 4.8, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- list(terms = tn, genes = genes, edges = e, layout = layout)
  .with_size(out, width, round(.header_mm(1) + 0.5 + 4.8 + 6.5 + 24 + panel_mm * (yhi - ylo), 1))
}

# ---- 5. One-page report -----------------------------------------------------------------------
#' One-page report of an L1000PAI result
#'
#' Arm plane, top calls, evidence card and, when the enrichment is given, the activity ladder of
#' the top pathways (otherwise the running-sum curves of the strongest call).
#'
#' @inheritParams pai_plot_params
#' @param e Optional result of [pai_enrich()].
#' @param file Optional file (`.pdf` or `.png`) to write the page to.
#' @param label Targets to label and to put on the evidence card.
#' @param contrast Text describing the contrast, shown under the title.
#' @param n_terms Number of pathways on the ladder.
#' @param width,height Page size in millimetres.
#' @return The page (a patchwork object), invisibly.
#' @family target plots
#' @export
pai_report <- function(x, e = NULL, file = NULL, label = NULL, contrast = NULL, n_terms = 8, width = 180, height = 250) {
  df <- .pai_df(x); lib <- .pai_library(x); md <- pai_metadata(x)
  sem <- .pai_semantics(lib, .pai_score_direction(x, lib))
  nc <- table(factor(as.character(df$activity_class), levels = c("Activated", "Inhibited", "Discordant", "Not_significant")))
  called <- df[as.character(df$activity_class) %in% c("Activated", "Inhibited"), ]
  gene <- called$perturbagen %||% called$term
  top <- c(utils::head(gene[called$activity_class == "Activated"][order(-called$activity_score[called$activity_class == "Activated"])], 2),
           utils::head(gene[called$activity_class == "Inhibited"][order(called$activity_score[called$activity_class == "Inhibited"])], 2))
  ev <- unique(c(label, top))
  tiles <- .pai_tiles(list(
    list(label = "Genes ranked", value = .fmt_n(md$n_ranked_genes)),
    list(label = "Pairs scored", value = .fmt_n(nrow(df))),
    list(label = sem$short[["Activated"]], value = .fmt_n(nc[["Activated"]]), cls = "Activated"),
    list(label = sem$short[["Inhibited"]], value = .fmt_n(nc[["Inhibited"]]), cls = "Inhibited"),
    list(label = "Discordant", value = .fmt_n(nc[["Discordant"]]), cls = "Discordant")))
  lettered <- function(p, l) {
    tag <- .wrap_grob(grid::textGrob(l, x = grid::unit(0, "mm"), y = grid::unit(1, "npc") - grid::unit(2.3, "mm"),
      just = c("left", "centre"), gp = .gp(9, face = "bold")))
    inner <- patchwork::wrap_elements(full = p & ggplot2::theme(plot.margin = ggplot2::margin(1, 2, 1, 2))) +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 2, 0, 0))
    patchwork::wrap_elements(full = (tag | inner) + patchwork::plot_layout(widths = grid::unit(c(4.2, 1), c("mm", "null")))) +
      ggplot2::theme(plot.margin = ggplot2::margin(2, 1, 2, 1))
  }
  a <- lettered(pai_armplot(x, label = label, n_label = 4), "a")
  b <- lettered(pai_dotplot(x, n = 8, highlight = label), "b")
  c_ <- lettered(pai_evidence_card(x, targets = ev), "c")
  has_terms <- !is.null(e) && any(.as_tt(e)$table$p.adjust < 0.05)
  if (!is.null(e) && !has_terms)
    message("No pathway term has adj. P < 0.05; panel d shows the running sum of the strongest call instead of the ladder.")
  d_ <- if (has_terms) lettered(pai_ladder(e, n = n_terms), "d") else
    lettered(pai_gseaplot(x, called$term[which.max(called$activity_score)]), "d")
  body <- patchwork::wrap_plots(list(a, b, c_, d_), design = "AB\nCC\nDD", heights = c(1.2, 0.62, 0.95))
  head <- .pai_header("L1000PAI report", contrast %||% "", size = 11)
  pills <- .pai_pills(list(list("L1000PAIResult", "object"), list(sprintf("library = \"%s\"", lib), "fn"),
    list(sprintf("mode = %s", md$library_mode_used), "fn"), list(sprintf("seed = %s", md$gsea_params$rng_seed %||% "none"), "fn")),
    note = sprintf("L1000PAI %s", as.character(utils::packageVersion("L1000PAI"))))
  p <- patchwork::wrap_plots(list(head, pills, tiles, body), ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(9, 5.5, 13, 1), c("mm", "mm", "mm", "null")))
  p <- .with_size(p, width, height)
  if (!is.null(file)) pai_save(p, file)
  invisible(p)
}
