# L1000PAI plotting style.
# All plotting functions share one visual language:
#   * colour means direction only: WARM = activation side, COOL = inhibition side;
#     Discordant is a lavender diamond, not-called targets a warm grey;
#   * the activation quadrant of the oriented arm plane is tinted warm, the inhibition quadrant cool;
#   * class names are written in their class colour (bold); every other text is ink;
#   * databases, objects and functions are shown as pills (monospace), not as colours;
#   * L-shaped axes (left and bottom) with short ticks, no grid unless asked.
# pai_use_palette("classic") switches to a darker Discordant (#6E5A8A), which fails the
# colour-vision check against COOL; the diamond shape then carries the difference.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a

.pai_palettes <- list(
  house = list(
    marks = c(Activated = "#B0402E", Inhibited = "#2E6E9E", Discordant = "#A97FD0", Not_significant = "#C9C6C0"),
    text  = c(Activated = "#B0402E", Inhibited = "#2E6E9E", Discordant = "#6E5A8A", Not_significant = "#52514E"),
    tint  = c(Activated = "#F3E1DC", Inhibited = "#DDE8F1", Discordant = "#E9E4F0", Not_significant = "#EFEFED")),
  classic = list(
    marks = c(Activated = "#B0402E", Inhibited = "#2E6E9E", Discordant = "#6E5A8A", Not_significant = "#C9C6C0"),
    text  = c(Activated = "#B0402E", Inhibited = "#2E6E9E", Discordant = "#6E5A8A", Not_significant = "#52514E"),
    tint  = c(Activated = "#F3E1DC", Inhibited = "#DDE8F1", Discordant = "#E9E4F0", Not_significant = "#EFEFED")))

.pai_shapes <- c(Activated = 21, Inhibited = 21, Discordant = 23, Not_significant = 21)   # Discordant = diamond

.tok <- list(ink = "#141414", ink2 = "#52514E", muted = "#8C8C8C", rule = "#DCDCDC",
             axis = "#B9B7B1",            # reference lines (zero, thresholds)
             frame = "#52514E",           # axis lines and ticks
             grid = "#EEEDEA", pill = "#F5F5F3", code = "#F3F5F7", code_ink = "#5F6B75",
             surface = "#FFFFFF")
# Fonts: options(pai.font = , pai.mono = ). "sans" is Arial on Windows and Helvetica on macOS.
.pai_font <- function() getOption("pai.font", "sans")
.pai_mono <- function() getOption("pai.mono", "mono")
.pai_mid <- "#EDECE8"                       # neutral midpoint of the diverging scale
# Module and term identity (hub plots): hues kept away from the WARM/COOL direction colours;
# adjacent pairs pass the colour-vision check; every module is also labelled M1, M2, ...
.pai_module_colors <- c("#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6")
.pai_other_module <- "#C9C6C0"

.pal <- new.env(parent = emptyenv())

#' Choose the colour palette of the L1000PAI plots
#'
#' Colour means direction only: Activated (warm) and Inhibited (cool) targets; Discordant targets
#' are drawn as lavender diamonds and targets that are not called in warm grey. `"house"` (the
#' default) uses a light Discordant (#A97FD0) that passes the colour-vision check against the
#' Inhibited colour; `"classic"` uses a darker one (#6E5A8A), which does not, so only the diamond
#' shape tells it apart for some readers.
#'
#' @param name `"house"` or `"classic"`.
#' @return The name of the palette in use before the call, invisibly, so that it can be restored.
#' @family plot style
#' @export
pai_use_palette <- function(name = c("house", "classic")) {
  name <- match.arg(name); p <- .pai_palettes[[name]]
  old <- .pal$name
  .pal$name <- name
  .pal$marks <- p$marks
  .pal$text <- p$text
  .pal$tint <- p$tint
  # The two arms of a pair: up arm on the warm side, down arm on the cool side.
  .pal$arm <- c(up = p$marks[["Activated"]], dn = p$marks[["Inhibited"]])
  invisible(old)
}
pai_use_palette("house")

.tint <- function(col, f = 0.6) {                          # blend towards white
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(t(m + (1 - m) * f))
}

# axes: "x" (bottom line only), "xy" (L-shaped), "none"
#' ggplot2 theme of the L1000PAI plots
#'
#' L-shaped axes with short ticks and no grid unless asked. Add it to a ggplot of your own to
#' match the L1000PAI plots.
#'
#' @param base_size Base font size in points; default `getOption("pai.base_size", 7)`.
#' @param base_family Font family; default `getOption("pai.font", "sans")` (Arial on Windows,
#'   Helvetica on macOS).
#' @param grid Grid lines to draw: `"none"`, `"x"`, `"y"` or `"both"`.
#' @param axes Axis lines: `"x"` (bottom only), `"xy"` (L-shaped) or `"none"`.
#' @return A ggplot2 theme.
#' @family plot style
#' @export
theme_pai <- function(base_size = getOption("pai.base_size", 7), base_family = .pai_font(),
                      grid = c("none", "x", "y", "both"), axes = c("x", "xy", "none")) {
  grid <- match.arg(grid); axes <- match.arg(axes)
  gl <- ggplot2::element_line(colour = .tok$grid, linewidth = 0.25)
  al <- ggplot2::element_line(colour = .tok$frame, linewidth = 0.25, lineend = "square")
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(colour = .tok$ink),
      axis.text = ggplot2::element_text(colour = .tok$ink2, size = ggplot2::rel(0.9)),
      axis.title = ggplot2::element_text(colour = .tok$ink, size = ggplot2::rel(0.92)),
      axis.line.x = if (axes %in% c("x", "xy")) al else ggplot2::element_blank(),
      axis.line.y = if (axes == "xy") al else ggplot2::element_blank(),
      axis.ticks = ggplot2::element_line(colour = .tok$frame, linewidth = 0.25),
      axis.ticks.y = if (axes == "xy") ggplot2::element_line(colour = .tok$frame, linewidth = 0.25) else ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(0.7, "mm"),
      panel.grid = ggplot2::element_blank(),
      panel.grid.major.x = if (grid %in% c("x", "both")) gl else ggplot2::element_blank(),
      panel.grid.major.y = if (grid %in% c("y", "both")) gl else ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.1), colour = .tok$ink),
      plot.subtitle = ggplot2::element_text(colour = .tok$ink2),
      plot.caption = ggplot2::element_text(colour = .tok$muted, hjust = 0, size = ggplot2::rel(0.85)),
      plot.title.position = "plot", plot.caption.position = "plot",
      strip.text = ggplot2::element_text(face = "bold", colour = .tok$ink, hjust = 0),
      legend.title = ggplot2::element_text(colour = .tok$ink2, size = ggplot2::rel(0.88)),
      legend.text = ggplot2::element_text(colour = .tok$ink2, size = ggplot2::rel(0.86)),
      legend.key.size = ggplot2::unit(3, "mm"),
      plot.margin = ggplot2::margin(2, 4, 2, 4)
    )
}

# ---- grid furniture: headers, pills, glyph legends, badges, stat tiles -------------------------
.gp <- function(size, col = .tok$ink, face = "plain", family = .pai_font(), ...)
  grid::gpar(fontsize = size, col = col, fontface = face, fontfamily = family, ...)
.wrap_grob <- function(g) patchwork::wrap_elements(full = g) + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))

.chip_grob <- function(x, y, col, w = 2.3) grid::roundrectGrob(x = x, y = y, width = grid::unit(w, "mm"),
  height = grid::unit(w, "mm"), just = c("left", "centre"), r = grid::unit(0.55, "mm"), gp = grid::gpar(fill = col, col = NA))

# A marker glyph as drawn in the plots: circle, diamond, square or line.
.glyph_grob <- function(x, y, col, shape = "circle", w = 2.2) {
  cx <- x + grid::unit(w / 2, "mm")
  switch(shape,
    circle = grid::circleGrob(cx, y, r = grid::unit(w / 2, "mm"), gp = grid::gpar(fill = col, col = NA)),
    diamond = grid::polygonGrob(grid::unit.c(cx - grid::unit(w * 0.62, "mm"), cx, cx + grid::unit(w * 0.62, "mm"), cx),
      grid::unit.c(y, y + grid::unit(w * 0.62, "mm"), y, y - grid::unit(w * 0.62, "mm")), gp = grid::gpar(fill = col, col = NA)),
    line = grid::segmentsGrob(x, y, x + grid::unit(4, "mm"), y, gp = grid::gpar(col = col, lwd = 2.2, lineend = "round")),
    dash = grid::segmentsGrob(x, y, x + grid::unit(4, "mm"), y, gp = grid::gpar(col = col, lwd = 2.2, lty = "22", lineend = "butt")),
    .chip_grob(x, y, col, w))
}

# A pill: rounded box with text. kind = "object" (black, white mono), "fn" (pale, mono),
# "class" (class tint with bold class-coloured text), "db" (pale, mono, small).
.pill_spec <- function(kind, cls = NULL) switch(kind,
  object = list(fill = .tok$ink, col = "white", family = .pai_mono(), face = "plain"),
  fn = list(fill = .tok$code, col = .tok$ink, family = .pai_mono(), face = "plain"),
  db = list(fill = .tok$pill, col = .tok$ink2, family = .pai_mono(), face = "plain"),
  class = list(fill = .pal$tint[[cls]], col = .pal$text[[cls]], family = .pai_font(), face = "bold"))
.pill_grob <- function(text, x, y, kind = "fn", cls = NULL, size = 6, pad = 1.2, h = 3.8) {
  s <- .pill_spec(kind, cls)
  tg <- grid::textGrob(text, x = x + grid::unit(pad, "mm"), y = y, just = c("left", "centre"),
    gp = .gp(size, s$col, s$face, s$family))
  box <- grid::roundrectGrob(x = x, y = y, width = grid::grobWidth(tg) + grid::unit(2 * pad, "mm"), height = grid::unit(h, "mm"),
    just = c("left", "centre"), r = grid::unit(0.7, "mm"), gp = grid::gpar(fill = s$fill, col = NA))
  list(grob = grid::grobTree(box, tg), width = grid::grobWidth(tg) + grid::unit(2 * pad, "mm"))
}
# A row of pills: list(list(text, kind, cls), ...), optional trailing note in ink2.
.pai_pills <- function(items, note = NULL, size = 6, gap = 1.2) {
  g <- list(); x <- grid::unit(0, "mm"); y <- grid::unit(0.5, "npc")
  for (it in items) {
    p <- .pill_grob(it[[1]], x, y, it[[2]] %||% "fn", it[[3]], size)
    g[[length(g) + 1]] <- p$grob; x <- x + p$width + grid::unit(gap, "mm")
  }
  if (!is.null(note)) g[[length(g) + 1]] <- grid::textGrob(note, x = x + grid::unit(0.8, "mm"), y = y,
    just = c("left", "centre"), gp = .gp(size, .tok$ink2))
  .wrap_grob(do.call(grid::grobTree, g))
}

# Panel header: bold title (optional class pill in front) and a grey subtitle.
# subtitle_parse = TRUE: `subtitle` is a vector of plotmath strings, one per line (for P values as 3.5 x 10^-9).
.pai_header <- function(title, subtitle = NULL, chip = NULL, size = 8, cls = NULL, pill = NULL, subtitle_parse = FALSE) {
  y1 <- grid::unit(1, "npc") - grid::unit(2.2, "mm")
  g <- list(); x0 <- grid::unit(0, "mm")
  if (!is.null(cls)) {
    p <- .pill_grob(pill %||% cls, x0, y1, "class", cls, size * 0.82, h = 3.6)
    g[[1]] <- p$grob; x0 <- p$width + grid::unit(1.6, "mm")
  } else if (!is.null(chip)) {
    g[[1]] <- .chip_grob(grid::unit(0, "mm"), y1, chip); x0 <- grid::unit(3.4, "mm")
  }
  g[[length(g) + 1]] <- grid::textGrob(title, x = x0, y = y1, just = c("left", "centre"), gp = .gp(size, face = "bold"))
  if (!is.null(subtitle) && !subtitle_parse) g[[length(g) + 1]] <- grid::textGrob(subtitle, x = grid::unit(0, "mm"), y = y1 - grid::unit(2.7, "mm"),
    just = c("left", "top"), gp = .gp(size * 0.78, .tok$ink2, lineheight = 1.1))
  if (!is.null(subtitle) && subtitle_parse) for (i in seq_along(subtitle))
    g[[length(g) + 1]] <- grid::textGrob(parse(text = subtitle[i])[[1]], x = grid::unit(0, "mm"),
      y = y1 - grid::unit(2.7 + 3 * (i - 1), "mm"), just = c("left", "top"), gp = .gp(size * 0.78, .tok$ink2))
  .wrap_grob(do.call(grid::grobTree, g))
}
.header_mm <- function(subtitle_lines = 0) 4.6 + 3 * subtitle_lines

# Legend row: c(label = colour, ...); shapes = "circle", "diamond", "square" or "line" (recycled);
# counts = optional counts printed bold after each label.
.pai_chips <- function(items, size = 6.2, line = FALSE, per_row = Inf, shapes = NULL, counts = NULL) {
  if (line) shapes <- "line"
  shapes <- rep_len(shapes %||% "square", length(items))
  rows <- split(seq_along(items), ceiling(seq_along(items) / per_row))
  g <- list(); nr <- length(rows)
  for (r in seq_along(rows)) {
    x <- grid::unit(0, "mm"); y <- grid::unit(1 - (r - 0.5) / nr, "npc")
    for (i in rows[[r]]) {
      g[[length(g) + 1]] <- .glyph_grob(x, y, items[[i]], shapes[i])
      tx <- x + grid::unit(if (shapes[i] %in% c("line", "dash")) 5.2 else 3.3, "mm")
      tg <- grid::textGrob(names(items)[i], x = tx, y = y, just = c("left", "centre"), gp = .gp(size, .tok$ink))
      g[[length(g) + 1]] <- tg; tx <- tx + grid::grobWidth(tg)
      if (!is.null(counts)) {
        cg <- grid::textGrob(counts[i], x = tx + grid::unit(1.2, "mm"), y = y, just = c("left", "centre"), gp = .gp(size, .tok$ink, "bold"))
        g[[length(g) + 1]] <- cg; tx <- tx + grid::unit(1.2, "mm") + grid::grobWidth(cg)
      }
      x <- tx + grid::unit(4, "mm")
    }
  }
  .wrap_grob(do.call(grid::grobTree, g))
}
# The four classes with their glyphs and counts.
.class_chips <- function(sem, counts, classes = c("Activated", "Inhibited", "Discordant", "Not_significant"), per_row = Inf, size = 6.2) {
  .pai_chips(stats::setNames(.pal$marks[classes], sem$short[classes]), shapes = ifelse(classes == "Discordant", "diamond", "circle"),
    counts = .fmt_n(as.integer(counts[classes])), per_row = per_row, size = size)
}

# A rounded badge: small-caps label plus body text, e.g. the call rule.
.pai_badge <- function(label, body, size = 6) {
  lab <- grid::textGrob(toupper(label), x = grid::unit(2, "mm"), y = grid::unit(0.5, "npc"), just = c("left", "centre"),
    gp = .gp(size * 0.9, .tok$muted, face = "bold"))
  txt <- grid::textGrob(body, x = grid::unit(2, "mm") + grid::grobWidth(lab) + grid::unit(2.2, "mm"),
    y = grid::unit(0.5, "npc"), just = c("left", "centre"), gp = .gp(size, .tok$ink2))
  box <- grid::roundrectGrob(x = grid::unit(0, "mm"), y = grid::unit(0.5, "npc"), just = c("left", "centre"),
    width = grid::unit(6.2, "mm") + grid::grobWidth(lab) + grid::grobWidth(txt), height = grid::unit(4.2, "mm"),
    r = grid::unit(0.9, "mm"), gp = grid::gpar(fill = .tok$pill, col = NA))
  .wrap_grob(grid::grobTree(box, lab, txt))
}

# Stat tiles: list(list(label, value, cls = NULL), ...). A class tile is tinted and its value coloured.
.pai_tiles <- function(tiles, gap = 2.5) {
  n <- length(tiles); g <- list()
  w <- (grid::unit(1, "npc") - grid::unit(gap * (n - 1), "mm")) * (1 / n)
  for (i in seq_len(n)) {
    x <- (w + grid::unit(gap, "mm")) * (i - 1)
    t <- tiles[[i]]; cls <- t$cls
    g[[length(g) + 1]] <- grid::roundrectGrob(x = x, y = grid::unit(0.5, "npc"), width = w, height = grid::unit(1, "npc"),
      just = c("left", "centre"), r = grid::unit(0.9, "mm"),
      gp = grid::gpar(fill = if (is.null(cls)) .tok$pill else .pal$tint[[cls]], col = NA))
    lx <- x + grid::unit(2.4, "mm")
    g[[length(g) + 1]] <- grid::textGrob(t$label, x = lx, y = grid::unit(1, "npc") - grid::unit(3, "mm"),
      just = c("left", "centre"), gp = .gp(6, if (is.null(cls)) .tok$ink2 else .pal$text[[cls]], if (is.null(cls)) "plain" else "bold"))
    g[[length(g) + 1]] <- grid::textGrob(t$value, x = lx, y = grid::unit(2.3, "mm"), just = c("left", "bottom"),
      gp = .gp(12, if (is.null(cls)) .tok$ink else .pal$text[[cls]], face = "bold"))
  }
  .wrap_grob(do.call(grid::grobTree, g))
}

.stack <- function(parts, heights_mm) {        # NA height = flexible
  h <- ifelse(is.na(heights_mm), 1, heights_mm)
  units <- ifelse(is.na(heights_mm), "null", "mm")
  patchwork::wrap_plots(parts, ncol = 1) + patchwork::plot_layout(heights = grid::unit(h, units))
}

# Database pills inside ggplot: short monospace tags instead of colours.
.db_short <- c(GO_BP = "GO:BP", GO_MF = "GO:MF", GO_CC = "GO:CC", KEGG = "KEGG", Reactome = "REAC", Elsevier = "ELSV")
.db_tag <- function(collection) unname(ifelse(collection %in% names(.db_short), .db_short[collection], collection))
geom_pai_pill <- function(mapping, data, size = 1.9, kind = "db", ...) {
  s <- .pill_spec(kind)
  ggplot2::geom_label(mapping = mapping, data = data, inherit.aes = FALSE, size = size, family = s$family,
    colour = s$col, fill = s$fill, linewidth = 0, label.r = grid::unit(0.6, "mm"),
    label.padding = grid::unit(0.35, "mm"), ...)
}


#' Arguments shared by the L1000PAI plotting functions
#'
#' @param x An L1000PAIResult (from [run_l1000_pai()]).
#' @param e Result of [pai_enrich()] (its term table from [pai_term_table()] also works).
#' @param n Number of terms to draw: the top `n` non-redundant terms by adjusted P.
#' @param terms Instead of `n`: terms given by key (`"Reactome|DNA Replication"`,
#'   `"GO_BP|GO:0006260"`) or by name. A name found in several collections takes the most
#'   significant one, with a message.
#' @param direction `"both"`, `"Activated"` or `"Inhibited"`: keep the terms whose better
#'   direction is this one.
#' @param z `"balanced"` (default) or `"goplot"`; see [pai_term_table()].
#' @param cutoff Adjusted P cut-off for the terms.
#' @param max_overlap Terms whose called targets overlap an already chosen term by this Jaccard
#'   index or more are skipped; `Inf` keeps every term.
#' @param collections Collections to draw from, e.g. `"Reactome"`; `NULL` uses all.
#' @param n_genes Number of genes or targets to draw.
#' @param library Library of `x` when it cannot be read from the result.
#' @param seed Seed of the random layout or jitter.
#' @param width Width in millimetres the layout is planned for; the recommended height follows.
#' @param title Title; `NULL` for the default.
#' @name pai_plot_params
#' @keywords internal
NULL

#' Plotting the results of L1000PAI
#'
#' The plotting functions turn an L1000PAI result into figures in one visual language: colour
#' means direction only (Activated warm, Inhibited cool), Discordant targets are lavender
#' diamonds, and explanations are drawn inside the figure as badges and glyph legends. Every
#' function returns a ggplot or patchwork object with the values drawn in
#' `attr(p, "plot_data")` and a recommended size used by [pai_save()].
#'
#' A typical session:
#' ```
#' x <- run_l1000_pai(rank_input, library = "crispr")
#' pai_dotplot(x)                  # top calls
#' pai_armplot(x)                  # both arms of every target
#' pai_evidence_card(x, n = 20)    # result table with a two-arm glyph
#'
#' e <- pai_enrich(x)              # pathways of the Activated and Inhibited targets
#' pai_gobubble(e); pai_gocircle(e, n = 10); pai_armmap(e, n = 8); pai_ladder(e, n = 20)
#' h <- pai_enrich_hubs(e)         # hub genes of the Inhibited pathways
#' pai_hubplot(h)
#'
#' p <- pai_gocircle(e, n = 20)
#' pai_save(p, "gocircle_top20.pdf")
#' ```
#'
#' Target plots take the result `x`; pathway plots take `e <- pai_enrich(x)`; hub plots take
#' `h <- pai_enrich_hubs(e)`. `n` always means the top n, `terms` / `targets` choose by name.
#'
#' Save PDF files with [pai_save()] or another cairo device: the base `pdf()` device replaces
#' the arrows, dagger and greater-or-equal sign of the labels with dots.
#'
#' @seealso [pai_dotplot()], [pai_enrich()], [pai_gocircle()], [pai_save()], [pai_use_palette()],
#'   [pai_style_guide()]
#' @name pai_plotting
NULL
