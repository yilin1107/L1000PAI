# One-page style sheet of the L1000PAI plotting style.
#' One-page sheet of the plotting style
#'
#' Colours, glyphs, pills and the rules the plotting functions follow, drawn with the palette in
#' use.
#'
#' @param file Optional file to write the page to.
#' @param width,height Page size in millimetres.
#' @return The page, invisibly.
#' @family plot style
#' @export
pai_style_guide <- function(file = NULL, width = 180, height = 214) {
  W <- width; H <- height
  sem <- .pai_semantics("crispr", "down_minus_up")
  cls <- c("Activated", "Inhibited", "Discordant", "Not_significant")
  # ---- canvas in millimetres (y from the top) ----
  tx <- function(x, y, s, size = 6, col = .tok$ink, face = "plain", h = 0, fam = .pai_font(), v = 0.5)
    ggplot2::annotate("text", x = x, y = H - y, label = s, size = size / ggplot2::.pt, colour = col, fontface = face,
      hjust = h, vjust = v, family = fam)
  box <- function(x, y, w, h, fill, r = NULL) ggplot2::annotate("rect", xmin = x, xmax = x + w, ymin = H - y - h, ymax = H - y, fill = fill)
  pill <- function(x, y, s, kind, cl = NULL, size = 6) {
    sp <- .pill_spec(kind, cl)
    ggplot2::annotate("label", x = x, y = H - y, label = s, hjust = 0, size = size / ggplot2::.pt, family = sp$family, fontface = sp$face,
      colour = sp$col, fill = sp$fill, linewidth = 0, label.r = grid::unit(0.7, "mm"), label.padding = grid::unit(0.8, "mm"))
  }
  sec <- function(y, n, title, sub = NULL) list(tx(4, y, n, 9, face = "bold"), tx(10, y, title, 7.2, face = "bold"),
    if (!is.null(sub)) tx(10, y + 4, sub, 5.6, .tok$ink2))
  L <- list()
  L <- c(L, tx(4, 6, "L1000PAI plotting style", 12, face = "bold"),
    tx(4, 12, sprintf("Every plotting function of L1000PAI %s draws in this language. Palette in use: \"%s\".",
      as.character(utils::packageVersion("L1000PAI")), .pal$name), 6.4, .tok$ink2))
  # 1. colour means direction
  L <- c(L, sec(22, "1", "Colour means direction, and nothing else", "Marks, text and tints per class; Discordant is a diamond, so shape repeats the difference."))
  for (i in seq_along(cls)) {
    x0 <- 10 + (i - 1) * 42; y0 <- 31
    L <- c(L, box(x0, y0, 38, 17, .pal$tint[[cls[i]]]),
      ggplot2::annotate("point", x = x0 + 5, y = H - y0 - 5.2, shape = if (cls[i] == "Discordant") 23 else 21,
        fill = .pal$marks[[cls[i]]], colour = "white", size = 4.2, stroke = 0.5),
      tx(x0 + 9.5, y0 + 4.2, sem$short[[cls[i]]], 6.8, .pal$text[[cls[i]]], "bold"),
      tx(x0 + 9.5, y0 + 8.0, sprintf("mark %s", .pal$marks[[cls[i]]]), 5, .tok$ink2, fam = .pai_mono()),
      tx(x0 + 9.5, y0 + 11.0, sprintf("text %s", .pal$text[[cls[i]]]), 5, .tok$ink2, fam = .pai_mono()),
      tx(x0 + 9.5, y0 + 14.0, sprintf("tint %s", .pal$tint[[cls[i]]]), 5, .tok$ink2, fam = .pai_mono()))
  }
  L <- c(L, box(10, 50.5, 164, 11.5, .tok$pill),
    tx(12.5, 53.5, "COLOUR CHECK", 5, .tok$muted, "bold"),
    tx(29, 53.5, "WARM \u00b7 COOL \u00b7 Discordant #A97FD0: all checks pass (colour-vision \u0394E 11.2 protan; normal vision 19.1; contrast \u2265 3:1).", 5.6, .tok$ink2),
    tx(29, 58.1, "Darker Discordant #6E5A8A vs COOL: \u0394E 3.4 (deutan) and 9.0 (normal vision), both fail; pai_use_palette(\"classic\") uses it.", 5.6, .tok$ink2))
  # 2. two arms, one plane (inset added below)
  L <- c(L, sec(68, "2", "Two arms, one plane", "Each arm is signed so that right/up = toward activation; the activity score is the sum of the two."))
  # 3. pills
  L <- c(L, tx(96, 68, "3", 9, face = "bold"), tx(102, 68, "Pills name things", 7.2, face = "bold"),
    tx(102, 72, "Objects, functions and databases are monospace pills, never colours.", 5.6, .tok$ink2))
  L <- c(L, tx(102, 80, "object", 5.4, .tok$muted), pill(118, 80, "L1000PAIResult", "object"),
    tx(102, 87, "function", 5.4, .tok$muted), pill(118, 87, "pai_armmap()", "fn"), pill(140, 87, "result_table()", "fn"),
    tx(102, 94, "class", 5.4, .tok$muted), pill(118, 94, "Activated", "class", "Activated"), pill(134, 94, "Inhibited", "class", "Inhibited"),
    pill(149.5, 94, "Discordant", "class", "Discordant"),
    tx(102, 101, "database", 5.4, .tok$muted))
  xs <- 118
  for (d in c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC", "ELSV")) { L <- c(L, pill(xs, 101, d, "db", size = 5.4)); xs <- xs + nchar(d) * 1.25 + 4 }
  L <- c(L, tx(102, 109, "Module hues appear only where identity is needed (hub plots), always with an M# label:", 5.6, .tok$ink2))
  for (i in seq_along(.pai_module_colors)) L <- c(L, box(102 + (i - 1) * 12, 112, 10, 4, .pai_module_colors[i]),
    tx(102 + (i - 1) * 12 + 5, 118.5, sprintf("M%d", i), 5.2, .tok$ink2, h = 0.5))
  L <- c(L, box(162, 112, 10, 4, .pai_other_module), tx(167, 118.5, "other", 5.2, .tok$ink2, h = 0.5),
    tx(102, 123, "Their hues avoid the WARM/COOL directions; five pass the colour-vision check for every pair.", 5.2, .tok$muted))
  # 4. type and furniture
  L <- c(L, sec(133, "4", "Type, numbers and furniture", "Sans-serif (Arial on Windows; options(pai.font =) to change); class names in their colour; true minus signs; P values as a \u00d7 10^b."))
  L <- c(L, tx(10, 143, "Panel title", 8, face = "bold"), tx(10, 147.3, "grey subtitle: what is drawn and how many", 6.2, .tok$ink2),
    tx(10, 154, "\u22122.56", 6.4, fam = .pai_font()), tx(22, 154, "axis and table numbers use U+2212, not a hyphen", 5.6, .tok$ink2),
    ggplot2::annotate("text", x = 10, y = H - 160, label = '"3.5" %*% 10^-9', parse = TRUE, hjust = 0, size = 6.4 / ggplot2::.pt, family = .pai_font()),
    tx(30, 160, "P values in table format", 5.6, .tok$ink2))
  # 5. rules
  L <- c(L, sec(170, "5", "Rules the functions follow"))
  rules <- c("Colour = direction (Activated warm, Inhibited cool); Discordant = lavender diamond; not called = warm grey.",
    "The activation quadrant is tinted warm and the inhibition quadrant cool.",
    "Arms: colour shows which way an arm points; line type (solid up arm, dashed down arm) shows which arm it is.",
    "Explanations are drawn as badges and glyph legends inside the figure, not left to the caption.",
    "L-shaped axes with short ticks; no grid unless a value must be read across.",
    "Every function returns ggplot/patchwork and attaches the values drawn as attr(p, \"plot_data\").")
  for (i in seq_along(rules)) L <- c(L, ggplot2::annotate("point", x = 11.5, y = H - (177 + (i - 1) * 5.6), shape = 15, size = 1.2, colour = .tok$ink2),
    tx(14, 177 + (i - 1) * 5.6, rules[i], 5.8))
  canvas <- ggplot2::ggplot() + L + ggplot2::scale_x_continuous(limits = c(0, W), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, H), expand = c(0, 0)) + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
  # arm schematic
  lim <- 4
  ex <- data.frame(x = c(2.4, -2.0, -1.6, 2.2, 1.0), y = c(2.0, -2.4, 1.6, -1.9, 0.5),
    cls = c("Activated", "Inhibited", "Discordant", "Discordant", "Not_significant"), lab = c("both arms agree", "both arms agree", "", "", "too weak"))
  arm <- ggplot2::ggplot(ex, ggplot2::aes(x, y)) + .arm_canvas(lim, sem, 1, label_size = 1.95) +
    ggplot2::geom_segment(data = ex[1, ], ggplot2::aes(x = 0, y = 0, xend = x, yend = 0), colour = .pal$marks[["Activated"]], linewidth = 0.9,
      arrow = grid::arrow(length = grid::unit(1, "mm"), type = "closed")) +
    ggplot2::geom_segment(data = ex[1, ], ggplot2::aes(x = x, y = 0, xend = x, yend = y), colour = .pal$marks[["Activated"]], linewidth = 0.9, linetype = "22",
      arrow = grid::arrow(length = grid::unit(1, "mm"), type = "closed")) +
    ggplot2::geom_point(ggplot2::aes(fill = cls), shape = ifelse(ex$cls == "Discordant", 23, 21), size = 2.4, colour = .tok$ink, stroke = 0.35) +
    ggplot2::annotate("text", x = 1.2, y = -0.35, label = "up arm", size = 1.9, colour = .pal$text[["Activated"]], family = .pai_font()) +
    ggplot2::annotate("text", x = 2.55, y = 1.0, label = "down arm", size = 1.9, colour = .pal$text[["Activated"]], family = .pai_font(), hjust = 0) +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_x_continuous(breaks = c(-4, -2, 0, 2, 4), labels = .lab_minus) + ggplot2::scale_y_continuous(breaks = c(-4, -2, 0, 2, 4), labels = .lab_minus) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    ggplot2::labs(x = "Up arm, oriented  (\u2212NES Up)", y = "Down arm  (NES Down)") + theme_pai(base_size = 5.6, axes = "xy")
  ins <- function(p, x0, y0, x1, y1) patchwork::inset_element(p, left = x0 / W, bottom = 1 - y1 / H, right = x1 / W, top = 1 - y0 / H,
    align_to = "full", clip = FALSE)
  page <- canvas + ins(arm, 8, 76, 66, 131) +
    ins(.pai_badge("Call rule", "coherent adj. P < 0.05  \u00b7  both arms aligned  \u00b7  |activity| \u2265 1"), 96, 139, 176, 145) +
    ins(.class_chips(sem, c(Activated = 1200, Inhibited = 800, Discordant = 400, Not_significant = 3000), per_row = 2), 96, 147, 176, 156) +
    ins(.pai_tiles(list(list(label = "Pairs scored (example)", value = "5,400"), list(label = "Activated", value = "1,200", cls = "Activated"),
      list(label = "Inhibited", value = "800", cls = "Inhibited"))), 96, 158, 176, 169)
  page <- .with_size(page, W, H)
  if (!is.null(file)) pai_save(page, file)
  invisible(page)
}
