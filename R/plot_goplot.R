# GOplot-style displays of the direction-specific enrichment.
# Colour means direction only; databases are shown as monospace pills.
#
# GOplot (Walter et al., Bioinformatics 2015) summarises each enriched term by a z-score,
# (up - down) / sqrt(count), over the term's genes. Here the genes are the called targets of a
# gene-level library and "up"/"down" are the Activated and Inhibited calls; the gene-level value
# shown (GOplot's logFC) is the activity score.

.collection_colors <- c(GO_BP = "#2a78d6", GO_MF = "#eb6834", GO_CC = "#1baf7a", KEGG = "#eda100",
                        Reactome = "#e87ba4", Elsevier = "#008300")          # validated order
.collection_names <- c(GO_BP = "GO BP", GO_MF = "GO MF", GO_CC = "GO CC", KEGG = "KEGG", Reactome = "Reactome", Elsevier = "Elsevier")

# Term table for GOplot-style displays
# 
# Significance comes from the direction-specific ORA of pai_enrich(): ORA of all called targets
# together lets the two directions cancel (cell-cycle targets are over-represented among the
# Inhibited calls and under-represented among the Activated ones), so it is not used. For every
# term tested in either direction the table gives the adjusted P of each direction, the better
# of the two, the called targets in the term (union of both directions), their numbers by
# direction, the GOplot z-score (activated - inhibited) / sqrt(count), and a balanced z-score
# (activated - n p) / sqrt(n p (1 - p)), with p the share of Activated among all called targets,
# so that z = 0 means "as many Activated as expected from the two list sizes".
#' Term table of the direction-specific enrichment
#'
#' For every term tested in either direction: the adjusted P of each direction and the better of
#' the two, the called targets in the term, their numbers by direction, the GOplot z-score
#' `(activated - inhibited) / sqrt(count)` and the balanced z-score
#' `(activated - n p) / sqrt(n p (1 - p))`, with `p` the share of Activated among all called
#' targets. The balanced z is 0 when a term holds as many Activated targets as expected from the
#' sizes of the two call lists; the GOplot z leans toward the larger list.
#'
#' [pai_enrich()] already contains this table; the plotting functions take either object.
#'
#' @param e Result of [pai_enrich()].
#' @return A `pai_term_table` object whose `table` has one row per term.
#' @family pathway plots
#' @export
pai_term_table <- function(e) {
  if (inherits(e, "pai_term_table")) return(e)
  if (!inherits(e, "pai_enrichment")) stop("`e` must be the result of pai_enrich().", call. = FALSE)
  e$terms %||% .build_term_table(e)
}

.build_term_table <- function(e) {
  lib <- e$library
  df <- e$targets; df <- df[!is.na(df$activity_score), ]
  gene <- df$gene
  cls <- stats::setNames(as.character(df$activity_class), gene)
  score <- stats::setNames(df$activity_score, gene)
  t <- e$table; t$key <- paste(t$collection, t$ID, sep = "|")
  keys <- unique(t$key)
  get <- function(d, col) { z <- t[t$direction == d, ]; stats::setNames(z[[col]], z$key)[keys] }
  pa <- get("Activated", "p.adjust"); pi_ <- get("Inhibited", "p.adjust")
  ga <- get("Activated", "geneID"); gi <- get("Inhibited", "geneID")
  first <- t[!duplicated(t$key), ]; first <- first[match(keys, first$key), ]
  genes <- mapply(function(a, b) sort(unique(c(if (is.na(a)) character() else strsplit(a, "/", fixed = TRUE)[[1]],
                                                 if (is.na(b)) character() else strsplit(b, "/", fixed = TRUE)[[1]]))),
                  ga, gi, SIMPLIFY = FALSE)
  up <- vapply(genes, function(z) sum(cls[z] == "Activated"), numeric(1))
  dn <- vapply(genes, function(z) sum(cls[z] == "Inhibited"), numeric(1))
  pa1 <- ifelse(is.na(pa), 1, pa); pi1 <- ifelse(is.na(pi_), 1, pi_)
  # share of Activated among all called targets: the expected share in a term if direction were
  # unrelated to membership. The balanced z standardises the Activated count against it.
  p_act <- sum(cls == "Activated") / sum(cls %in% c("Activated", "Inhibited"))
  nn <- up + dn
  tab <- data.frame(collection = first$collection, ID = first$ID, Description = first$Description, key = keys,
    padj_activated = pa1, padj_inhibited = pi1, p.adjust = pmin(pa1, pi1),
    direction = ifelse(pa1 <= pi1, "Activated", "Inhibited"),
    count = lengths(genes), n_activated = up, n_inhibited = dn, zscore = (up - dn) / sqrt(lengths(genes)),
    zscore_balanced = (up - nn * p_act) / sqrt(nn * p_act * (1 - p_act)),
    mean_activity = vapply(genes, function(z) mean(score[z]), numeric(1)),
    geneID = vapply(genes, paste, character(1), collapse = "/"), stringsAsFactors = FALSE)
  rownames(tab) <- NULL
  structure(list(table = tab, scores = score, classes = cls, called = gene[cls %in% c("Activated", "Inhibited")],
    background = gene, library = lib, semantics = e$semantics, params = e$params, sets = e$sets, p_activated = p_act,
    score_direction = e$score_direction, threshold = e$threshold, targets = e$targets),
    class = "pai_term_table")
}

# Every pathway-level plot takes the enrichment (or its term table).
.as_tt <- function(e) {
  if (inherits(e, "pai_term_table")) return(e)
  if (inherits(e, "pai_enrichment")) return(e$terms %||% .build_term_table(e))
  stop("Give the result of pai_enrich(): e <- pai_enrich(x).", call. = FALSE)
}

#' @export
print.pai_term_table <- function(x, ...) {
  cat(sprintf("<pai_term_table>  %s terms tested in %s; %s with adj. P < 0.05\n", .fmt_n(nrow(x$table)),
    paste(unique(x$table$collection), collapse = ", "), .fmt_n(sum(x$table$p.adjust < 0.05))))
  invisible(x)
}

# Non-redundant pick, as GOplot's reduce_overlap: walk terms by adjusted P and skip any term whose
# called-target set overlaps an already chosen one by Jaccard >= max_overlap (Inf keeps every term).
# direction = "Activated" or "Inhibited" keeps the terms whose better ORA direction is that one.
.pick_terms <- function(tt, n = 10, cutoff = 0.05, max_overlap = 0.5, collections = NULL, direction = NULL) {
  t <- tt$table[tt$table$p.adjust < cutoff, ]
  if (!is.null(collections)) t <- t[t$collection %in% collections, ]
  if (!is.null(direction) && !identical(direction, "both")) t <- t[t$direction == direction, ]
  t <- t[order(t$p.adjust, t$key), ]
  pick <- integer(); gl <- strsplit(t$geneID, "/", fixed = TRUE)
  for (i in seq_len(nrow(t))) {
    ok <- all(vapply(pick, function(j) {
      a <- gl[[i]]; b <- gl[[j]]; length(intersect(a, b)) / length(union(a, b)) < max_overlap
    }, logical(1)))
    if (ok) pick <- c(pick, i)
    if (length(pick) == n) break
  }
  out <- t[pick, ]; out$tid <- sprintf("T%d", seq_len(nrow(out))); rownames(out) <- NULL
  out
}

.zcol <- function(z) switch(z, goplot = list(col = "zscore", name = "z-score"), balanced = list(col = "zscore_balanced", name = "balanced z"))
.go_badge <- function(tt, z = "balanced", short = FALSE) .pai_badge("Terms", if (short) {
  if (identical(z, "balanced")) sprintf("ORA by direction  \u00b7  balanced z,  p = %.2f", tt$p_activated) else
    "ORA by direction  \u00b7  z = (\u2191 \u2212 \u2193) / \u221an"
} else if (identical(z, "balanced"))
  sprintf("direction-specific ORA  \u00b7  balanced z = (\u2191 \u2212 n\u00b7p) / \u221a(n\u00b7p\u00b7(1 \u2212 p)),  p = %.2f = share of Activated among all calls", tt$p_activated) else
  "direction-specific ORA  \u00b7  z = (\u2191 \u2212 \u2193) / \u221an  (GOplot)")
.dir_counts <- function(n_inh, n_act) sprintf("%d\u2193 %d\u2191", n_inh, n_act)
# Term key used beside several plots: ID, database pill, name, counts and q. Laid out in millimetres
# for a key `key_mm` wide: names are cut to the room left and q is dropped when the key is narrow.
.term_key <- function(sel, key_mm, q = key_mm >= 62, top = NULL, swatch = NULL) {
  k <- nrow(sel); top <- top %||% k; o <- if (is.null(swatch)) 0 else 3
  x_cnt <- key_mm - if (q) 17.5 else 1
  name_width <- max(8, floor((x_cnt - 10 - (o + 13.5)) / 1.12))
  key <- sel; key$y <- top - seq_len(k) + 1; key$db <- .db_tag(key$collection); key$lab <- .short_label(key$Description, name_width)
  key$cnt <- .dir_counts(key$n_inhibited, key$n_activated); key$q <- .fmt_p10_expr(key$p.adjust)
  g <- ggplot2::ggplot(key, ggplot2::aes(y = y))
  if (!is.null(swatch)) g <- g + ggplot2::geom_point(ggplot2::aes(x = 0.7, fill = tid), shape = 22, size = 2.3, colour = "white", stroke = 0.3) +
    ggplot2::scale_fill_manual(values = swatch, guide = "none")
  g <- g + ggplot2::geom_text(ggplot2::aes(x = o, label = tid), hjust = 0, size = 2.25, fontface = "bold", colour = .tok$ink, family = .pai_font()) +
    geom_pai_pill(ggplot2::aes(x = o + 5.2, y = y, label = db), data = key, hjust = 0, size = 1.7) +
    ggplot2::geom_text(ggplot2::aes(x = o + 13.5, label = lab), hjust = 0, size = 2.05, colour = .tok$ink, family = .pai_font()) +
    ggplot2::geom_text(ggplot2::aes(x = x_cnt, label = cnt), hjust = 1, size = 1.95, colour = .tok$ink2, family = .pai_font())
  if (q) g <- g + ggplot2::geom_text(ggplot2::aes(x = x_cnt + 1.2, label = q), parse = TRUE, hjust = 0, size = 1.95, colour = .tok$ink2, family = .pai_font())
  g + ggplot2::scale_x_continuous(limits = c(-0.8, key_mm), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0.4, top + 0.6), expand = c(0, 0)) + ggplot2::theme_void()
}
.div_scale <- function(lim, name, ...) ggplot2::scale_fill_gradient2(low = .pal$marks[["Inhibited"]], mid = .pai_mid,
  high = .pal$marks[["Activated"]], limits = c(-lim, lim), name = name, labels = .lab_minus, ...)

# Every display below takes the same term choice as the other functions: `n` (top n non-redundant
# terms by adjusted P), or `terms` (keys or names), `direction`, `collections`, `cutoff`,
# `max_overlap` (Inf keeps redundant terms); and `z = "balanced"` or `"goplot"`.

# GOBubble ------------------------------------------------------------------------------------
# Every tested term by z-score and significance. The top `n_label` non-redundant terms are
# numbered on the plot and listed beside it (label = "id"), or named on the plot (label = "name").
# facet = TRUE draws one panel per database.
#' Every tested term by z-score and significance (GOBubble)
#'
#' The top `n_label` non-redundant terms are numbered on the plot and listed beside it
#' (`label = "id"`) or named on the plot (`label = "name"`).
#'
#' @inheritParams pai_plot_params
#' @param n_label Number of terms to number or name.
#' @param label `"id"` or `"name"`.
#' @param facet `TRUE` draws one panel per collection.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_gobubble <- function(e, n_label = 10, label = c("id", "name"), facet = FALSE, z = c("balanced", "goplot"),
                         cutoff = 0.05, max_overlap = 0.5, collections = NULL, width = 180, title = NULL) {
  tt <- .as_tt(e)
  label <- match.arg(label); z <- match.arg(z); zc <- .zcol(z)
  t <- tt$table; if (!is.null(collections)) t <- t[t$collection %in% collections, ]
  t$zval <- t[[zc$col]]; t$neglogq <- .neglog10(t$p.adjust, 1e-30); t$sig <- t$p.adjust < cutoff
  t$db <- factor(.db_tag(t$collection), unique(.db_short))
  lab <- if (n_label > 0 && any(t$sig)) .select_terms(tt, NULL, n_label, cutoff, max_overlap, collections) else
    transform(t[0, ], tid = character())
  lab$zval <- lab[[zc$col]]; lab$neglogq <- .neglog10(lab$p.adjust, 1e-30); lab$db <- factor(.db_tag(lab$collection), levels(t$db))
  lab$name <- sprintf("%s \u00b7 %s", .short_label(lab$Description, 26), .db_tag(lab$collection))
  thr <- -log10(cutoff); xl <- max(abs(t$zval), 2) * 1.05
  p <- ggplot2::ggplot(t, ggplot2::aes(zval, neglogq)) +
    ggplot2::geom_hline(yintercept = thr, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::geom_point(data = t[!t$sig, ], ggplot2::aes(size = count), colour = "#E4E2DD", alpha = 0.7) +
    ggplot2::geom_point(data = t[t$sig, ], ggplot2::aes(size = count, fill = direction), shape = 21, colour = "white", stroke = 0.45, alpha = 0.9) +
    ggplot2::annotate("text", x = -xl, y = thr, label = sprintf("adj. P = %s", cutoff), hjust = 0, vjust = -0.5, size = 2.0,
      colour = .tok$muted, family = .pai_font()) +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_size_continuous(range = c(0.8, 6), name = "Called targets in the term", breaks = scales::breaks_pretty(3),
      guide = ggplot2::guide_legend(title.position = "left", nrow = 1, override.aes = list(colour = .tok$muted))) +
    ggplot2::scale_x_continuous(limits = c(-xl, xl), labels = .lab_minus) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.12))) +
    ggplot2::labs(x = sprintf("%s  (\u2190 more %s  \u00b7  more %s \u2192)", zc$name, tolower(tt$semantics$short[["Inhibited"]]),
      tolower(tt$semantics$short[["Activated"]])), y = "\u2212log10 adjusted P") +
    theme_pai(axes = "xy") + ggplot2::theme(legend.position = "bottom", legend.justification = "left")
  if (nrow(lab)) p <- p + if (label == "name")
    do.call(ggrepel::geom_text_repel, c(list(data = lab, mapping = ggplot2::aes(label = name), size = 2.1, family = .pai_font(),
      colour = .tok$ink, min.segment.length = 0, segment.size = 0.25, segment.colour = .tok$ink2, box.padding = 0.3,
      nudge_x = 1.6, hjust = 0, direction = "y", max.overlaps = Inf, seed = 1), .halo)) else
    do.call(ggrepel::geom_text_repel, c(list(data = lab, mapping = ggplot2::aes(label = tid), size = 2.2, fontface = "bold",
      family = .pai_font(), colour = .tok$ink, min.segment.length = 0.15, segment.size = 0.25, segment.colour = .tok$ink2,
      box.padding = 0.25, point.padding = 0.15, max.overlaps = Inf, seed = 1), .halo))
  if (facet) p <- p + ggplot2::facet_wrap(~db, ncol = 3, drop = TRUE) +
    ggplot2::theme(strip.text = ggplot2::element_text(family = .pai_mono(), size = 6.2, colour = .tok$ink2, face = "plain"),
      panel.spacing = ggplot2::unit(3, "mm"))
  below <- width < 130
  body <- if (label == "id" && nrow(lab)) {
    if (below) patchwork::wrap_elements(full = (p / .term_key(lab, width - 6)) +
      patchwork::plot_layout(heights = grid::unit(c(1, nrow(lab) * 4 + 1), c("null", "mm"))))
    else patchwork::wrap_elements(full = (p | .term_key(lab, width / (if (facet) 2.6 else 2.15) - 4, top = max(nrow(lab), 10))) +
      patchwork::plot_layout(widths = c(if (facet) 1.6 else 1.15, 1)))
  } else p
  chips <- .pai_chips(c(stats::setNames(.pal$marks[c("Inhibited", "Activated")], sprintf(if (below) "%s" else if (cutoff <= 0.05) "significant among %s targets" else
      paste0("adj. P < ", cutoff, " among %s targets"),
    tt$semantics$short[c("Inhibited", "Activated")])), stats::setNames("#E4E2DD", sprintf("adj. P \u2265 %s", cutoff))),
    shapes = "circle", size = 6)
  out <- .stack(list(.pai_header(title %||% "Enriched terms by direction and significance",
      sprintf("One bubble per tested term; %s%s", if (nrow(lab)) sprintf("numbered: top %d non-redundant terms", nrow(lab)) else
        sprintf("no term has adj. P < %s, so none is numbered", cutoff), if (facet) "; one panel per database" else "")),
    chips, body, .go_badge(tt, z, short = below)), c(.header_mm(1), 5, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 5, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- t[, c("key", "zscore", "zscore_balanced", "p.adjust", "count", "sig", "direction")]
  attr(out, "labelled") <- lab[, c("tid", "key")]
  .with_size(out, width, (if (facet) 175 else 120) + (if (below && label == "id") nrow(lab) * 4 + 1 else 0))
}

# GOBar ---------------------------------------------------------------------------------------
# Top terms as bars of -log10 adjusted P, coloured by z-score, with the Inhibited / Activated
# counts at the end of each bar. order = "z" puts the most Activated-leaning term on top;
# "padj" the most significant. by = "collection" splits the bars by database.
#' Top terms as bars (GOBar)
#'
#' Bars of -log10 adjusted P coloured by z-score, with the numbers of Inhibited and Activated
#' targets at the end of each bar.
#'
#' @inheritParams pai_plot_params
#' @param order `"z"` puts the most Activated-leaning term on top, `"padj"` the most significant.
#' @param by `"collection"` splits the bars by collection.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_gobar <- function(e, n = 15, terms = NULL, direction = c("both", "Activated", "Inhibited"), z = c("balanced", "goplot"),
                      order = c("z", "padj"), by = c("none", "collection"), cutoff = 0.05, max_overlap = 0.5, collections = NULL,
                      width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction); z <- match.arg(z); order <- match.arg(order); by <- match.arg(by); zc <- .zcol(z)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  k <- nrow(sel)
  sel$zval <- sel[[zc$col]]; sel$neglogq <- .neglog10(sel$p.adjust, 1e-30)
  sel$lab <- sprintf("%s  %s", sel$tid, .short_label(sel$Description, 46)); sel$db <- .db_tag(sel$collection)
  rk <- if (order == "z") base::order(sel$zval, -sel$p.adjust) else base::order(-sel$p.adjust, sel$key)
  sel$lab <- factor(sel$lab, sel$lab[rk])                           # first level at the bottom
  zl <- max(abs(sel$zval)); xm <- max(sel$neglogq)
  sel$cnt <- .dir_counts(sel$n_inhibited, sel$n_activated)
  p <- ggplot2::ggplot(sel, ggplot2::aes(y = lab)) +
    ggplot2::geom_col(ggplot2::aes(x = neglogq, fill = zval), width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(x = neglogq + xm * 0.015, label = cnt), hjust = 0, size = 2.0, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::geom_vline(xintercept = -log10(cutoff), colour = .tok$muted, linewidth = 0.25, linetype = "22") +
    .div_scale(zl, zc$name, guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(24, "mm"), barheight = ggplot2::unit(2, "mm"),
      title.position = "left")) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    ggplot2::labs(x = "\u2212log10 adjusted P", y = NULL) +
    theme_pai(axes = "x", grid = "none") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(colour = .tok$ink, size = ggplot2::rel(0.92)),
      legend.position = "bottom", legend.justification = "left")
  if (by == "collection") p <- p + ggplot2::facet_grid(db ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0, family = .pai_mono(), size = 6, colour = .tok$ink2, face = "plain", hjust = 0),
      panel.spacing.y = ggplot2::unit(1.5, "mm"))
  else {
    left <- ggplot2::ggplot(sel, ggplot2::aes(y = lab)) +
      geom_pai_pill(ggplot2::aes(x = 0, y = lab, label = db), data = sel, hjust = 0, size = 1.7) +
      ggplot2::geom_text(ggplot2::aes(x = 0.16, label = lab), hjust = 0, size = 2.1, colour = .tok$ink, family = .pai_font()) +
      ggplot2::scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) + ggplot2::theme_void() +
      ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0))
    p <- (left | (p + ggplot2::theme(axis.text.y = ggplot2::element_blank()))) + patchwork::plot_layout(widths = c(1.15, 2))
  }
  out <- .stack(list(.pai_header(title %||% sprintf("Top terms by %s", if (order == "z") "direction" else "significance"),
      sprintf("%s; bar length: \u2212log10 adj. P; colour: %s; dashed: adj. P = %s", .sel_note(sel, n, terms, max_overlap, direction), zc$name, cutoff)),
    p, .go_badge(tt, z)), c(.header_mm(1), NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 1, 6.5), c("mm", "null", "mm")))
  sel$row <- match(as.character(sel$lab), levels(sel$lab))          # 1 = bottom bar
  attr(out, "plot_data") <- sel[, c("tid", "key", "neglogq", "zscore", "zscore_balanced", "zval", "n_inhibited", "n_activated", "row")]
  .with_size(out, width, round(30 + k * 5.2 + (if (by == "collection") 2 * length(unique(sel$db)) else 0), 1))
}

# GOCircle ------------------------------------------------------------------------------------
# Top terms on a circle. From the centre outwards:
#   bars       -log10 adjusted P (reference circles: adj. P = 0.05 and round values), coloured by z;
#   split ring the term's called targets that are Inhibited (cool) and Activated (warm);
#   outer ring every member at its activity score, sorted from lowest to highest within the
#              term, with the zero line and |score| = threshold drawn.
# A gap at 12 o'clock carries two radial axes: -log10 adj. P for the bars (dashed tick: the cut-off) and
# activity score for the outer ring (dashed ticks: +/-threshold). members = "all" also draws the scored members that were
# not called (grey); label_genes = k names the k members with the largest |score| in each term.
#' Top terms on a circle (GOCircle)
#'
#' From the centre outwards: bars of -log10 adjusted P coloured by z-score; a ring split into the
#' shares of Inhibited and Activated called targets; and every member at its activity score,
#' sorted from lowest to highest within the term. Two radial axes in the gap at the top give the
#' scales.
#'
#' @inheritParams pai_plot_params
#' @param members `"all"` (default) draws every scored member of the term, the members that were
#'   not called in grey; `"called"` draws the called members only.
#' @param label_genes Name the `label_genes` lowest-scoring Inhibited and highest-scoring
#'   Activated members of each term (default 1; 0 names none).
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_gocircle <- function(e, n = 10, terms = NULL, direction = c("both", "Activated", "Inhibited"), z = c("balanced", "goplot"),
                         members = c("all", "called"), label_genes = 1, cutoff = 0.05, max_overlap = 0.5, collections = NULL,
                         width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction); z <- match.arg(z); members <- match.arg(members); zc <- .zcol(z)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  k <- nrow(sel); sel$zval <- sel[[zc$col]]; sel$neglogq <- .neglog10(sel$p.adjust, 1e-30)
  sel$f_inh <- sel$n_inhibited / (sel$n_inhibited + sel$n_activated)
  mem <- if (members == "all") .term_members(tt, sel$key) else lapply(sel$geneID, function(s) strsplit(s, "/", fixed = TRUE)[[1]])
  thr <- tt$threshold %||% 1
  # geometry on the angle axis: an axis gap of width `ag` at 12 o'clock, then k sectors of width 1
  ag <- 2.1; gp <- 0.14; X <- ag + k * (1 + gp)                 # the gap holds two radial axes
  sel$x0 <- ag + (seq_len(k) - 1) * (1 + gp) + gp / 2; sel$x1 <- sel$x0 + 1; sel$xm <- sel$x0 + 0.5
  genes <- do.call(rbind, lapply(seq_len(k), function(i) {
    g <- mem[[i]]; d <- data.frame(tid = sel$tid[i], i = i, gene = g, score = unname(tt$scores[g]), cls = unname(tt$classes[g]),
      stringsAsFactors = FALSE)
    d <- d[!is.na(d$score), ]; d <- d[base::order(d$score, d$gene), ]; m <- nrow(d)
    d$rank <- seq_len(m); d$x <- sel$x0[i] + 0.05 + 0.9 * (if (m == 1) 0.5 else (d$rank - 1) / (m - 1)); d
  }))
  # radii
  r_hole <- 1.0; r_bar <- c(r_hole, 2.45); r_split <- c(2.6, 2.74); r_ring <- c(2.96, 4.24); r_mid <- mean(r_ring); half <- diff(r_ring) / 2 * 0.92
  qmax <- max(sel$neglogq, -log10(cutoff)) * 1.05
  rq <- function(v) r_bar[1] + diff(r_bar) * v / qmax
  smax <- max(abs(genes$score)); rs <- function(v) r_mid + half * v / smax
  genes$r <- rs(genes$score)
  qtick <- pretty(c(0, qmax), 4); qtick <- qtick[qtick <= qmax * 0.999]      # axis ticks, -log10 adj. P
  if (length(qtick) > 5) qtick <- qtick[c(TRUE, FALSE)]
  qref <- c(if (cutoff < 1) -log10(cutoff), qtick[qtick > 0])                              # circles: the cut-off (dashed) and the ticks
  stick <- pretty(c(-smax, smax), 4); stick <- stick[abs(stick) <= smax]     # axis ticks, activity score
  # radial axes in the gap; labels sit a fixed distance from the axis, so the angle offset shrinks with radius
  ax <- ag / 2; off <- function(r, d) d / r * X / (2 * pi)
  axis_q <- data.frame(r = rq(qtick), lab = .lab_minus(qtick))
  axis_s <- data.frame(r = rs(stick), lab = ifelse(stick > 0, sprintf("+%s", stick), .lab_minus(stick)))
  axr <- function(d, side, gap = 0.16) transform(d, x = ax + side * off(r, gap), x0 = ax + side * off(r, 0.09), x1 = ax)
  zl <- max(abs(sel$zval))
  called <- genes$cls %in% c("Activated", "Inhibited")
  # label_genes = k: the k lowest-scoring Inhibited and the k highest-scoring Activated members of each term
  lab <- if (label_genes > 0) do.call(rbind, lapply(split(genes[called, ], genes$tid[called]), function(d) rbind(
    utils::head(d[d$cls == "Inhibited", ][base::order(d$score[d$cls == "Inhibited"], d$gene[d$cls == "Inhibited"]), ], label_genes),
    utils::head(d[d$cls == "Activated", ][base::order(-d$score[d$cls == "Activated"], d$gene[d$cls == "Activated"]), ], label_genes)))) else genes[0, ]
  if (nrow(lab)) {
    lab$ang <- 90 - 360 * lab$x / X; lab$flip <- lab$ang < -90
    lab$r_lab <- r_ring[2] + 0.08
  }
  r_tid <- r_ring[2] + if (nrow(lab)) 0.95 else 0.32
  sect <- sel
  qa <- axr(axis_q, -1); sa <- axr(axis_s, -1)
  qc <- if (cutoff < 1) axr(data.frame(r = rq(-log10(cutoff)), lab = as.character(cutoff)), 1) else NULL
  sc <- axr(data.frame(r = rs(c(-thr, thr)), lab = ""), 1)
  sc_lab <- axr(data.frame(r = r_mid, lab = sprintf("\u00b1%s", thr)), 1)
  p <- ggplot2::ggplot() +
    # scale references
    ggplot2::geom_segment(data = data.frame(r = rq(qref)), ggplot2::aes(x = ag, xend = X, y = r, yend = r), colour = .tok$rule, linewidth = 0.25,
      linetype = ifelse(qref == -log10(cutoff), "22", "solid")) +
    ggplot2::geom_rect(data = sect, ggplot2::aes(xmin = x0, xmax = x1, ymin = r_ring[1], ymax = r_ring[2]), fill = .tok$pill) +
    ggplot2::geom_segment(data = sect, ggplot2::aes(x = x0, xend = x1, y = r_mid, yend = r_mid), colour = .tok$ink2, linewidth = 0.25) +
    ggplot2::geom_segment(data = merge(sect[, c("tid", "x0", "x1")], data.frame(r = rs(c(-thr, thr)))),
      ggplot2::aes(x = x0, xend = x1, y = r, yend = r), colour = .tok$muted, linewidth = 0.22, linetype = "22") +
    # bars and split ring
    ggplot2::geom_rect(data = sect, ggplot2::aes(xmin = x0 + 0.06, xmax = x1 - 0.06, ymin = r_bar[1], ymax = rq(neglogq), fill = zval)) +
    .div_scale(zl, zc$name, guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(22, "mm"), barheight = ggplot2::unit(2, "mm"),
      title.position = "left")) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(data = sect, ggplot2::aes(xmin = x0, xmax = x0 + f_inh, ymin = r_split[1], ymax = r_split[2]), fill = .pal$marks[["Inhibited"]]) +
    ggplot2::geom_rect(data = sect, ggplot2::aes(xmin = x0 + f_inh, xmax = x1, ymin = r_split[1], ymax = r_split[2]), fill = .pal$marks[["Activated"]]) +
    # members, sorted by activity score
    ggplot2::geom_point(data = genes[!called, ], ggplot2::aes(x = x, y = r), colour = .pal$marks[["Not_significant"]], size = 0.45, stroke = 0) +
    ggplot2::geom_point(data = genes[called, ], ggplot2::aes(x = x, y = r, colour = cls), size = 0.75, stroke = 0) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::geom_text(data = sect, ggplot2::aes(x = xm, y = r_tid, label = tid), size = 2.5, fontface = "bold", colour = .tok$ink, family = .pai_font()) +
    # two radial axes in the gap: -log10 adj. P for the bars, activity score for the outer ring.
    # Ticks and numbers on the left; the dashed references (cut-off, \u00b1threshold) marked on the right.
    ggplot2::annotate("segment", x = ax, xend = ax, y = c(r_bar[1], rs(-smax)), yend = c(r_bar[2], rs(smax)), colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::geom_segment(data = rbind(qa, sa), ggplot2::aes(x = x0, xend = x1, y = r, yend = r), colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::geom_text(data = rbind(qa, sa), ggplot2::aes(x = x, y = r, label = lab), hjust = 1, size = 1.8, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::geom_segment(data = rbind(qc, sc), ggplot2::aes(x = x0, xend = x1, y = r, yend = r), colour = .tok$muted, linewidth = 0.3) +
    ggplot2::geom_text(data = rbind(qc, sc_lab), ggplot2::aes(x = x, y = r, label = lab), hjust = 0, size = 1.65, colour = .tok$muted,
      family = .pai_font()) +
    ggplot2::annotate("text", x = ax, y = c(r_bar[2], rs(smax)) + 0.1, label = c("\u2212log10 adj. P", "activity score"), vjust = 0, size = 1.7,
      colour = .tok$ink2, family = .pai_font()) +
    ggplot2::scale_x_continuous(limits = c(0, X), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, r_tid + 0.3), expand = c(0, 0)) +
    ggplot2::coord_polar(theta = "x", start = -2 * pi * ax / X, clip = "off") + ggplot2::theme_void(base_family = .pai_font()) +
    # coord_polar leaves a 10% band round the circle: the colour bar goes in the band below it, and the empty band above is trimmed
    ggplot2::theme(legend.position = "inside", legend.position.inside = c(0.5, 0), legend.justification = c(0.5, 0), legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 6.4, colour = .tok$ink2), legend.text = ggplot2::element_text(size = 6, colour = .tok$ink2),
      legend.margin = ggplot2::margin(0, 0, 0, 0), plot.margin = ggplot2::margin(-5, 2, 0, 2, "mm"))
  if (nrow(lab)) {
    lab$ang <- 90 - 360 * (lab$x - ax) / X; lab$flip <- lab$ang < -90
    p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x = x, y = r_lab, label = gene, angle = ifelse(flip, ang + 180, ang),
      hjust = ifelse(flip, 1, 0), colour = cls), size = 1.85, fontface = "bold", family = .pai_font(), show.legend = FALSE)
  }
  body <- (p | .term_key(sel, width / 2.45 - 4)) + patchwork::plot_layout(widths = c(1.45, 1))
  chips <- .pai_chips(c(stats::setNames(.pal$marks[c("Inhibited", "Activated")], sprintf("%s target", tt$semantics$short[c("Inhibited", "Activated")])),
    if (members == "all") c(`member not called` = .pal$marks[["Not_significant"]])), shapes = "circle", size = 6)
  key <- .pai_badge("Rings", sprintf("inside out: \u2212log10 adj. P (colour: %s%s)  \u00b7  Inhibited | Activated share  \u00b7  members sorted by activity score (line: 0; dashed: \u00b1%s)",
    zc$name, if (cutoff < 1) sprintf("; dashed circle: adj. P = %s", cutoff) else "", thr))
  out <- .stack(list(.pai_header(title %||% "Top terms: significance, direction and member targets",
      sprintf("%s; %s members of each term on the outer ring, lowest to highest activity score", .sel_note(sel, n, terms, max_overlap, direction),
        if (members == "all") "all scored" else "called")),
    chips, patchwork::wrap_elements(full = body), key), c(.header_mm(1), 5, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 5, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- list(terms = sel, genes = genes[, c("tid", "gene", "score", "cls", "rank")],
    labelled = lab[, intersect(c("tid", "gene", "score"), names(lab))])
  .with_size(out, width, round(max(118, 50 + k * 5.5) + (if (label_genes > 0) 14 else 0), 1))   # gene names need a wider margin
}

# GOChord with an activity ring ----------------------------------------------------------------
#' Called targets shared by the top terms (GOChord)
#'
#' @inheritParams pai_plot_params
#' @param n_genes Number of targets, those in most terms first.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_gochord <- function(e, n = 5, n_genes = 24, terms = NULL, direction = c("both", "Activated", "Inhibited"), z = c("balanced", "goplot"),
                        cutoff = 0.05, max_overlap = 0.5, collections = NULL, width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction); z <- match.arg(z); zc <- .zcol(z)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  memb <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i)
    data.frame(gene = strsplit(sel$geneID[i], "/", fixed = TRUE)[[1]], tid = sel$tid[i], stringsAsFactors = FALSE)))
  cnt <- table(memb$gene)
  gdf <- data.frame(gene = names(cnt), n = as.integer(cnt), score = unname(tt$scores[names(cnt)]))
  gdf <- gdf[order(-gdf$n, -abs(gdf$score), gdf$gene), ]
  g <- utils::head(gdf$gene, n_genes)
  w <- memb[memb$gene %in% g, ]; tids <- sel$tid
  pal <- stats::setNames(c(.pai_module_colors, rep(.pai_other_module, max(0, length(tids) - length(.pai_module_colors))))[seq_along(tids)], tids)
  W <- nrow(w); gap <- 0.008 * 2 * pi; side_gap <- 0.08 * 2 * pi; span <- pi - side_gap
  gw <- table(factor(w$gene, levels = g)); tw <- table(factor(w$tid, levels = tids))
  ga <- (span - gap * (length(g) - 1)) * as.numeric(gw) / W
  ta <- (span - gap * (length(tids) - 1)) * as.numeric(tw) / W
  gstart <- stats::setNames(pi / 2 + side_gap / 2 + c(0, cumsum(ga + gap))[seq_along(g)], g)
  top <- pi / 2 - side_gap / 2
  tend <- stats::setNames(top - c(0, cumsum(ta + gap))[seq_along(tids)], tids); tstart <- tend - ta
  r <- 1; ring <- 0.06; gpos <- gstart; tcur <- tend; polys <- list(); kk <- 0
  for (gn in g) for (tid in tids) {
    if (!any(w$gene == gn & w$tid == tid)) next
    da <- (span - gap * (length(g) - 1)) / W; db <- (span - gap * (length(tids) - 1)) / W
    a1 <- gpos[[gn]]; a2 <- a1 + da; gpos[[gn]] <- a2
    b2 <- tcur[[tid]]; b1 <- b2 - db; tcur[[tid]] <- b1
    pp <- rbind(.arc(r, a1, a2), .bezier(c(r * cos(a2), r * sin(a2)), c(r * cos(b1), r * sin(b1))),
                .arc(r, b1, b2), .bezier(c(r * cos(b2), r * sin(b2)), c(r * cos(a1), r * sin(a1))))
    kk <- kk + 1; polys[[kk]] <- data.frame(x = pp[, 1], y = pp[, 2], id = kk, tid = tid, gene = gn)
  }
  rib <- do.call(rbind, polys); rib$cls <- unname(tt$classes[rib$gene])
  sect <- function(start, len, lab, r0, r1) do.call(rbind, lapply(seq_along(start), function(i) {
    a <- .arc(r0, start[i], start[i] + len[i]); b <- .arc(r1, start[i] + len[i], start[i])
    data.frame(x = c(a[, 1], b[, 1]), y = c(a[, 2], b[, 2]), id = lab[i], name = lab[i])
  }))
  gs <- sect(gstart, ga, g, r + 0.015, r + 0.015 + ring); gs$score <- tt$scores[gs$name]
  ts <- sect(tstart, ta, tids, r + 0.015, r + 0.015 + ring)
  glab <- data.frame(gene = g, a = gstart + ga / 2)
  glab$x <- (r + ring + 0.05) * cos(glab$a); glab$y <- (r + ring + 0.05) * sin(glab$a); glab$deg <- glab$a * 180 / pi - 180
  tlab <- data.frame(tid = tids, a = tstart + ta / 2)
  tlab$x <- (r + ring + 0.06) * cos(tlab$a); tlab$y <- (r + ring + 0.06) * sin(tlab$a)
  lim <- max(abs(tt$scores[g]))
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = rib, ggplot2::aes(x, y, group = id, fill = tid), alpha = 0.5, colour = "white", linewidth = 0.12) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_polygon(data = gs, ggplot2::aes(x, y, group = id, fill = score), colour = "white", linewidth = 0.3) +
    .div_scale(lim, "Activity score", guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(20, "mm"),
      barheight = ggplot2::unit(2, "mm"), title.position = "top")) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_polygon(data = ts, ggplot2::aes(x, y, group = id, fill = name), colour = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::geom_text(data = glab, ggplot2::aes(x, y, label = gene, angle = deg), hjust = 1, size = 2.1, colour = .tok$ink, family = .pai_font()) +
    ggplot2::geom_text(data = tlab, ggplot2::aes(x, y, label = tid), hjust = 0, size = 2.5, fontface = "bold", colour = .tok$ink, family = .pai_font()) +
    ggplot2::coord_equal(xlim = c(-1.6, 1.35), ylim = c(-1.25, 1.25), clip = "off") +
    ggplot2::theme_void(base_family = .pai_font()) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(size = 6.5, colour = .tok$ink2),
      legend.text = ggplot2::element_text(size = 6, colour = .tok$ink2))
  key <- .term_key(sel, width / 2.2 - 4, swatch = pal)
  body <- (p | key) + patchwork::plot_layout(widths = c(1.2, 1))
  out <- .stack(list(.pai_header(title %||% "Called targets shared by the top terms",
      sprintf("Left: the %d targets in most of the %d terms, ringed by activity score. Right: terms (%s)", length(g), nrow(sel),
        tolower(.sel_note(sel, n, terms, max_overlap, direction)))),
    patchwork::wrap_elements(full = body), .go_badge(tt, z)), c(.header_mm(1), NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 1, 6.5), c("mm", "null", "mm")))
  attr(out, "plot_data") <- list(links = unique(rib[, c("gene", "tid")]), terms = sel, genes = g)
  .with_size(out, width, 110)
}

# GOHeat: targets x terms, with a direction bar on top of each term ----------------------------------
# order = "count": targets by the number of shown terms they belong to, then |activity score|;
# "cluster": targets and terms by average-linkage clustering of the membership matrix (binary,
# i.e. Jaccard, distance). The bar above each term is its share of Inhibited (cool) and Activated
# (warm) called targets, over all its called targets.
#' Activity of the targets in the top terms (GOHeat)
#'
#' Targets by terms, coloured by activity score where the target is a member; above each term a
#' bar of its Inhibited and Activated shares.
#'
#' @inheritParams pai_plot_params
#' @param n_genes Number of targets, those in most terms first.
#' @param order `"count"` orders the targets by the number of terms they belong to; `"cluster"`
#'   clusters targets and terms by their membership (average linkage, binary distance).
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_goheat <- function(e, n = 10, n_genes = 30, terms = NULL, direction = c("both", "Activated", "Inhibited"), order = c("count", "cluster"),
                       z = c("balanced", "goplot"), cutoff = 0.05, max_overlap = 0.5, collections = NULL, width = 150, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction); order <- match.arg(order); z <- match.arg(z); zc <- .zcol(z)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  memb <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i)
    data.frame(gene = strsplit(sel$geneID[i], "/", fixed = TRUE)[[1]], tid = sel$tid[i], stringsAsFactors = FALSE)))
  cnt <- table(memb$gene)
  gdf <- data.frame(gene = names(cnt), n = as.integer(cnt), score = unname(tt$scores[names(cnt)]))
  gdf <- gdf[base::order(-gdf$n, -abs(gdf$score), gdf$gene), ]
  g <- utils::head(gdf, n_genes)
  M <- matrix(0, nrow(g), nrow(sel), dimnames = list(g$gene, sel$tid))
  M[cbind(match(memb$gene[memb$gene %in% g$gene], g$gene), match(memb$tid[memb$gene %in% g$gene], sel$tid))] <- 1
  gene_order <- g$gene; term_order <- sel$tid; hc_g <- hc_t <- NULL
  if (order == "cluster") {
    dg <- stats::dist(M, method = "binary"); dg[is.na(dg)] <- 1
    dt <- stats::dist(t(M), method = "binary"); dt[is.na(dt)] <- 1
    hc_g <- stats::hclust(dg, "average"); gene_order <- g$gene[hc_g$order]
    if (nrow(sel) > 2) { hc_t <- stats::hclust(dt, "average"); term_order <- sel$tid[hc_t$order] }
  }
  d <- expand.grid(gene = g$gene, tid = sel$tid, stringsAsFactors = FALSE)
  d$member <- M[cbind(match(d$gene, g$gene), match(d$tid, sel$tid))] == 1
  d$score <- ifelse(d$member, tt$scores[d$gene], NA)
  d$gene <- factor(d$gene, levels = rev(gene_order)); d$tid <- factor(d$tid, levels = term_order)
  lim <- max(abs(g$score))
  sel$zval <- sel[[zc$col]]; sel$f_inh <- sel$n_inhibited / (sel$n_inhibited + sel$n_activated)
  collab <- stats::setNames(sprintf("%s  %s", sel$tid, .short_label(sel$Description, 26)), sel$tid)
  xs <- ggplot2::scale_x_discrete(limits = term_order, expand = c(0, 0.5))
  tb <- sel; tb$x <- match(tb$tid, term_order)
  topbar <- ggplot2::ggplot(tb) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.4, xmax = x + 0.4, ymin = 0, ymax = f_inh), fill = .pal$marks[["Inhibited"]]) +
    ggplot2::geom_rect(ggplot2::aes(xmin = x - 0.4, xmax = x + 0.4, ymin = f_inh, ymax = 1), fill = .pal$marks[["Activated"]]) +
    ggplot2::geom_text(ggplot2::aes(x = x, y = 1.12, label = sprintf("%s\n%s", tid, .fmt_num(zval, 1))), vjust = 0, size = 1.95,
      lineheight = 0.95, colour = .tok$ink, family = .pai_font()) +
    ggplot2::scale_x_continuous(limits = c(0.5, nrow(sel) + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, 1.75), expand = c(0, 0)) + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 4, 0, 4))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = tid, y = gene)) +
    ggplot2::geom_tile(fill = .tok$pill, colour = "white", linewidth = 0.6) +
    ggplot2::geom_tile(data = d[d$member, ], ggplot2::aes(fill = score), colour = "white", linewidth = 0.6) +
    .div_scale(lim, "Activity score", guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(20, "mm"),
      barheight = ggplot2::unit(2, "mm"), title.position = "left")) +
    ggplot2::scale_x_discrete(limits = term_order, labels = collab) +
    ggplot2::labs(x = NULL, y = NULL) + theme_pai(grid = "none", axes = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1, colour = .tok$ink, size = ggplot2::rel(0.85)),
      axis.ticks = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(colour = .tok$ink, size = ggplot2::rel(0.88)),
      legend.position = "bottom", legend.justification = "left", plot.margin = ggplot2::margin(0, 4, 2, 4))
  body <- (topbar / p) + patchwork::plot_layout(heights = grid::unit(c(11, 1), c("mm", "null")))
  chips <- .pai_chips(c(`share Inhibited` = .pal$marks[["Inhibited"]], `share Activated` = .pal$marks[["Activated"]]), size = 6)
  out <- .stack(list(.pai_header(title %||% "Activity of the targets in the top terms",
      sprintf("Top %d called targets in %d terms%s; grey: not in the term. Above each term: Inhibited/Activated share and %s",
        nrow(g), nrow(sel), if (order == "cluster") ", clustered" else "", zc$name)),
    chips, body, .go_badge(tt, z)), c(.header_mm(1), 4.6, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 4.6, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- d
  attr(out, "order") <- list(genes = gene_order, terms = term_order, heights_genes = if (!is.null(hc_g)) sort(hc_g$height),
    merges_genes = if (!is.null(hc_g)) data.frame(a = hc_g$merge[, 1], b = hc_g$merge[, 2], height = hc_g$height), input = g$gene,
    heights_terms = if (!is.null(hc_t)) sort(hc_t$height), bars = sel[, c("tid", "n_inhibited", "n_activated", "f_inh")])
  .with_size(out, width, round(40 + 11 + nrow(g) * 3.9 + 32, 1))
}

# GOVenn: overlap of two or three terms, split by direction ------------------------------------
# Each region gives the number of called targets, how many are Inhibited and Activated, and a
# small bar of that split. Terms: the top n (2 or 3) non-redundant terms, or `terms`.
#' Overlap of two or three terms, split by direction (GOVenn)
#'
#' @inheritParams pai_plot_params
#' @param n 2 or 3 terms.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_govenn <- function(e, terms = NULL, n = 3, direction = c("both", "Activated", "Inhibited"), cutoff = 0.05, max_overlap = 0.5,
                       collections = NULL, width = 150, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  k <- nrow(sel); if (!k %in% 2:3) stop("pai_govenn() draws 2 or 3 terms.", call. = FALSE)
  sets <- lapply(sel$geneID, function(s) strsplit(s, "/", fixed = TRUE)[[1]]); names(sets) <- sel$tid
  allg <- sort(unique(unlist(sets)))
  code <- vapply(allg, function(gn) paste(as.integer(vapply(sets, function(s) gn %in% s, logical(1))), collapse = ""), character(1))
  if (k == 3) {
    ctr <- data.frame(tid = sel$tid, x0 = c(-0.55, 0.55, 0), y0 = c(0.32, 0.32, -0.63))
    pos <- data.frame(code = c("100", "010", "001", "110", "101", "011", "111"),
      x = c(-1.02, 1.02, 0, 0, -0.62, 0.62, 0), y = c(0.62, 0.62, -1.18, 0.80, -0.36, -0.36, 0.02))
    tl <- data.frame(tid = sel$tid, x = c(-1.55, 1.55, 0), y = c(1.32, 1.32, -1.78), h = c(0, 1, 0.5))
  } else {
    ctr <- data.frame(tid = sel$tid, x0 = c(-0.5, 0.5), y0 = c(0, 0))
    pos <- data.frame(code = c("10", "01", "11"), x = c(-1.0, 1.0, 0), y = c(0, 0, 0))
    tl <- data.frame(tid = sel$tid, x = c(-1.5, 1.5), y = c(1.12, 1.12), h = c(0, 1))
  }
  reg <- do.call(rbind, lapply(pos$code, function(cd) {
    gs <- allg[code == cd]
    data.frame(code = cd, n = length(gs), n_inh = sum(tt$classes[gs] == "Inhibited"), n_act = sum(tt$classes[gs] == "Activated"))
  }))
  reg <- merge(pos, reg, by = "code", sort = FALSE)
  reg$f_inh <- ifelse(reg$n > 0, reg$n_inh / reg$n, 0)
  pal <- stats::setNames(.pai_module_colors[seq_len(k)], sel$tid)
  bw <- 0.2; by <- 0.035
  p <- ggplot2::ggplot() +
    ggforce::geom_circle(data = ctr, ggplot2::aes(x0 = x0, y0 = y0, r = 1, colour = tid, fill = tid), alpha = 0.06, linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") + ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::geom_text(data = reg, ggplot2::aes(x = x, y = y + 0.13, label = n), size = 3.2, fontface = "bold", colour = .tok$ink, family = .pai_font()) +
    ggplot2::geom_text(data = reg, ggplot2::aes(x = x - 0.03, y = y - 0.02, label = sprintf("%d\u2193", n_inh)), hjust = 1, size = 2.1,
      colour = .pal$text[["Inhibited"]], fontface = "bold", family = .pai_font()) +
    ggplot2::geom_text(data = reg, ggplot2::aes(x = x + 0.03, y = y - 0.02, label = sprintf("%d\u2191", n_act)), hjust = 0, size = 2.1,
      colour = .pal$text[["Activated"]], fontface = "bold", family = .pai_font()) +
    ggplot2::geom_rect(data = reg[reg$n > 0, ], ggplot2::aes(xmin = x - bw, xmax = x - bw + 2 * bw * f_inh, ymin = y - 0.13 - by, ymax = y - 0.13 + by),
      fill = .pal$marks[["Inhibited"]]) +
    ggplot2::geom_rect(data = reg[reg$n > 0, ], ggplot2::aes(xmin = x - bw + 2 * bw * f_inh, xmax = x + bw, ymin = y - 0.13 - by, ymax = y - 0.13 + by),
      fill = .pal$marks[["Activated"]]) +
    ggplot2::geom_text(data = tl, ggplot2::aes(x = x, y = y, label = tid, hjust = h), size = 2.8, fontface = "bold", colour = .tok$ink, family = .pai_font()) +
    ggplot2::coord_equal(clip = "off") + ggplot2::theme_void()
  chips <- .pai_chips(stats::setNames(pal, sprintf("%s  %s  (%s)", sel$tid, .short_label(sel$Description, 40),
    .dir_counts(sel$n_inhibited, sel$n_activated))), per_row = 1, size = 6.1)
  out <- .stack(list(.pai_header(title %||% "Overlap of the top terms, by direction",
      sprintf("Called targets of %s; each region: number of targets, Inhibited \u2193 and Activated \u2191, and their split", tolower(.sel_note(sel, n, terms, max_overlap, direction)))),
    chips, p), c(.header_mm(1), 4.2 * k, NA)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 4.3 * k, 1), c("mm", "mm", "null")))
  attr(out, "plot_data") <- reg[, c("code", "n", "n_inh", "n_act")]
  attr(out, "terms") <- sel[, c("tid", "key")]
  .with_size(out, width, round(25 + 4.3 * k + width * 0.62, 1))
}

# GOCluster: circular dendrogram of the targets by term membership ------------------------------
.dendro_segments <- function(hc, leaf = 0) {
  m <- length(hc$order); pos <- numeric(m); pos[hc$order] <- seq_len(m)
  nx <- numeric(nrow(hc$merge)); nh <- hc$height; seg <- list()
  xy <- function(j) if (j < 0) c(pos[-j], leaf) else c(nx[j], nh[j])
  for (i in seq_len(nrow(hc$merge))) {
    a <- xy(hc$merge[i, 1]); b <- xy(hc$merge[i, 2]); nx[i] <- (a[1] + b[1]) / 2
    seg[[length(seg) + 1]] <- data.frame(x = c(a[1], b[1], a[1]), xend = c(a[1], b[1], b[1]), y = c(a[2], b[2], nh[i]), yend = c(nh[i], nh[i], nh[i]))
  }
  do.call(rbind, seg)
}

# Targets of the top terms clustered by which terms they belong to (average linkage on binary,
# i.e. Jaccard, distance), drawn as a circular dendrogram. Rings outside the leaves: activity
# score, then one ring per term (filled when the target is a member).
#' Targets of the top terms clustered by membership (GOCluster)
#'
#' A circular dendrogram of the called targets by the terms they belong to; rings outside the
#' leaves give the activity score and the membership of each term.
#'
#' @inheritParams pai_plot_params
#' @param n_genes Number of targets, those in most terms first.
#' @param method Linkage passed to [stats::hclust()].
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_gocluster <- function(e, n = 5, n_genes = 40, terms = NULL, direction = c("both", "Activated", "Inhibited"), method = "average",
                          cutoff = 0.05, max_overlap = 0.5, collections = NULL, width = 180, title = NULL) {
  tt <- .as_tt(e)
  direction <- match.arg(direction)
  sel <- .select_terms(tt, terms, n, cutoff, max_overlap, collections, direction)
  k <- nrow(sel)
  memb <- do.call(rbind, lapply(seq_len(k), function(i)
    data.frame(gene = strsplit(sel$geneID[i], "/", fixed = TRUE)[[1]], tid = sel$tid[i], stringsAsFactors = FALSE)))
  cnt <- table(memb$gene)
  gdf <- data.frame(gene = names(cnt), n = as.integer(cnt), score = unname(tt$scores[names(cnt)]), stringsAsFactors = FALSE)
  gdf <- utils::head(gdf[base::order(-gdf$n, -abs(gdf$score), gdf$gene), ], n_genes)
  M <- matrix(0, nrow(gdf), k, dimnames = list(gdf$gene, sel$tid))
  mm <- memb[memb$gene %in% gdf$gene, ]; M[cbind(match(mm$gene, gdf$gene), match(mm$tid, sel$tid))] <- 1
  dg <- stats::dist(M, method = "binary"); dg[is.na(dg)] <- 1
  hc <- stats::hclust(dg, method)
  m <- nrow(gdf); H <- max(hc$height, 1e-9)
  seg <- .dendro_segments(hc, leaf = -1)
  rr <- function(h) ifelse(h < 0, 1, 0.12 + 0.80 * (1 - h / H))                           # leaves at radius 1, root near the centre
  seg$r <- rr(seg$y); seg$rend <- rr(seg$yend)
  leaf <- data.frame(gene = gdf$gene[hc$order], pos = seq_len(m))
  leaf$score <- gdf$score[match(leaf$gene, gdf$gene)]; leaf$cls <- unname(tt$classes[leaf$gene])
  r_act <- c(1.04, 1.14); ring_w <- 0.075
  rings <- do.call(rbind, lapply(seq_len(k), function(j) data.frame(gene = leaf$gene, pos = leaf$pos, tid = sel$tid[j],
    member = M[leaf$gene, j] == 1, r0 = r_act[2] + 0.02 + (j - 1) * (ring_w + 0.012), r1 = r_act[2] + 0.02 + (j - 1) * (ring_w + 0.012) + ring_w)))
  r_lab <- max(rings$r1) + 0.05
  leaf$ang <- 90 - 360 * (leaf$pos - 0.5) / m
  leaf$flip <- leaf$ang < -90
  pal <- stats::setNames(c(.pai_module_colors, rep(.pai_other_module, max(0, k - length(.pai_module_colors))))[seq_len(k)], sel$tid)
  lim <- max(abs(gdf$score))
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = x, xend = xend, y = r, yend = rend), colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::geom_rect(data = leaf, ggplot2::aes(xmin = pos - 0.5, xmax = pos + 0.5, ymin = r_act[1], ymax = r_act[2], fill = score),
      colour = "white", linewidth = 0.2) +
    .div_scale(lim, "Activity score", guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(22, "mm"), barheight = ggplot2::unit(2, "mm"),
      title.position = "left")) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_rect(data = rings, ggplot2::aes(xmin = pos - 0.5, xmax = pos + 0.5, ymin = r0, ymax = r1,
      fill = ifelse(member, tid, "none")), colour = "white", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c(pal, none = .tok$pill), guide = "none") +
    ggplot2::geom_text(data = leaf, ggplot2::aes(x = pos, y = r_lab, label = gene, angle = ifelse(flip, ang + 180, ang),
      hjust = ifelse(flip, 1, 0)), size = 2.0, colour = .tok$ink, family = .pai_font()) +
    ggplot2::scale_x_continuous(limits = c(0.5, m + 0.5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, r_lab + 0.42), expand = c(0, 0)) +
    ggplot2::coord_polar(theta = "x", clip = "off") + ggplot2::theme_void(base_family = .pai_font()) +
    ggplot2::theme(legend.position = "bottom", legend.title = ggplot2::element_text(size = 6.4, colour = .tok$ink2),
      legend.text = ggplot2::element_text(size = 6, colour = .tok$ink2))
  key <- .term_key(sel, width / 2.35 - 4, swatch = pal)
  body <- (p | (patchwork::plot_spacer() / key / patchwork::plot_spacer() + patchwork::plot_layout(heights = c(1, 0.25 * k, 1)))) +
    patchwork::plot_layout(widths = c(1.35, 1))
  out <- .stack(list(.pai_header(title %||% "Targets of the top terms, clustered by membership",
      sprintf("%d called targets in most of the %d terms; rings from inside: activity score, then T1\u2013T%d membership", m, k, k)),
    patchwork::wrap_elements(full = body), .pai_badge("Clustering", sprintf("%s linkage on binary (Jaccard) distance of term membership", method))),
    c(.header_mm(1), NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 1, 6.5), c("mm", "null", "mm")))
  attr(out, "plot_data") <- list(genes = leaf[, c("gene", "pos", "score", "cls")], membership = rings[rings$member, c("gene", "tid")],
    heights = sort(hc$height), terms = sel[, c("tid", "key")],
    merges = data.frame(a = hc$merge[, 1], b = hc$merge[, 2], height = hc$height), input = gdf$gene)
  .with_size(out, width, 130)
}
