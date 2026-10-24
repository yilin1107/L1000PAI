# Enrichment hub analysis.
#
# Significant terms of one direction are pooled across collections and grouped into modules by
# the overlap of their genes (Jaccard index, average-linkage clustering, cut at a fixed
# similarity). Redundant terms (a GO parent and child, the same pathway in KEGG and Reactome)
# then count once, and a hub gene is a target that sits in many modules, not merely in many
# overlapping terms.

# .pai_module_colors and .pai_other_module come from plot_style.R.

# Hub analysis of an ORA (pai_enrich) or GSEA (pai_gsea_pathways) result
# @param res Result of pai_enrich() or pai_gsea_pathways().
# @param x Optional L1000PAIResult, to attach activity scores to the genes.
# @param direction "Activated" or "Inhibited" (for GSEA: the end of the ranking).
# @param similarity Jaccard similarity at which terms join the same module.
#' Hub genes of the enriched pathways
#'
#' Significant terms of one direction are pooled across collections and grouped into modules by
#' the overlap of their genes (Jaccard index, average-linkage clustering cut at `similarity`), so
#' that redundant terms count once. A hub gene is a target in many modules, not merely in many
#' overlapping terms. Targets annotated in more tested terms than 95% of all targets are flagged as
#' broad.
#'
#' @param res Result of [pai_enrich()] or [pai_gsea_pathways()].
#' @param direction `"Inhibited"` or `"Activated"` (for GSEA: the end of the ranking).
#' @param cutoff Adjusted P cut-off for the terms.
#' @param collections Collections to pool; `NULL` uses all.
#' @param similarity Jaccard similarity at which terms join the same module.
#' @return A `pai_hubs` object with `terms`, `modules`, `genes` (hub ranking), `membership` and
#'   `jaccard`.
#' @family hub plots
#' @export
pai_enrich_hubs <- function(res, direction = c("Inhibited", "Activated"), cutoff = 0.05,
                            collections = NULL, similarity = 0.25) {
  direction <- match.arg(direction)
  if (inherits(res, "pai_enrichment")) {
    t <- res$table[res$table$direction == direction & res$table$p.adjust < cutoff, ]
    genes_col <- "geneID"; source <- "ORA"
  } else if (inherits(res, "pai_pathway_gsea")) {
    s <- if (direction == "Activated") 1 else -1
    t <- res$table[sign(res$table$NES) == s & res$table$p.adjust < cutoff, ]
    genes_col <- "core_enrichment"; source <- "GSEA leading edge"
  } else stop("`res` must come from pai_enrich() or pai_gsea_pathways().")
  if (!is.null(collections)) t <- t[t$collection %in% collections, ]
  if (nrow(t) < 2) stop(sprintf(paste0("Fewer than two terms with adj. P < %s in this direction; nothing to group. ",
    "A larger cutoff (e.g. cutoff = 1) groups the top terms descriptively."), cutoff), call. = FALSE)
  t$key <- paste(t$collection, t$ID, sep = "|")
  gl <- stats::setNames(strsplit(t[[genes_col]], "/", fixed = TRUE), t$key)
  keys <- names(gl); n <- length(keys)
  inter <- matrix(0, n, n, dimnames = list(keys, keys))
  allg <- sort(unique(unlist(gl)))
  inc <- vapply(gl, function(g) allg %in% g, logical(length(allg)))          # genes x terms
  inter <- crossprod(inc * 1)
  sz <- colSums(inc)
  jac <- inter / (outer(sz, sz, "+") - inter)
  hc <- stats::hclust(stats::as.dist(1 - jac), method = "average")
  cl <- stats::cutree(hc, h = 1 - similarity)
  t$cluster <- cl[t$key]
  # Modules ranked by their best adjusted P (ties: more terms, then name); label = best term.
  mods <- do.call(rbind, lapply(split(t, t$cluster), function(z) {
    z <- z[order(z$p.adjust, z$Description), ]
    data.frame(cluster = z$cluster[1], label = z$Description[1], label_source = z$collection[1],
      best_padj = z$p.adjust[1], n_terms = nrow(z), n_genes = length(unique(unlist(gl[z$key]))),
      collections = paste(sort(unique(z$collection)), collapse = ", "), stringsAsFactors = FALSE)
  }))
  mods <- mods[order(mods$best_padj, -mods$n_terms, mods$label), ]
  mods$module <- sprintf("M%d", seq_len(nrow(mods)))
  t$module <- mods$module[match(t$cluster, mods$cluster)]
  memb <- do.call(rbind, lapply(seq_len(nrow(t)), function(i)
    data.frame(gene = gl[[t$key[i]]], key = t$key[i], module = t$module[i], stringsAsFactors = FALSE)))
  genes <- do.call(rbind, lapply(split(memb, memb$gene), function(z) data.frame(gene = z$gene[1],
    n_terms = nrow(z), n_modules = length(unique(z$module)),
    modules = paste(mods$module[mods$module %in% z$module], collapse = ";"), stringsAsFactors = FALSE)))
  sc <- if (inherits(res, "pai_enrichment")) stats::setNames(res$targets$activity_score, res$targets$gene) else res$stats
  genes$activity_score <- unname(sc[genes$gene])
  # Annotation breadth: tested terms (all collections used, background genes, size limits) containing the gene.
  bg <- if (inherits(res, "pai_enrichment")) res$background else names(res$stats)
  cols_used <- unique(t$collection); br <- stats::setNames(numeric(length(bg)), bg)
  for (cl in cols_used) {
    gs <- .read_gmt2(res$sets[[cl]]$file)
    sets <- split(gs$TERM2GENE$gene, gs$TERM2GENE$term)
    sets <- lapply(sets, function(z) intersect(unique(z), bg))
    sets <- sets[lengths(sets) >= res$params$min_size & lengths(sets) <= res$params$max_size]
    tab <- table(unlist(sets)); br[names(tab)] <- br[names(tab)] + as.numeric(tab)
  }
  thr <- stats::quantile(br[br > 0], 0.95, names = FALSE, type = 7)
  genes$breadth <- as.numeric(br[genes$gene]); genes$broad <- genes$breadth >= thr
  genes <- genes[order(-genes$n_modules, -genes$n_terms, -abs(genes$activity_score), genes$gene), ]
  genes$hub_rank <- seq_len(nrow(genes)); rownames(genes) <- NULL
  structure(list(terms = t, modules = mods, genes = genes, membership = memb, jaccard = jac,
    direction = direction, source = source, library = res$library, semantics = res$semantics,
    params = list(cutoff = cutoff, similarity = similarity, broad_threshold = thr)), class = "pai_hubs")
}

# Five module hues pass the colour-vision check for all pairs; a sixth shown module is drawn grey
# (it keeps its M# label), as are the remaining modules ("Other").
.module_palette <- function(mods, n_modules) {
  shown <- utils::head(mods$module, n_modules)
  cols <- c(.pai_module_colors, rep(.pai_other_module, max(0, length(shown) - length(.pai_module_colors))))[seq_along(shown)]
  stats::setNames(c(cols, .pai_other_module), c(shown, "Other"))
}
.short_label <- function(s, n = 34) ifelse(nchar(s) > n, paste0(substr(s, 1, n - 1), "\u2026"), s)
.hub_title <- function(h, what) sprintf("%s of the %s targets", what, tolower(h$semantics$short[[h$direction]]))
.hub_badge <- function(h) .pai_badge("Modules", sprintf("%s terms, adj. P < %s, all collections  \u00b7  join at Jaccard \u2265 %s",
  h$source, h$params$cutoff, h$params$similarity))

# Membership matrix: hub genes x modules --------------------------------------------------------
#' Hub genes by module membership
#'
#' @param h Result of [pai_enrich_hubs()].
#' @param n_genes Number of hub genes.
#' @param n_modules Number of modules shown as columns (five have their own colour).
#' @inheritParams pai_plot_params
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family hub plots
#' @export
pai_hubplot <- function(h, n_genes = 15, n_modules = 6, title = NULL) {
  pal <- .module_palette(h$modules, n_modules); mods <- names(pal)[names(pal) != "Other"]
  g <- utils::head(h$genes, n_genes)
  cnt <- stats::aggregate(key ~ gene + module, data = h$membership[h$membership$gene %in% g$gene, ], FUN = length)
  names(cnt)[3] <- "n_terms"
  cnt <- cnt[cnt$module %in% mods, ]
  grid_df <- expand.grid(gene = g$gene, module = mods, stringsAsFactors = FALSE)
  d <- merge(grid_df, cnt, all.x = TRUE); d$in_module <- !is.na(d$n_terms)
  d$gene <- factor(d$gene, levels = rev(g$gene)); d$col <- match(d$module, mods)
  bars <- data.frame(gene = factor(g$gene, levels = rev(g$gene)), n_modules = g$n_modules, n_terms = g$n_terms, broad = g$broad)
  k <- length(mods); bx <- k + 0.9; bw <- 2.2
  mlab <- h$modules[match(mods, h$modules$module), ]
  head_df <- data.frame(col = seq_len(k), lab = sprintf("%s  %s", mods, .short_label(mlab$label, 30)),
    n = mlab$n_terms, module = mods)
  maxm <- max(h$genes$n_modules)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = col, y = gene)) +
    ggplot2::geom_point(data = d[!d$in_module, ], colour = .tok$grid, size = 1.1) +
    ggplot2::geom_point(data = d[d$in_module, ], ggplot2::aes(fill = module, size = n_terms), shape = 21, colour = "white", stroke = 0.5) +
    ggplot2::geom_segment(data = bars, inherit.aes = FALSE, ggplot2::aes(x = bx, xend = bx + bw, y = gene, yend = gene),
      colour = .tok$pill, linewidth = 2.2) +
    ggplot2::geom_segment(data = bars, inherit.aes = FALSE, ggplot2::aes(x = bx, xend = bx + bw * n_modules / maxm, y = gene, yend = gene),
      colour = .tok$ink2, linewidth = 2.2) +
    ggplot2::geom_text(data = bars, inherit.aes = FALSE, ggplot2::aes(x = bx + bw + 0.15, y = gene,
      label = sprintf("%d modules \u00b7 %d terms", n_modules, n_terms)), hjust = 0, size = 2.05, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::scale_y_discrete(labels = function(z) ifelse(z %in% g$gene[g$broad], paste0(z, " \u2020"), z)) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.6, 4.2), name = "Terms of the module containing the gene",
      breaks = scales::breaks_pretty(3), guide = ggplot2::guide_legend(title.position = "left", nrow = 1,
        override.aes = list(fill = .tok$muted))) +
    ggplot2::scale_x_continuous(breaks = seq_len(k), labels = head_df$lab, limits = c(0.5, bx + bw + 3.2), expand = c(0, 0),
      position = "top") +
    ggplot2::labs(x = NULL, y = NULL) + theme_pai(grid = "none") +
    ggplot2::theme(axis.text.x.top = ggplot2::element_text(angle = 40, hjust = 0, vjust = 0, colour = .tok$ink, size = ggplot2::rel(0.85)),
      axis.line.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(colour = .tok$ink, face = "bold"),
      legend.position = "bottom", legend.justification = "left")
  # module colour chips above the column labels
  p <- p + ggplot2::geom_point(data = head_df, inherit.aes = FALSE, ggplot2::aes(x = col, y = length(g$gene) + 0.9, fill = module),
    shape = 22, size = 2.6, colour = "white", stroke = 0.3) + ggplot2::coord_cartesian(clip = "off")
  out <- patchwork::wrap_plots(list(.pai_header(title %||% .hub_title(h, "Hub genes"),
      sprintf("Targets ranked by the modules they belong to (%d modules, %d terms)\n\u2020 broad: in more tested terms than 95%% of targets",
        nrow(h$modules), nrow(h$terms))),
    p, .hub_badge(h)), ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(2) + 1, 1, 6.5), c("mm", "null", "mm")))
  attr(out, "plot_data") <- d
  .with_size(out, 92, round(60 + length(g$gene) * 3.4, 1))
}

# Chord diagram: hub genes (left) to modules (right) ----------------------------------------------
.arc <- function(r, a1, a2, n = 40) { a <- seq(a1, a2, length.out = n); cbind(r * cos(a), r * sin(a)) }
.bezier <- function(p0, p2, n = 40) {               # quadratic, control point at the centre
  s <- seq(0, 1, length.out = n)
  cbind((1 - s)^2 * p0[1] + s^2 * p2[1], (1 - s)^2 * p0[2] + s^2 * p2[2])
}

#' Hub genes linked to their modules (chord diagram)
#'
#' @inheritParams pai_hubplot
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family hub plots
#' @export
pai_chordplot <- function(h, n_genes = 12, n_modules = 6, title = NULL) {
  pal <- .module_palette(h$modules, n_modules); mods <- names(pal)[names(pal) != "Other"]
  g <- utils::head(h$genes[h$genes$n_modules > 0, ], n_genes)$gene
  w <- unique(h$membership[h$membership$gene %in% g & h$membership$module %in% mods, c("gene", "module")])
  w$weight <- 1
  g <- g[g %in% w$gene]
  W <- sum(w$weight); gap <- 0.012 * 2 * pi; side_gap <- 0.08 * 2 * pi
  span <- pi - side_gap                                     # each side
  # genes: left side, top to bottom (angles decreasing from pi/2 + side_gap/2 ... going CCW)
  gw <- tapply(w$weight, w$gene, sum)[g]; mw <- tapply(w$weight, w$module, sum)[mods]; mw <- mw[!is.na(mw)]
  mods <- names(mw)
  ga <- (span - gap * (length(g) - 1)) * gw / W
  ma <- (span - gap * (length(mods) - 1)) * mw / W
  gstart <- pi / 2 + side_gap / 2 + c(0, cumsum(ga + gap))[seq_along(g)]
  top <- pi / 2 - side_gap / 2
  mend <- top - c(0, cumsum(ma + gap))[seq_along(mods)]
  mstart <- mend - ma
  names(gstart) <- g; names(mstart) <- mods; names(mend) <- mods
  r <- 1; ring <- 0.06
  gpos <- gstart; mcur <- mend
  polys <- list(); k <- 0
  for (gn in g) for (m in mods) {
    ww <- w$weight[w$gene == gn & w$module == m]; if (!length(ww)) next
    da <- (span - gap * (length(g) - 1)) * ww / W
    db <- (span - gap * (length(mods) - 1)) * ww / W
    a1 <- gpos[[gn]]; a2 <- a1 + da; gpos[[gn]] <- a2
    b2 <- mcur[[m]]; b1 <- b2 - db; mcur[[m]] <- b1
    p <- rbind(.arc(r, a1, a2), .bezier(c(r * cos(a2), r * sin(a2)), c(r * cos(b1), r * sin(b1))),
               .arc(r, b1, b2), .bezier(c(r * cos(b2), r * sin(b2)), c(r * cos(a1), r * sin(a1))))
    k <- k + 1; polys[[k]] <- data.frame(x = p[, 1], y = p[, 2], id = k, module = m, gene = gn, weight = ww)
  }
  rib <- do.call(rbind, polys)
  sect <- function(start, len, lab, kind) do.call(rbind, lapply(seq_along(start), function(i) {
    a <- .arc(r + 0.015, start[i], start[i] + len[i]); b <- .arc(r + 0.015 + ring, start[i] + len[i], start[i])
    data.frame(x = c(a[, 1], b[, 1]), y = c(a[, 2], b[, 2]), id = paste(kind, i), name = lab[i], kind = kind)
  }))
  gs <- sect(gstart, ga, g, "gene"); ms <- sect(mstart, ma, mods, "module")
  glab <- data.frame(gene = g, a = gstart + ga / 2); glab$x <- (r + ring + 0.05) * cos(glab$a); glab$y <- (r + ring + 0.05) * sin(glab$a)
  glab$deg <- glab$a * 180 / pi - 180
  hub <- h$genes[match(g, h$genes$gene), ]
  glab$lab <- sprintf("%s  %d", g, hub$n_modules)
  mlab <- data.frame(module = mods, a = mstart + ma / 2)
  mlab$x0 <- (r + ring + 0.02) * cos(mlab$a); mlab$y0 <- (r + ring + 0.02) * sin(mlab$a)
  ml <- h$modules[match(mods, h$modules$module), ]
  mlab$lab <- sprintf("%s  %s", mods, vapply(ml$label, function(z) paste(strwrap(.short_label(z, 30), 19), collapse = "\n"), character(1)))
  nl <- lengths(strsplit(mlab$lab, "\n"))
  # spread labels vertically (top to bottom) so wrapped labels never overlap; leader lines join them to their arcs
  need <- 0.13 * nl + 0.04
  yy <- mlab$y0
  for (i in seq_along(yy)[-1]) yy[i] <- min(yy[i], yy[i - 1] - (need[i - 1] + need[i]) / 2)
  if (min(yy) < -1.2) yy <- yy + (-1.2 - min(yy))
  for (i in rev(seq_along(yy))[-1]) yy[i] <- max(yy[i], yy[i + 1] + (need[i] + need[i + 1]) / 2)
  mlab$y <- yy; mlab$x <- pmax(mlab$x0, 0) + 0.28
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = rib, ggplot2::aes(x, y, group = id, fill = module), alpha = 0.5, colour = "white", linewidth = 0.15) +
    ggplot2::geom_polygon(data = gs, ggplot2::aes(x, y, group = id), fill = .tok$ink2, colour = "white", linewidth = 0.3) +
    ggplot2::geom_polygon(data = ms, ggplot2::aes(x, y, group = id, fill = name), colour = "white", linewidth = 0.3) +
    ggplot2::geom_text(data = glab, ggplot2::aes(x, y, label = lab, angle = deg), hjust = 1, size = 2.3, colour = .tok$ink,
      family = .pai_font()) +
    ggplot2::geom_segment(data = mlab, ggplot2::aes(x = x0, y = y0, xend = x - 0.04, yend = y), colour = .tok$axis, linewidth = 0.3) +
    ggplot2::geom_text(data = mlab, ggplot2::aes(x, y, label = lab), hjust = 0, vjust = 0.5, size = 2.25, colour = .tok$ink,
      family = .pai_font(), lineheight = 0.9) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::coord_equal(xlim = c(-1.75, 2.45), ylim = c(-1.25, 1.25), clip = "off") +
    ggplot2::theme_void(base_family = .pai_font())
  out <- patchwork::wrap_plots(list(.pai_header(title %||% .hub_title(h, "Hub genes and modules"),
      "Hub genes (left, number = modules) linked to their modules (right)"),
    p, .hub_badge(h)), ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 1, 6.5), c("mm", "null", "mm")))
  attr(out, "plot_data") <- unique(rib[, c("gene", "module", "weight")])
  .with_size(out, 92, 110)
}

# Module-gene network ----------------------------------------------------------------------------
#' Module-gene network of the hub analysis
#'
#' @inheritParams pai_hubplot
#' @inheritParams pai_plot_params
#' @param n_genes Number of targets, those in most modules first.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family hub plots
#' @export
pai_enrich_cnetplot <- function(h, n_modules = 5, n_genes = 20, seed = 1, title = NULL) {
  pal <- .module_palette(h$modules, n_modules); mods <- names(pal)[names(pal) != "Other"]
  memb <- unique(h$membership[h$membership$module %in% mods, c("gene", "module")])
  hubs <- h$genes[h$genes$gene %in% memb$gene, ]
  g <- utils::head(hubs, n_genes)$gene
  e <- memb[memb$gene %in% g, ]
  ml <- h$modules[match(mods, h$modules$module), ]
  vert <- data.frame(name = c(mods, g), type = c(rep("module", length(mods)), rep("gene", length(g))))
  gr <- igraph::graph_from_data_frame(e[, c("module", "gene")], directed = FALSE, vertices = vert)
  set.seed(seed); lay <- igraph::layout_with_fr(gr, niter = 2000)
  nodes <- cbind(vert, x = lay[, 1], y = lay[, 2])
  nodes$n_modules <- h$genes$n_modules[match(nodes$name, h$genes$gene)]
  nodes$n_terms_mod <- ml$n_terms[match(nodes$name, mods)]
  nodes$lab <- nodes$name
  ee <- merge(merge(e, nodes[, c("name", "x", "y")], by.x = "module", by.y = "name"),
              nodes[, c("name", "x", "y")], by.x = "gene", by.y = "name", suffixes = c("", "end"))
  gn <- nodes[nodes$type == "gene", ]; mn <- nodes[nodes$type == "module", ]
  gn$hub <- gn$n_modules >= max(2, stats::quantile(gn$n_modules, 0.75))
  gn$score <- h$genes$activity_score[match(gn$name, h$genes$gene)]
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = ee, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = module), curvature = 0.1,
      linewidth = 0.45, alpha = 0.6) +
    ggplot2::geom_point(data = gn, ggplot2::aes(x, y, size = n_modules, fill = score), shape = 21, colour = "white", stroke = 0.6) +
    ggplot2::scale_fill_gradient2(low = .pal$marks[["Inhibited"]], mid = .pai_mid, high = .pal$marks[["Activated"]],
      limits = c(-1, 1) * max(abs(h$genes$activity_score), na.rm = TRUE), labels = .lab_minus, name = "Activity",
      guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(20, "mm"), barheight = ggplot2::unit(2, "mm"), title.position = "left")) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(data = mn, ggplot2::aes(x, y, fill = name), shape = 22, size = 7.5, colour = "white", stroke = 1) +
    ggplot2::geom_text(data = mn, ggplot2::aes(x, y, label = name), size = 2.3, fontface = "bold", colour = "white", family = .pai_font()) +
    do.call(ggrepel::geom_text_repel, c(list(data = gn, mapping = ggplot2::aes(x, y, label = lab, fontface = ifelse(hub, "bold", "plain")),
      size = 2.2, family = .pai_font(), colour = .tok$ink, segment.size = 0.2, segment.colour = .tok$axis, box.padding = 0.2,
      max.overlaps = Inf, seed = seed), .halo)) +
    ggplot2::scale_colour_manual(values = pal, guide = "none") + ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.4, 4), name = "Modules per gene", breaks = scales::breaks_pretty(3),
      guide = ggplot2::guide_legend(title.position = "left", nrow = 1, override.aes = list(shape = 16, colour = .tok$muted))) +
    ggplot2::coord_cartesian(clip = "off") + theme_pai(grid = "none") +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(), legend.position = "bottom", legend.justification = "left")
  chips <- .pai_chips(stats::setNames(pal[mods], sprintf("%s %s", mods, .short_label(ml$label, 30))), size = 6.2, per_row = 2)
  out <- patchwork::wrap_plots(list(.pai_header(title %||% .hub_title(h, "Module\u2013gene network"),
      sprintf("Squares: top %d modules  \u00b7  dots: %d targets in most modules  \u00b7  bold: hub", length(mods), length(g))),
    chips, p, .hub_badge(h)), ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 4.2 * ceiling(length(mods) / 2), 1, 6.5), c("mm", "mm", "null", "mm")))
  e$activity_score <- h$genes$activity_score[match(e$gene, h$genes$gene)]
  attr(out, "plot_data") <- e
  .with_size(out, 92, 110)
}

# Enrichment map: term similarity network --------------------------------------------------------
#' Enrichment map of the hub analysis
#'
#' The most significant terms as nodes, linked when their gene sets overlap by at least the
#' module similarity; colour shows the module.
#'
#' @inheritParams pai_hubplot
#' @inheritParams pai_plot_params
#' @param n Number of terms.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family hub plots
#' @export
pai_emapplot <- function(h, n = 40, n_modules = 6, seed = 1, title = NULL) {
  pal <- .module_palette(h$modules, n_modules)
  t <- h$terms[order(h$terms$p.adjust, h$terms$key), ]
  t <- utils::head(t, n)
  t$colour_group <- ifelse(t$module %in% names(pal), t$module, "Other")
  j <- h$jaccard[t$key, t$key, drop = FALSE]
  idx <- which(upper.tri(j) & j >= h$params$similarity, arr.ind = TRUE)
  ed <- data.frame(from = t$key[idx[, 1]], to = t$key[idx[, 2]], jaccard = j[idx])
  gr <- igraph::graph_from_data_frame(ed, directed = FALSE, vertices = data.frame(name = t$key))
  set.seed(seed); lay <- igraph::layout_with_fr(gr, weights = if (nrow(ed)) ed$jaccard else NULL, niter = 3000)
  nodes <- cbind(t, x = lay[, 1], y = lay[, 2])
  nodes$size <- if ("k" %in% names(nodes)) nodes$k else lengths(strsplit(nodes$core_enrichment, "/"))
  e2 <- merge(merge(ed, nodes[, c("key", "x", "y")], by.x = "from", by.y = "key"), nodes[, c("key", "x", "y")],
    by.x = "to", by.y = "key", suffixes = c("", "end"))
  rep_lab <- do.call(rbind, lapply(split(nodes, nodes$module), function(z) z[which.min(z$p.adjust), ]))
  rep_lab <- rep_lab[rep_lab$module %in% names(pal), ]
  rep_lab$lab <- sprintf("%s  %s", rep_lab$module, vapply(rep_lab$Description, function(s) paste(strwrap(.short_label(s, 44), 22), collapse = "\n"), character(1)))
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = e2, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = jaccard), colour = "#DCDAD5") +
    ggplot2::geom_point(data = nodes, ggplot2::aes(x, y, size = size, fill = colour_group), shape = 21, colour = "white", stroke = 0.6) +
    do.call(ggrepel::geom_text_repel, c(list(data = rep_lab, mapping = ggplot2::aes(x, y, label = lab), size = 2.25,
      family = .pai_font(), colour = .tok$ink, lineheight = 0.9, box.padding = 0.6, segment.size = 0.25, segment.colour = .tok$ink2,
      max.overlaps = Inf, seed = seed), .halo)) +
    ggplot2::scale_fill_manual(values = pal, guide = "none") +
    ggplot2::scale_linewidth_continuous(range = c(0.2, 1.2), guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.2, 4.5), name = "Targets in the term", breaks = function(l) unique(round(seq(l[1], l[2], length.out = 3))),
      guide = ggplot2::guide_legend(title.position = "left", nrow = 1, override.aes = list(fill = .tok$muted))) +
    ggplot2::coord_cartesian(clip = "off") + theme_pai(grid = "none") +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(), axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(), legend.position = "bottom", legend.justification = "left")
  out <- patchwork::wrap_plots(list(.pai_header(title %||% .hub_title(h, "Enrichment map"),
      sprintf("%d most significant terms  \u00b7  lines: Jaccard \u2265 %s  \u00b7  colour: module", nrow(t), h$params$similarity)),
    p, .hub_badge(h)), ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 1, 6.5), c("mm", "null", "mm")))
  attr(out, "plot_data") <- list(nodes = nodes, edges = ed)
  .with_size(out, 92, 110)
}
