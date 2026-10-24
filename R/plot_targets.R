# Target-level plots of an L1000PAIResult.
# Every function returns a ggplot or patchwork object; the values drawn are attached as
# attr(p, "plot_data") and a recommended size in millimetres as attr(p, "size_mm") (see pai_save()).

# ---- shared helpers ----------------------------------------------------------------------
.pai_df <- function(x) {
  if (methods::is(x, "L1000PAIResult")) as.data.frame(result_table(x)) else as.data.frame(x)
}
.pai_library <- function(x, library = NULL) {
  if (!is.null(library)) return(library)
  if (methods::is(x, "L1000PAIResult")) return(pai_metadata(x)$library %||% "unknown")
  df <- as.data.frame(x)
  if ("library_id" %in% names(df)) return(as.character(df$library_id[1]))
  "unknown"
}
.pai_score_direction <- function(x, library) {
  if (methods::is(x, "L1000PAIResult")) {
    d <- pai_metadata(x)$pairing_params$score_direction
    if (!is.null(d) && !is.na(d)) return(d)
  }
  resolve_score_direction(library = library)
}
.pai_semantics <- function(library, score_direction) {
  what <- switch(library, crispr = "knockout", oe = "overexpression", chem = "compound", single_drug = "drug",
    ligand = "ligand", kinase_gpcr = , kinase_gpcr_up = , kinase_gpcr_down = "perturbation", "perturbation")
  resembles <- paste0("resembles the ", what); reverse <- paste0("is the reverse of the ", what)
  drug_like <- library %in% c("chem", "single_drug")
  act <- if (identical(score_direction, "down_minus_up")) reverse else resembles
  inh <- if (identical(score_direction, "down_minus_up")) resembles else reverse
  short <- if (drug_like) c(Activated = "Concordant", Inhibited = "Opposing") else c(Activated = "Activated", Inhibited = "Inhibited")
  list(what = what,
    short = c(short, Discordant = "Discordant", Not_significant = "Not called"),
    quadrant = if (drug_like) c(Activated = "concordant", Inhibited = "opposing") else c(Activated = "activation", Inhibited = "inhibition"),
    explain = c(Activated = paste("your contrast", act), Inhibited = paste("your contrast", inh),
                Discordant = "both arms significant, same sign", Not_significant = "not called"),
    long = c(Activated = paste0(short[["Activated"]], ": your contrast ", act),
             Inhibited = paste0(short[["Inhibited"]], ": your contrast ", inh),
             Discordant = "Discordant: both arms significant, same sign", Not_significant = "Not called"),
    score_label = if (identical(score_direction, "down_minus_up")) "Activity score  (NES Down \u2212 NES Up)" else
      "Activity score  (NES Up \u2212 NES Down)",
    up_arm = paste0("genes up after the ", what), dn_arm = paste0("genes down after the ", what))
}
.pai_libname <- function(library) switch(library, crispr = "CRISPR knockout", oe = "overexpression",
  chem = "compound", single_drug = "drug", ligand = "ligand", kinase_gpcr = "kinase/GPCR", library)
.pai_rule_body <- function(x) {
  cp <- if (methods::is(x, "L1000PAIResult")) pai_metadata(x)$classify_params else NULL
  sprintf("coherent adj. P < %s  \u00b7  both arms aligned  \u00b7  |activity| \u2265 %s",
    cp$classify_p_value_threshold %||% 0.05, cp$activity_threshold %||% 1)
}
.fmt_n <- function(n) formatC(n, big.mark = ",", format = "d")
.fmt_p <- function(p) { s <- formatC(p, format = "e", digits = 1); sub("e-0?", "e\u2212", s) }
.neglog10 <- function(p, floor = 1e-10) -log10(pmax(p, floor))
.lab_minus <- function(v) sub("^-", "\u2212", format(v, trim = TRUE, drop0trailing = TRUE))
.fmt_num <- function(v, d = 2) sub("^-", "\u2212", formatC(v, format = "f", digits = d))
.sup <- function(s) chartr("0123456789-", "\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079\u207b", s)
# plotmath version, for geom_text(parse = TRUE): Liberation/Arial lack superscript 0 and minus.
.fmt_p10_expr <- function(p, prefix = NULL) {
  pre <- if (is.null(prefix)) rep("", length(p)) else sprintf('"%s" ~ ', rep_len(prefix, length(p)))
  unname(mapply(function(v, pr) {
    if (is.na(v)) return(paste0(pr, '"NA"'))
    if (v >= 0.01) return(sprintf('%s"%s"', pr, formatC(v, format = "f", digits = 2)))
    parts <- strsplit(formatC(v, format = "e", digits = 1), "e")[[1]]
    sprintf('%s"%s" %%*%% 10^%d', pr, parts[1], as.integer(parts[2]))
  }, p, pre))
}
.fmt_p10 <- function(p) vapply(p, function(v) {           # 0.21, 3.5\u00d710\u207b\u2079 (table format)
  if (is.na(v)) return("NA")
  if (v >= 0.01) return(formatC(v, format = "f", digits = 2))
  m <- formatC(v, format = "e", digits = 1); parts <- strsplit(m, "e")[[1]]
  paste0(parts[1], "\u00d710", .sup(as.character(as.integer(parts[2]))))
}, character(1))
.halo <- list(bg.color = "white", bg.r = 0.12)

# a. Top calls -----------------------------------------------------------------------------
#' Top calls in each direction
#'
#' The strongest Activated and Inhibited targets of an L1000PAI result, ranked by |activity
#' score| (ties: coherent adjusted P, then name), drawn as lollipops in two panels.
#'
#' @inheritParams pai_plot_params
#' @param n Number of targets per direction.
#' @param highlight Targets to shade, e.g. `"TP53"`.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_dotplot <- function(x, n = 10, highlight = NULL, library = NULL, title = NULL) {
  lib <- .pai_library(x, library)
  sem <- .pai_semantics(lib, .pai_score_direction(x, lib))
  df <- .pai_df(x); df$activity_class <- as.character(df$activity_class)
  called <- df[df$activity_class %in% c("Activated", "Inhibited"), ]
  top <- do.call(rbind, lapply(c("Activated", "Inhibited"), function(cl) {
    d <- called[called$activity_class == cl, ]
    utils::head(d[order(-abs(d$activity_score), d$padj_coherent, d$term), ], n)
  }))
  top$label <- top$perturbagen %||% top$term
  n_cls <- table(factor(df$activity_class, levels = c("Activated", "Inhibited")))
  xmax <- max(abs(top$activity_score)) * 1.06
  panel <- function(cl, show_x) {
    d <- top[top$activity_class == cl, ]
    d$label <- factor(d$label, levels = d$label[order(d$activity_score)])
    col <- .pal$marks[[cl]]
    p <- ggplot2::ggplot(d, ggplot2::aes(x = activity_score, y = label))
    if (any(d$label %in% highlight)) p <- p + ggplot2::geom_tile(data = d[d$label %in% highlight, ],
      ggplot2::aes(x = 0, width = Inf, height = 0.86), fill = .tok$pill)
    p <- p + ggplot2::geom_vline(xintercept = 0, colour = .tok$axis, linewidth = 0.35) +
      ggplot2::geom_segment(ggplot2::aes(x = 0, xend = activity_score, yend = label), colour = .tint(col, 0.45),
        linewidth = 0.65, lineend = "round") +
      ggplot2::geom_point(shape = 21, fill = col, colour = "white", size = 2.9, stroke = 0.7) +
      ggplot2::scale_x_continuous(limits = c(-xmax, xmax), breaks = scales::breaks_pretty(5), labels = .lab_minus) +
      ggplot2::labs(x = if (show_x) sem$score_label else NULL, y = NULL) +
      theme_pai(grid = "x") + ggplot2::theme(axis.line.x = if (show_x) ggplot2::element_line(colour = .tok$axis, linewidth = 0.3) else ggplot2::element_blank())
    if (!show_x) p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    p + ggplot2::theme(axis.ticks.y = ggplot2::element_blank(), plot.margin = ggplot2::margin(0, 4, 1, 4))
  }
  vary <- diff(range(.neglog10(top$padj_coherent))) >= 0.5
  rule <- .pai_rule_body(x)
  sub1 <- sprintf('"%s of %s %s signatures called"', .fmt_n(sum(n_cls)), .fmt_n(nrow(df)), .pai_libname(lib))
  sub2 <- .fmt_p10_expr(max(top$padj_coherent), prefix = "All shown have coherent adj. P \u2264")
  sub2 <- paste0(sub2, ' ~ "(the P-value floor)"')
  na <- sum(top$activity_class == "Activated"); ni <- sum(top$activity_class == "Inhibited")
  parts <- list(
    .pai_header(title %||% sprintf("Top %d calls in each direction", n), if (vary) sub1 else c(sub1, sub2), subtitle_parse = TRUE),
    .pai_header(sprintf("%s targets", .fmt_n(n_cls[["Activated"]])), sem$explain[["Activated"]],
      cls = "Activated", pill = sem$short[["Activated"]], size = 7.2),
    panel("Activated", FALSE),
    .pai_header(sprintf("%s targets", .fmt_n(n_cls[["Inhibited"]])), sem$explain[["Inhibited"]],
      cls = "Inhibited", pill = sem$short[["Inhibited"]], size = 7.2),
    panel("Inhibited", TRUE),
    .pai_badge("Call rule", rule))
  out <- patchwork::wrap_plots(parts, ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(if (vary) 1 else 2) + 1.5, 8.6, na, 8.6, ni + 1.6, 6),
      c("mm", "mm", "null", "mm", "null", "mm")))
  attr(out, "plot_data") <- top
  .with_size(out, 92, round(30 + (na + ni) * 4.7, 1))
}

# b. Why one target was called ----------------------------------------------------------------
.running_es <- function(stats, set, exponent = 1) {
  hit <- names(stats) %in% set
  w <- abs(stats)^exponent; nr <- sum(w[hit])
  res <- cumsum(ifelse(hit, w / nr, 0)) - cumsum(ifelse(hit, 0, 1 / (length(stats) - sum(hit))))
  list(position = seq_along(stats), running = res, hit = which(hit), es = res[which.max(abs(res))],
       peak = which.max(abs(res)))
}

#' Why one target was called: the two running-sum curves
#'
#' Running enrichment scores of the target's up arm (solid line) and down arm (dashed line) along
#' your ranked genes. Colour shows which way each arm points (toward activation or inhibition).
#'
#' @inheritParams pai_plot_params
#' @param term The target, as in the `term` column of [result_table()].
#' @param term2gene The signature library used for the run; by default it is read again with
#'   [get_l1000_library()].
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_gseaplot <- function(x, term, term2gene = NULL, library = NULL, title = NULL) {
  stopifnot(methods::is(x, "L1000PAIResult"))
  lib <- .pai_library(x, library)
  sem <- .pai_semantics(lib, .pai_score_direction(x, lib))
  df <- .pai_df(x)
  row <- df[df$term == term, , drop = FALSE]
  if (nrow(row) != 1L) stop("`term` must match exactly one row of the result table.")
  if (is.null(term2gene)) term2gene <- L1000PAI::get_l1000_library(lib)
  stats <- sort(rank_table(x), decreasing = TRUE)
  exponent <- pai_metadata(x)$gsea_params$exponent %||% 1
  arms <- list(up = row$term_Up, dn = row$term_Down)
  cv <- lapply(arms, function(a) .running_es(stats, term2gene$gene[term2gene$term == a], exponent))
  n <- length(stats); cls <- as.character(row$activity_class); who <- row$perturbagen %||% term
  up_txt <- sub("the ", paste0(who, " "), sem$up_arm); dn_txt <- sub("the ", paste0(who, " "), sem$dn_arm)
  curve <- do.call(rbind, lapply(names(cv), function(a) data.frame(arm = a, position = cv[[a]]$position, running = cv[[a]]$running)))
  hits <- do.call(rbind, lapply(names(cv), function(a) data.frame(arm = a, position = cv[[a]]$hit)))
  peaks <- do.call(rbind, lapply(names(cv), function(a) data.frame(arm = a, position = cv[[a]]$peak, running = cv[[a]]$es)))
  peaks$lab <- sprintf("ES %s%.2f", ifelse(peaks$running < 0, "\u2212", ""), abs(peaks$running))
  curve$arm <- factor(curve$arm, c("up", "dn")); hits$arm <- factor(hits$arm, c("up", "dn"))
  # House style: colour = which way the arm points (toward activation or inhibition); line type = which arm.
  o <- .pai_orient(row, .pai_score_direction(x, lib))
  arm_dir <- c(up = if (o$x >= 0) "Activated" else "Inhibited", dn = if (o$y >= 0) "Activated" else "Inhibited")
  arm_col <- c(up = .pal$marks[[arm_dir[["up"]]]], dn = .pal$marks[[arm_dir[["dn"]]]])
  xs <- ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(1, n), labels = scales::comma)
  p1 <- ggplot2::ggplot(curve, ggplot2::aes(position, running)) +
    ggplot2::geom_hline(yintercept = 0, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::geom_area(ggplot2::aes(fill = arm), alpha = 0.10, position = "identity") +
    ggplot2::geom_line(ggplot2::aes(colour = arm, linetype = arm), linewidth = 0.75, lineend = "butt", linejoin = "round") +
    ggplot2::geom_point(data = peaks, ggplot2::aes(fill = arm), shape = 21, colour = "white", size = 2.6, stroke = 0.7) +
    ggplot2::geom_text(data = peaks, ggplot2::aes(label = lab, vjust = ifelse(running > 0, -0.9, 1.8)),
      size = 2.3, colour = .tok$ink, family = .pai_font(), fontface = "bold") +
    ggplot2::scale_colour_manual(values = arm_col, guide = "none") +
    ggplot2::scale_fill_manual(values = arm_col, guide = "none") +
    ggplot2::scale_linetype_manual(values = c(up = "solid", dn = "22"), guide = "none") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.16, 0.16)), labels = .lab_minus) + xs +
    ggplot2::labs(x = NULL, y = "Running enrichment score") + theme_pai(grid = "y") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank())
  p2 <- ggplot2::ggplot(hits, ggplot2::aes(x = position, y = arm, colour = arm)) +
    ggplot2::geom_point(shape = "|", size = 2.1, alpha = 0.85) +
    ggplot2::scale_colour_manual(values = arm_col, guide = "none") +
    ggplot2::scale_y_discrete(limits = c("dn", "up"), labels = c(up = "Up arm", dn = "Down arm")) + xs +
    ggplot2::labs(x = NULL, y = NULL) + theme_pai(grid = "none") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank())
  rk <- data.frame(position = seq_len(n), stat = unname(stats))
  zc <- which(diff(sign(stats)) != 0)[1]
  p3 <- ggplot2::ggplot(rk, ggplot2::aes(position, stat)) +
    ggplot2::geom_area(fill = "#CFCCC6") +
    ggplot2::geom_hline(yintercept = 0, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::annotate("text", x = n * 0.05, y = max(stats) * 0.55, label = "up in your contrast", hjust = 0,
      size = 2.2, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::annotate("text", x = n * 0.95, y = min(stats) * 0.55, label = "down in your contrast", hjust = 1,
      size = 2.2, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::annotate("segment", x = zc, xend = zc, y = min(stats) * 0.2, yend = max(stats) * 0.2,
      colour = .tok$muted, linewidth = 0.3) +
    xs + ggplot2::scale_y_continuous(breaks = scales::breaks_extended(3), labels = .lab_minus) +
    ggplot2::labs(x = sprintf("Rank in your contrast  (%s genes)", .fmt_n(n)), y = NULL) + theme_pai(grid = "none")
  head <- .pai_header(title %||% who,
    .fmt_p10_expr(row$padj_coherent, prefix = sprintf("NES Up %s  \u00b7  NES Down %s  \u00b7  activity %s  \u00b7  adj. P =",
      .fmt_num(row$NES_Up), .fmt_num(row$NES_Down), .fmt_num(row$activity_score))), subtitle_parse = TRUE,
    cls = cls, pill = sem$short[[cls]])
  where <- function(nes) if (nes > 0) "rise" else "fall"
  read <- .pai_badge("In your contrast", sprintf("%s %s up-genes %s, down-genes %s",
    who, sem$what, where(row$NES_Up), where(row$NES_Down)))
  toward <- c(Activated = paste("toward", sem$quadrant[["Activated"]]), Inhibited = paste("toward", sem$quadrant[["Inhibited"]]))
  chips <- .pai_chips(stats::setNames(arm_col, c(sprintf("Up arm (up after %s): %s", sem$what, toward[[arm_dir[["up"]]]]),
    sprintf("Down arm (down after %s): %s", sem$what, toward[[arm_dir[["dn"]]]]))), shapes = c("line", "dash"), size = 6, per_row = 1)
  out <- .stack(list(head, read, chips, p1, p2, p3), c(.header_mm(1), 6, 9, NA, NA, NA)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 6, 9, 3, 0.62, 1.15),
      c("mm", "mm", "mm", "null", "null", "null")))
  attr(out, "running_es") <- c(up = unname(cv$up$es), dn = unname(cv$dn$es))
  attr(out, "arm_direction") <- arm_dir
  attr(out, "plot_data") <- list(curves = curve, hits = hits, peaks = peaks, ranks = rk)
  .with_size(out, 92, 100)
}

# c. Both arms for every target: the oriented arm plane ------------------------------
# Each arm is signed so that "towards activation" is positive. The activity score is then the
# sum of the two coordinates, activation sits top-right, inhibition bottom-left, and the two
# off-diagonal quadrants hold pairs whose arms disagree.
.pai_orient <- function(df, sd) {
  if (identical(sd, "down_minus_up"))
    list(x = -df$NES_Up, y = df$NES_Down, xlab = "Up arm, oriented  (\u2212NES Up)", ylab = "Down arm  (NES Down)")
  else
    list(x = df$NES_Up, y = -df$NES_Down, xlab = "Up arm  (NES Up)", ylab = "Down arm, oriented  (\u2212NES Down)")
}
.pai_threshold <- function(x) {
  cp <- if (methods::is(x, "L1000PAIResult")) pai_metadata(x)$classify_params else NULL
  cp$activity_threshold %||% 1
}
# Invisible points along the four quadrant captions drawn by .arm_canvas(), added to the data of
# a repelled-label layer (with empty labels) so that the labels avoid the captions.
.arm_label_obstacles <- function(lim, sem, label_size = 2.05) {
  txt <- c(paste("both arms:", sem$quadrant[["Activated"]]), paste("both arms:", sem$quadrant[["Inhibited"]]),
           "arms disagree", "arms disagree")
  x1 <- c(0.96, -0.96, -0.96, 0.96) * lim; dirn <- c(-1, 1, 1, -1); yc <- c(0.94, -0.06, 0.94, -0.06) * lim
  cw <- 0.028 * lim * label_size / 2.05                       # about one character of the caption font
  do.call(rbind, lapply(1:4, function(i) {
    xs <- seq(0, nchar(txt[i]) * cw, by = 0.04 * lim)
    data.frame(x_oriented = x1[i] + dirn[i] * xs, y_oriented = yc[i], lab = "", stringsAsFactors = FALSE)
  }))
}
.arm_canvas <- function(lim, sem, thr, label_size = 2.05, iso = TRUE, labels = TRUE) {
  out <- list(
    ggplot2::annotate("rect", xmin = 0, xmax = lim, ymin = 0, ymax = lim, fill = .pal$tint[["Activated"]], alpha = 0.75),
    ggplot2::annotate("rect", xmin = -lim, xmax = 0, ymin = -lim, ymax = 0, fill = .pal$tint[["Inhibited"]], alpha = 0.75),
    ggplot2::geom_hline(yintercept = 0, colour = .tok$ink2, linewidth = 0.2),
    ggplot2::geom_vline(xintercept = 0, colour = .tok$ink2, linewidth = 0.2))
  if (iso) {
    iy <- 0.70 * lim * c(1, -1)
    out <- c(out, list(
      ggplot2::geom_abline(intercept = c(thr, -thr), slope = -1, colour = .tok$muted, linewidth = 0.22, linetype = "22"),
      ggplot2::annotate("text", x = c(thr, -thr) - iy, y = iy, label = c(sprintf("score = +%s", thr), sprintf("score = \u2212%s", thr)),
        angle = -45, vjust = -0.45, size = label_size * 0.88, colour = .tok$muted, family = .pai_font())))
  }
  # the two lower labels sit just under the x axis, leaving the corners to the data
  if (labels) out <- c(out, list(ggplot2::annotate("text", x = c(0.96, -0.96, -0.96, 0.96) * lim, y = c(0.96, -0.04, 0.96, -0.04) * lim,
    hjust = c(1, 0, 0, 1), vjust = c(1, 1, 1, 1),
    label = c(paste("both arms:", sem$quadrant[["Activated"]]), paste("both arms:", sem$quadrant[["Inhibited"]]), "arms disagree", "arms disagree"),
    colour = c(.pal$text[["Activated"]], .pal$text[["Inhibited"]], .pal$text[["Discordant"]], .pal$text[["Discordant"]]),
    size = label_size, family = .pai_font())))
  out
}

#' Both arms of every target: the oriented arm plane
#'
#' Each arm is signed so that right/up means toward activation; the activity score is then the
#' sum of the two coordinates. Activated targets sit in the warm top-right quadrant, Inhibited
#' targets in the cool bottom-left one, and pairs whose arms disagree in the other two.
#'
#' @inheritParams pai_plot_params
#' @param label Targets to label in addition to the strongest calls.
#' @param n_label Number of the strongest calls to label (half Activated, half Inhibited).
#' @param subtitle Subtitle; `NULL` for the default.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_armplot <- function(x, label = NULL, n_label = 6, library = NULL, title = NULL, subtitle = NULL) {
  lib <- .pai_library(x, library); sd <- .pai_score_direction(x, lib); sem <- .pai_semantics(lib, sd)
  df <- .pai_df(x); df <- df[!is.na(df$NES_Up) & !is.na(df$NES_Down), ]
  o <- .pai_orient(df, sd); df$x_oriented <- o$x; df$y_oriented <- o$y
  df$activity_class <- factor(as.character(df$activity_class), levels = c("Not_significant", "Discordant", "Inhibited", "Activated"))
  df <- df[order(df$activity_class), ]
  counts <- table(df$activity_class)
  called <- df[df$activity_class %in% c("Activated", "Inhibited"), ]
  lab <- called[order(-abs(called$activity_score)), ]
  lab <- rbind(utils::head(lab[lab$activity_class == "Activated", ], ceiling(n_label / 2)),
               utils::head(lab[lab$activity_class == "Inhibited", ], floor(n_label / 2)))
  gene <- df$perturbagen %||% df$term
  if (!is.null(label)) lab <- unique(rbind(lab, df[df$term %in% label | gene %in% label, ]))
  lab$lab <- lab$perturbagen %||% lab$term
  thr <- .pai_threshold(x)
  lim <- ceiling(max(abs(c(df$x_oriented, df$y_oriented))) * 1.02 * 2) / 2
  shapes <- c(Not_significant = 16, Discordant = 18, Inhibited = 16, Activated = 16)
  sizes <- c(Not_significant = 0.5, Discordant = 0.85, Inhibited = 0.62, Activated = 0.62)
  p <- ggplot2::ggplot(df, ggplot2::aes(x_oriented, y_oriented)) +
    .arm_canvas(lim, sem, thr) +
    ggplot2::geom_point(ggplot2::aes(colour = activity_class, shape = activity_class, size = activity_class),
      alpha = 0.8, stroke = 0) +
    do.call(ggrepel::geom_text_repel, c(list(data = rbind(transform(lab[, c("x_oriented", "y_oriented", "lab")], ps = 2.4),
        transform(.arm_label_obstacles(lim, sem), ps = 0.3)),
      mapping = ggplot2::aes(label = lab, point.size = ps), size = 2.3, family = .pai_font(),
      fontface = "bold", colour = .tok$ink, min.segment.length = 0, segment.size = 0.25, segment.colour = .tok$ink,
      box.padding = 0.4, point.padding = 0.2, force = 4, force_pull = 0.3, max.time = 3, max.iter = 1e5,
      max.overlaps = Inf, seed = 1), .halo)) +
    ggplot2::geom_point(data = lab, ggplot2::aes(fill = activity_class), shape = ifelse(lab$activity_class == "Discordant", 23, 21),
      size = 1.9, stroke = 0.4, colour = .tok$ink) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_shape_manual(values = shapes, guide = "none") +
    ggplot2::scale_size_manual(values = sizes, guide = "none") +
    ggplot2::scale_x_continuous(breaks = scales::breaks_width(2), labels = .lab_minus) +
    ggplot2::scale_y_continuous(breaks = scales::breaks_width(2), labels = .lab_minus) +
    ggplot2::coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    ggplot2::labs(x = o$xlab, y = o$ylab) +
    theme_pai(axes = "xy")
  chips <- .class_chips(sem, counts, per_row = 2)
  out <- .stack(list(.pai_header(title %||% "Directional NES geometry",
      subtitle %||% sprintf("%s %s pairs  \u00b7  score = sum of the two oriented arms", .fmt_n(nrow(df)), .pai_libname(lib))),
    chips, p, .pai_badge("Call rule", .pai_rule_body(x))), c(.header_mm(1), 9, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 0.5, 9, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- df
  .with_size(out, 92, 118)
}

# d. Several contrasts side by side ---------------------------------------------------------
#' The same targets in several results
#'
#' Targets called in the most results, one column per result; the bar on the right gives how many
#' results make the same call.
#'
#' @inheritParams pai_plot_params
#' @param results A named list of L1000PAIResult objects from the same library; the names label
#'   the columns.
#' @param terms Targets to show; `NULL` takes the `n` targets called most often in each direction.
#' @param n Number of targets per direction.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_compareplot <- function(results, terms = NULL, n = 8, library = NULL, title = NULL) {
  stopifnot(is.list(results), !is.null(names(results)))
  lib <- .pai_library(results[[1]], library); sem <- .pai_semantics(lib, .pai_score_direction(results[[1]], lib))
  long <- do.call(rbind, lapply(names(results), function(nm) {
    d <- .pai_df(results[[nm]])
    data.frame(study = nm, term = d$term, label = d$perturbagen %||% d$term, activity_score = d$activity_score,
      padj = d$padj_coherent, cls = as.character(d$activity_class))
  }))
  long$called <- long$cls %in% c("Activated", "Inhibited"); k <- length(results)
  agg <- do.call(rbind, lapply(split(long, long$term), function(d) data.frame(term = d$term[1],
    n_act = sum(d$cls == "Activated"), n_inh = sum(d$cls == "Inhibited"), mean_score = mean(d$activity_score, na.rm = TRUE))))
  if (is.null(terms)) {
    a <- agg[agg$n_act > 0, ]; a <- a[order(-a$n_act, -a$mean_score), ]
    i <- agg[agg$n_inh > 0, ]; i <- i[order(-i$n_inh, i$mean_score), ]
    terms <- c(utils::head(a$term, n), utils::head(i$term, n))
  }
  agg <- agg[agg$term %in% terms, ]
  agg$consensus <- ifelse(agg$n_act >= agg$n_inh, "Activated", "Inhibited")
  agg$n_same <- ifelse(agg$consensus == "Activated", agg$n_act, agg$n_inh)
  d <- merge(long[long$term %in% terms, ], agg[, c("term", "consensus", "n_same", "mean_score")], by = "term")
  d$term <- factor(d$term, levels = agg$term[order(agg$mean_score)])
  study_lab <- stats::setNames(names(results), names(results))
  d$study <- factor(d$study, levels = names(results))
  d$panel <- factor(d$consensus, levels = c("Activated", "Inhibited"))
  d$neglogp <- .neglog10(d$padj)
  lim <- max(abs(d$activity_score), na.rm = TRUE)
  cnt <- unique(d[, c("term", "n_same", "panel", "consensus")])
  bx <- k + 0.75; bw <- 0.9
  p <- ggplot2::ggplot(d, ggplot2::aes(x = as.numeric(study), y = term)) +
    ggplot2::geom_segment(data = cnt, inherit.aes = FALSE, ggplot2::aes(x = bx, xend = bx + bw, y = term, yend = term),
      colour = .tok$pill, linewidth = 1.9, lineend = "round") +
    ggplot2::geom_segment(data = cnt, inherit.aes = FALSE, ggplot2::aes(x = bx, xend = bx + bw * n_same / k, y = term, yend = term,
      colour = consensus), linewidth = 1.9, lineend = "round") +
    ggplot2::geom_text(data = cnt, inherit.aes = FALSE, ggplot2::aes(x = bx + bw + 0.12, y = term, label = sprintf("%d/%d", n_same, k)),
      hjust = 0, size = 2.2, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::geom_point(data = d[!d$called, ], ggplot2::aes(size = neglogp), shape = 21, fill = "white", colour = .tok$axis, stroke = 0.45) +
    ggplot2::geom_point(data = d[d$called, ], ggplot2::aes(size = neglogp, fill = activity_score), shape = 21, colour = "white", stroke = 0.5) +
    ggplot2::facet_grid(panel ~ ., scales = "free_y", space = "free_y", labeller = ggplot2::as_labeller(sem$short)) +
    ggplot2::scale_fill_gradient2(low = .pal$marks[["Inhibited"]], mid = .pai_mid, high = .pal$marks[["Activated"]],
      limits = c(-lim, lim), name = "Activity score", labels = .lab_minus, guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(22, "mm"),
      barheight = ggplot2::unit(2, "mm"), title.position = "top")) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_size_continuous(name = "\u2212log10 coherent adj. P", range = c(0.9, 3.4), breaks = scales::breaks_pretty(3),
      guide = ggplot2::guide_legend(title.position = "top", nrow = 1, override.aes = list(fill = .tok$muted, colour = "white"))) +
    ggplot2::scale_x_continuous(breaks = seq_len(k), labels = study_lab[names(results)], limits = c(0.5, bx + bw + 0.75),
      expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) +
    theme_pai(grid = "y") +
    ggplot2::theme(legend.position = "bottom", legend.justification = "left", legend.box.just = "left",
      axis.line.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(colour = .tok$ink, lineheight = 0.95, angle = 30, hjust = 1, vjust = 1),
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold", hjust = 0, colour = .tok$ink),
      panel.spacing.y = ggplot2::unit(2.5, "mm"))
  chips <- .pai_chips(c(stats::setNames(.pal$marks[c("Activated", "Inhibited")], sprintf("Mostly %s", tolower(sem$short[c("Activated", "Inhibited")]))),
    `Not called (open)` = .tok$axis))
  out <- .stack(list(.pai_header(title %||% "Targets called in the most studies",
      "Bars on the right: studies giving the same call"), chips, p, .pai_badge("Call rule", .pai_rule_body(results[[1]]))),
    c(.header_mm(1), 5, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1), 5.2, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- d
  .with_size(out, max(92, 60 + 8 * k), round(50 + length(unique(d$term)) * 5, 1))
}

# e. Which genes drive the calls --------------------------------------------------------------
#' Leading-edge genes behind the calls
#'
#' Targets (squares) linked to the genes of your contrast that drive their two arms (leading
#' edge); solid links come from the up arm, dashed links from the down arm, and genes are coloured
#' by their statistic in your contrast.
#'
#' @inheritParams pai_plot_params
#' @param terms Targets to draw; `NULL` takes the `n` strongest calls (half in each direction).
#' @param n Number of targets when `terms` is `NULL`.
#' @param n_genes Leading-edge genes per target, largest |statistic| first.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_cnetplot <- function(x, terms = NULL, n = 4, n_genes = 6, library = NULL, title = NULL, seed = 1) {
  stopifnot(methods::is(x, "L1000PAIResult"))
  lib <- .pai_library(x, library); sem <- .pai_semantics(lib, .pai_score_direction(x, lib))
  df <- .pai_df(x); stats <- rank_table(x)
  if (is.null(terms)) {
    called <- df[as.character(df$activity_class) %in% c("Activated", "Inhibited"), ]
    called <- called[order(-abs(called$activity_score)), ]
    terms <- c(utils::head(called$term[called$activity_class == "Activated"], ceiling(n / 2)),
               utils::head(called$term[called$activity_class == "Inhibited"], floor(n / 2)))
  }
  edges <- do.call(rbind, lapply(terms, function(tm) {
    r <- df[df$term == tm, ]
    up <- strsplit(r$core_enrichment_Up, "/")[[1]]; dn <- strsplit(r$core_enrichment_Down, "/")[[1]]
    le <- unique(c(up, dn)); le <- le[le %in% names(stats)]
    le <- utils::head(le[order(-abs(stats[le]))], n_genes)
    data.frame(from = r$perturbagen %||% tm, to = le, cls = as.character(r$activity_class),
      arm = ifelse(le %in% up, "up", "dn"))
  }))
  tnodes <- unique(edges[, c("from", "cls")]); gnodes <- unique(edges$to)
  g <- igraph::graph_from_data_frame(edges[, c("from", "to")], directed = FALSE,
    vertices = data.frame(name = c(tnodes$from, setdiff(gnodes, tnodes$from))))
  set.seed(seed); lay <- igraph::layout_with_fr(g)
  nodes <- data.frame(name = igraph::V(g)$name, x = lay[, 1], y = lay[, 2])
  nodes$type <- ifelse(nodes$name %in% tnodes$from, "target", "gene")
  nodes$cls <- tnodes$cls[match(nodes$name, tnodes$from)]
  nodes$stat <- ifelse(nodes$type == "gene", stats[nodes$name], NA)
  nodes$degree <- igraph::degree(g)[nodes$name]
  e <- merge(merge(edges, nodes[, c("name", "x", "y")], by.x = "from", by.y = "name"),
             nodes[, c("name", "x", "y")], by.x = "to", by.y = "name", suffixes = c("", "end"))
  gn <- nodes[nodes$type == "gene", ]; tn <- nodes[nodes$type == "target", ]
  gn$abs <- abs(gn$stat); gn$shared <- gn$degree > 1
  p <- ggplot2::ggplot() +
    ggplot2::geom_curve(data = e, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = cls, linetype = arm),
      curvature = 0.12, linewidth = 0.45, alpha = 0.8) +
    ggplot2::geom_point(data = gn, ggplot2::aes(x, y, size = abs, fill = stat), shape = 21, colour = "white", stroke = 0.6) +
    ggplot2::scale_fill_gradient2(low = .pal$marks[["Inhibited"]], mid = .pai_mid, high = .pal$marks[["Activated"]],
      limits = c(-1, 1) * max(abs(gn$stat)), labels = .lab_minus, name = "Gene in your contrast  (rank statistic)",
      guide = ggplot2::guide_colourbar(barwidth = ggplot2::unit(22, "mm"), barheight = ggplot2::unit(2, "mm"), title.position = "top",
        order = 1)) +
    ggnewscale::new_scale_fill() +
    ggplot2::geom_point(data = tn, ggplot2::aes(x, y, fill = cls), shape = 22, size = 5.2, colour = "white", stroke = 0.9) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_linetype_manual(values = c(up = "solid", dn = "22"), guide = "none") +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.4, 3.6), name = "|gene rank statistic|", breaks = scales::breaks_pretty(3),
      guide = ggplot2::guide_legend(title.position = "top", nrow = 1, override.aes = list(shape = 16, colour = .tok$muted),
        order = 2)) +

    do.call(ggrepel::geom_text_repel, c(list(data = rbind(
        data.frame(x = gn$x, y = gn$y, name = gn$name, ff = ifelse(gn$shared, "bold", "plain"), al = 0.72,
          ps = 3, stringsAsFactors = FALSE),
        data.frame(x = tn$x, y = tn$y, name = tn$name, ff = "bold", al = 1, ps = 5.2, stringsAsFactors = FALSE)),
      # one layer, one font size and colour (a second size or colour scale would break the legends);
      # targets bold and fully opaque, genes lighter
      mapping = ggplot2::aes(x, y, label = name, fontface = ff, alpha = al, point.size = ps), size = 2.3, colour = .tok$ink,
      family = .pai_font(), segment.size = 0.22, segment.colour = .tok$axis, box.padding = 0.35, point.padding = 0.2,
      force = 3, max.time = 3, max.iter = 1e5, max.overlaps = Inf, seed = seed), .halo)) +
    ggplot2::scale_alpha_identity() +
    ggplot2::coord_cartesian(clip = "off") +
    theme_pai(grid = "none") +
    ggplot2::theme(axis.line.x = ggplot2::element_blank(), axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(), axis.title = ggplot2::element_blank(),
      legend.position = "bottom", legend.justification = "left")
  chips1 <- .pai_chips(stats::setNames(.pal$marks[c("Activated", "Inhibited")], sprintf("%s target and its links", sem$short[c("Activated", "Inhibited")])))
  chips2 <- .pai_chips(stats::setNames(rep(.tok$ink2, 2), c("via up-arm leading edge", "via down-arm leading edge")), shapes = c("line", "dash"))
  out <- .stack(list(.pai_header(title %||% "Leading-edge genes behind the calls",
      sprintf("Up to %d genes per target; bold = shared; colour: up/down in your contrast", n_genes)),
    chips1, chips2, p), c(.header_mm(1), 5, 5, NA)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1), 4.8, 4.8, 1), c("mm", "mm", "mm", "null")))
  e$stat <- unname(stats[e$to])
  attr(out, "plot_data") <- e
  .with_size(out, 92, 100)
}

# f. Drug families (Chem) ----------------------------------------------------------------------
#' Drug families among the calls (Chem libraries)
#'
#' For each compound family, the share of its compounds called Concordant or Opposing. Counts
#' only; no family-level test is made.
#'
#' @inheritParams pai_plot_params
#' @param n Number of families (half most Concordant, half most Opposing).
#' @param min_members Families with fewer compounds in the library are left out.
#' @param family_level,evidence_tier Passed to [annotate_chem_families()].
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family target plots
#' @export
pai_familyplot <- function(x, n = 12, min_members = 3, family_level = "direct_moa", evidence_tier = c("A", "B"), title = NULL) {
  lib <- .pai_library(x); sem <- .pai_semantics(lib, .pai_score_direction(x, lib))
  if (!lib %in% c("chem")) stop("pai_familyplot() needs a Chem result.")
  ann <- annotate_chem_families(.pai_df(x), family_level = family_level, evidence_tier = evidence_tier, include_unmapped = FALSE)
  ann <- ann[!is.na(ann$family_id), ]
  unit <- if ("analysis_unit_id" %in% names(ann)) "analysis_unit_id" else "term"
  fam_label <- if ("family_label" %in% names(ann)) "family_label" else "family_id"
  per <- do.call(rbind, lapply(split(ann, ann$family_id), function(d) {
    if ("preferred_for_analysis" %in% names(d)) d <- d[order(!as.logical(d$preferred_for_analysis)), ]
    u <- d[!duplicated(d[[unit]]), ]
    data.frame(family = u[[fam_label]][1], n = nrow(u), conc = sum(u$activity_class == "Activated"), opp = sum(u$activity_class == "Inhibited"))
  }))
  per <- per[per$n >= min_members, ]
  per$frac_conc <- per$conc / per$n; per$frac_opp <- per$opp / per$n
  top <- unique(rbind(utils::head(per[order(-per$frac_conc, -per$n), ], ceiling(n / 2)),
                      utils::head(per[order(-per$frac_opp, -per$n), ], floor(n / 2))))
  top$family <- factor(top$family, levels = top$family[order(top$frac_conc - top$frac_opp)])
  bars <- rbind(data.frame(family = top$family, cls = "Activated", value = top$frac_conc, k = top$conc, n = top$n),
                data.frame(family = top$family, cls = "Inhibited", value = -top$frac_opp, k = top$opp, n = top$n))
  p <- ggplot2::ggplot(bars, ggplot2::aes(y = family)) +
    ggplot2::geom_segment(data = top, ggplot2::aes(x = -1, xend = 1, y = family, yend = family), colour = .tok$pill,
      linewidth = 2.6, lineend = "butt") +
    ggplot2::geom_segment(data = bars[bars$k > 0, ], ggplot2::aes(x = 0, xend = value, yend = family, colour = cls),
      linewidth = 2.6, lineend = "butt") +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::geom_text(data = bars[bars$k > 0, ], ggplot2::aes(x = sign(value) * 1.05, label = sprintf("%d/%d", k, n),
      hjust = ifelse(value > 0, 0, 1)), size = 2.2, colour = .tok$ink, family = .pai_font()) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_x_continuous(labels = function(v) scales::percent(abs(v)), limits = c(-1.3, 1.3), breaks = c(-1, -0.5, 0, 0.5, 1)) +
    ggplot2::labs(x = "Share of the family's compounds called", y = NULL) +
    theme_pai(grid = "x") + ggplot2::theme(axis.text.y = ggplot2::element_text(colour = .tok$ink), axis.ticks.y = ggplot2::element_blank())
  chips <- .pai_chips(stats::setNames(.pal$marks[c("Inhibited", "Activated")],
    c(sprintf("%s \u00b7 the reverse", sem$short[["Inhibited"]]), sprintf("%s \u00b7 resembles the compound", sem$short[["Activated"]]))))
  out <- .stack(list(.pai_header(title %||% "Drug families among the calls",
      sprintf("Families with \u2265 %d compounds in the library; one vote per compound", min_members)),
    chips, p, .pai_badge("Counts only", "no family-level test")), c(.header_mm(1), 5, NA, 6.5)) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1), 5.2, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- bars
  .with_size(out, 92, round(30 + nrow(top) * 6, 1))
}

