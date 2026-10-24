# One call from a differential-expression table to a PDF report and an Excel workbook: pai_deg_report().
# The PDF is drawn page by page with grid; every size is in millimetres on an A4 page with 12 mm margins.

# Columns used inside ggplot2 aesthetics of the report plots.
utils::globalVariables(c("v", "side", "nes", "txt", "label", "k", "e", "lab"))

# =========================================================================================================
# 1. Input
# =========================================================================================================
.deg_cols <- list(
  gene = c("gene", "Gene", "symbol", "Symbol", "SYMBOL", "gene_symbol", "GeneSymbol", "hgnc_symbol", "gene_name",
           "Gene.symbol", "external_gene_name", "ID", "gene_id", "GeneID", "ensembl", "ensembl_gene_id", "entrez", "EntrezID"),
  stat = c("t", "stat", "t_stat", "statistic", "z", "zscore", "Wald"),
  fc = c("logFC", "log2FoldChange", "log2FC", "avg_log2FC", "logFC_log10ratio", "lfc"),
  p = c("P.Value", "pvalue", "PValue", "p_val", "p.value", "pval"),
  padj = c("adj.P.Val", "padj", "FDR", "p_val_adj", "qvalue", "q.value"),
  expr = c("AveExpr", "baseMean", "logCPM", "mean_expression"))
.rp_first_col <- function(d, key) { v <- intersect(.deg_cols[[key]], names(d))[1]; if (is.na(v)) NULL else v }

.deg_prepare <- function(deg, gene_col = NULL, stat_col = NULL, id_type = c("auto", "symbol", "ensembl", "entrez")) {
  id_type <- match.arg(id_type)
  deg <- as.data.frame(deg, stringsAsFactors = FALSE, check.names = FALSE)
  notes <- character()
  gene_col <- gene_col %||% .rp_first_col(deg, "gene")
  if (is.null(gene_col)) {
    deg$.row_names <- rownames(deg); gene_col <- ".row_names"; notes <- c(notes, "gene identifiers taken from the row names")
  }
  ids <- trimws(as.character(deg[[gene_col]]))
  if (id_type == "auto")
    id_type <- if (mean(grepl("^ENSG\\d+", ids)) > 0.5) "ensembl" else if (mean(grepl("^\\d+$", ids)) > 0.5) "entrez" else "symbol"
  fc_col <- .rp_first_col(deg, "fc"); p_col <- .rp_first_col(deg, "p")
  stat_col <- stat_col %||% .rp_first_col(deg, "stat")
  if (!is.null(stat_col)) {
    stat <- suppressWarnings(as.numeric(deg[[stat_col]])); stat_method <- sprintf("the column '%s'", stat_col)
  } else {
    if (is.null(fc_col) || is.null(p_col)) stop("No signed statistic (t, stat) and no fold change plus P value column found; set stat_col.", call. = FALSE)
    stat <- sign(as.numeric(deg[[fc_col]])) * -log10(pmax(as.numeric(deg[[p_col]]), .Machine$double.xmin))
    stat_col <- paste(fc_col, p_col, sep = " + "); stat_method <- sprintf("sign(%s) \u00d7 \u2212log10(%s)", fc_col, p_col)
  }
  symbol <- ids
  if (id_type != "symbol") {
    keytype <- if (id_type == "ensembl") "ENSEMBL" else "ENTREZID"
    key <- if (id_type == "ensembl") sub("\\..*$", "", ids) else ids
    map <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(key), column = "SYMBOL",
      keytype = keytype, multiVals = "first"))
    symbol <- unname(map[key]); notes <- c(notes, sprintf("%s identifiers mapped to gene symbols with org.Hs.eg.db", id_type))
  }
  qc <- data.frame(input_row = seq_along(ids), input_id = ids, gene = symbol, statistic = stat, stringsAsFactors = FALSE)
  padj_col <- .rp_first_col(deg, "padj"); expr_col <- .rp_first_col(deg, "expr")
  if (!is.null(fc_col)) qc$log_fold_change <- suppressWarnings(as.numeric(deg[[fc_col]]))
  if (!is.null(padj_col)) qc$adj_p <- suppressWarnings(as.numeric(deg[[padj_col]]))
  if (!is.null(expr_col)) qc$expression <- suppressWarnings(as.numeric(deg[[expr_col]]))
  qc$status <- ifelse(is.na(qc$gene) | qc$gene == "", "no gene symbol", ifelse(!is.finite(qc$statistic), "no statistic", "used"))
  used <- which(qc$status == "used")
  o <- used[order(-abs(qc$statistic[used]))]
  dup <- o[duplicated(qc$gene[o])]
  qc$status[dup] <- "duplicate gene, smaller |statistic| dropped"
  keep <- setdiff(o, dup)
  rank_input <- data.frame(gene = qc$gene[keep], stat = qc$statistic[keep], stringsAsFactors = FALSE)
  rank_input <- rank_input[order(-rank_input$stat), , drop = FALSE]; rownames(rank_input) <- NULL
  list(rank_input = rank_input, qc = qc, gene_col = gene_col, stat_col = stat_col, stat_method = stat_method,
       id_type = id_type, fc_col = fc_col, padj_col = padj_col, expr_col = expr_col, notes = notes,
       counts = c(rows = nrow(qc), used = length(keep), no_symbol = sum(qc$status == "no gene symbol"),
                  no_stat = sum(qc$status == "no statistic"), duplicates = length(dup)))
}

# =========================================================================================================
# 2. Guard modules
# =========================================================================================================
# Proliferation component of the input ranking (Hallmark E2F targets and G2M checkpoint, fgsea).
# Evaluate expr with a fixed seed and restore the caller's random-number state afterwards.
.rp_with_seed <- function(seed, expr) {
  env <- globalenv()
  old <- if (exists(".Random.seed", envir = env, inherits = FALSE)) get(".Random.seed", envir = env) else NULL
  on.exit(if (is.null(old)) suppressWarnings(rm(".Random.seed", envir = env)) else assign(".Random.seed", old, envir = env))
  set.seed(seed)
  expr
}
.deg_proliferation <- function(rank_input, sets, seed = 20260813) {
  r <- stats::setNames(rank_input$stat, rank_input$gene)
  g <- .rp_with_seed(seed, suppressWarnings(fgsea::fgseaMultilevel(sets, r, minSize = 10, maxSize = 500, eps = 0)))
  if (!nrow(g)) return(data.frame(gene_set = names(sets), size = 0L, nes = NA_real_, pvalue = NA_real_, padj = NA_real_,
                                  strong = FALSE, leading_edge = "", stringsAsFactors = FALSE))
  out <- data.frame(gene_set = g$pathway, size = g$size, nes = g$NES, pvalue = g$pval, padj = g$padj, stringsAsFactors = FALSE)
  out$strong <- !is.na(out$padj) & out$padj < 0.05 & abs(out$nes) >= 2
  out$leading_edge <- vapply(g$leadingEdge, function(z) paste(utils::head(z, 30), collapse = "/"), "")
  out[order(out$gene_set), , drop = FALSE]
}

# =========================================================================================================
# 3. Tables (the same tables fill the PDF and the Excel sheets)
# =========================================================================================================
.rp_call_label <- function(cls, lib) {
  cls <- as.character(cls)
  lab <- if (lib == "chem") c(Activated = "Concordant", Inhibited = "Opposing", Discordant = "Discordant", Not_significant = "Not called")
    else c(Activated = "Activated", Inhibited = "Inhibited (KO-like)", Discordant = "Discordant", Not_significant = "Not called")
  unname(ifelse(cls %in% names(lab), lab[cls], cls))
}
.deg_target_table <- function(x, prep) {
  d <- as.data.frame(result_table(x))
  ri <- prep$rank_input; ri$pct <- rank(ri$stat) / nrow(ri)
  target <- strip_direction_suffix(d$perturbagen %||% d$term)
  m <- match(target, ri$gene)
  out <- data.frame(target = target, call = .rp_call_label(d$activity_class, "crispr"), activity_score = d$activity_score,
    activity_percentile = rank(d$activity_score) / nrow(d), padj_coherent = d$padj_coherent,
    nes_up_arm = d$NES_Up, nes_down_arm = d$NES_Down, padj_up_arm = d$p.adjust_Up, padj_down_arm = d$p.adjust_Down,
    p_activation = d$p_activation, p_inhibition = d$p_inhibition,
    target_measured = !is.na(m), target_statistic = ri$stat[m], target_percentile = ri$pct[m], stringsAsFactors = FALSE)
  out <- out[order(-out$activity_score), , drop = FALSE]; rownames(out) <- NULL
  out
}
.rp_conf_text <- c(corroborated = "2+ sources", single = "1 source", conflict = "sources disagree", author_curated = "author decision")
.deg_compound_table <- function(xc, reference) {
  d <- as.data.frame(result_table(xc))
  compound <- strip_direction_suffix(d$perturbagen %||% d$term)
  pct <- rank(d$activity_score) / nrow(d)
  ann <- get_l1000_annotation_catalog("chem")
  a <- ann[match(compound, ann$compound_name), , drop = FALSE]
  r <- reference[match(compound, reference$compound), , drop = FALSE]
  blank <- function(v) ifelse(is.na(v), "", v)
  tier <- a$EvidenceTier
  mech <- blank(a$TargetClass)
  low <- mech == "" & tier %in% c("C", "D")
  mech[low] <- ifelse(blank(a$LegacyClass1[low]) == "", "unknown", paste0("legacy: ", blank(a$LegacyClass1[low])))
  mech[mech == "" & tier %in% c("U", "B")] <- "unknown"
  mech[mech == "" & tier == "X"] <- "not a small molecule"
  out <- data.frame(compound = compound, call = .rp_call_label(d$activity_class, "chem"), activity_score = d$activity_score,
    activity_percentile = pct, padj_coherent = d$padj_coherent, nes_up_arm = d$NES_Up, nes_down_arm = d$NES_Down,
    mechanism = mech, target_genes = blank(a$TargetGenes), mechanism_family = blank(a$FamilyClass),
    mechanism_class = blank(a$BroadClass), evidence_tier = tier, evidence_confidence = blank(a$EvidenceConfidence),
    evidence_sources = blank(a$EvidenceSources), evidence_records = blank(a$EvidenceRecordIds),
    common_reference_percentile = r$reference_percentile, percentile_minus_reference = pct - r$reference_percentile,
    common_response = blank(r$common_response), stringsAsFactors = FALSE)
  out$common_hit <- (out$call == "Concordant" & out$common_response == "usually_top") |
                    (out$call == "Opposing" & out$common_response == "usually_bottom")
  out <- out[order(-out$activity_score), , drop = FALSE]; rownames(out) <- NULL
  out
}
# Mechanism classes among the called compounds: one-sided Fisher test per class and direction, BH over classes.
.rp_no_class <- "No class (tier below B)"
.deg_class_table <- function(cc, min_n = 5) {
  cls <- ifelse(cc$mechanism_class == "", .rp_no_class, cc$mechanism_class)
  N <- nrow(cc); kc <- sum(cc$call == "Concordant"); ko <- sum(cc$call == "Opposing")
  s <- lapply(split(cc, cls), function(d) c(n = nrow(d), concordant = sum(d$call == "Concordant"), opposing = sum(d$call == "Opposing"),
                                           discordant = sum(d$call == "Discordant")))
  t <- data.frame(mechanism_class = names(s), do.call(rbind, s), stringsAsFactors = FALSE, row.names = NULL)
  t$expected_concordant <- t$n * kc / N; t$expected_opposing <- t$n * ko / N
  fp <- function(k, n, K) stats::phyper(k - 1, K, N - K, n, lower.tail = FALSE)
  t$p_concordant <- fp(t$concordant, t$n, kc); t$p_opposing <- fp(t$opposing, t$n, ko)
  test <- t$n >= min_n & t$mechanism_class != .rp_no_class
  t$padj_concordant <- NA_real_; t$padj_opposing <- NA_real_
  t$padj_concordant[test] <- stats::p.adjust(t$p_concordant[test], "BH"); t$padj_opposing[test] <- stats::p.adjust(t$p_opposing[test], "BH")
  t$top_concordant <- vapply(t$mechanism_class, function(k) paste(utils::head(cc$compound[cls == k & cc$call == "Concordant"], 5), collapse = "; "), "")
  t$top_opposing <- vapply(t$mechanism_class, function(k) {
    d <- cc[cls == k & cc$call == "Opposing", ]; paste(utils::head(d$compound[order(d$activity_score)], 5), collapse = "; ") }, "")
  t <- t[order(-(t$concordant + t$opposing), t$mechanism_class), , drop = FALSE]; rownames(t) <- NULL
  t
}

# =========================================================================================================
# 4. Checks and templated text
# =========================================================================================================
.rp_hall_name <- function(s) unname(ifelse(s == "HALLMARK_E2F_TARGETS", "E2F targets", ifelse(s == "HALLMARK_G2M_CHECKPOINT", "G2M checkpoint", s)))
.rp_called <- function(d) d[d$call %in% c("Activated", "Inhibited (KO-like)", "Concordant", "Opposing"), , drop = FALSE]
.deg_checks <- function(o) {
  ct <- o$tabs$crispr; cc <- o$tabs$chem
  n_formal <- nrow(.rp_called(ct))
  disc_c <- mean(ct$call == "Discordant"); disc_x <- mean(cc$call == "Discordant")
  top_c <- utils::head(.rp_called(cc)[order(-abs(.rp_called(cc)$activity_score)), ], 20); n_common <- sum(top_c$common_hit)
  top_g <- utils::head(.rp_called(ct)[order(-abs(.rp_called(ct)$activity_score)), ], 20); n_unm <- sum(!top_g$target_measured)
  sp <- o$prolif[o$prolif$strong, , drop = FALSE]
  lv <- function(v, ok, chk) if (is.na(v)) "Check" else if (ok(v)) "OK" else if (chk(v)) "Check" else "Caution"
  cov <- o$coverage[["crispr"]]
  data.frame(
    check = c("Gene coverage", "Formal calls", "Discordant share", "Proliferation signal", "Common hits", "Unmeasured targets", "Pathways"),
    status = c(lv(cov, function(v) v >= 0.9, function(v) v >= 0.7), lv(n_formal, function(v) v >= 10, function(v) v >= 1),
               lv(max(disc_c, disc_x), function(v) v <= 0.10, function(v) v <= 0.30), if (nrow(sp)) "Check" else "OK",
               lv(n_common, function(v) v <= 3, function(v) v <= 6), lv(n_unm, function(v) v <= 2, function(v) v <= 5),
               if (o$pathway_mode == "ora") "OK" else "Check"),
    result = c(if (is.na(cov)) "not available" else sprintf("%.0f%% of the signature genes are in your table", 100 * cov),
               sprintf("%s gene knockouts and %s compounds called", .fmt_n(n_formal), .fmt_n(nrow(.rp_called(cc)))),
               sprintf("%.0f%% of gene knockouts, %.0f%% of compounds", 100 * disc_c, 100 * disc_x),
               if (nrow(sp)) paste(sprintf("%s NES %s", .rp_hall_name(sp$gene_set), .fmt_num(sp$nes, 1)), collapse = "; ") else "no strong E2F or G2M signal",
               sprintf("%d of the 20 strongest compound calls", n_common),
               sprintf("%d of the 20 strongest gene calls", n_unm),
               switch(o$pathway_mode, ora = "pathways found among the called genes", gsea = "none among the called genes; ranked test used",
                      "pathway analysis did not run")),
    meaning = c("Genes of the L1000 signatures that your table measures. Low coverage weakens every score.",
                "A call passes the score, the coherent test and both halves. With none, rankings are descriptive only.",
                "Both halves of a signature moved the same way, so no direction can be read. Common in tissue mixtures.",
                "A strong cell-cycle signal makes many top calls reflect proliferation rather than a specific target.",
                "Compounds that sit at the same end in unrelated drug experiments; their match is less specific.",
                "A call can come from the signature alone when the gene itself is absent from your table.",
                "Pathways of the called genes or, when none is found, of all genes ranked by activity score."),
    value = c(cov, n_formal, max(disc_c, disc_x), nrow(sp), n_common, n_unm, as.numeric(o$pathway_mode == "ora")),
    stringsAsFactors = FALSE)
}
.rp_top_names <- function(d, col, n = 3) {
  d <- utils::head(d, n)
  if (!nrow(d)) return("none")
  paste0(d[[col]], if ("target_measured" %in% names(d)) ifelse(d$target_measured, "", "*") else "", collapse = ", ")
}
.deg_text <- function(o) {
  ct <- o$tabs$crispr; cc <- o$tabs$chem; p <- o$prep
  act <- ct[ct$call == "Activated", ]; inh <- ct[ct$call == "Inhibited (KO-like)", ]; inh <- inh[order(inh$activity_score), ]
  conc <- cc[cc$call == "Concordant", ]; opp <- cc[cc$call == "Opposing", ]; opp <- opp[order(opp$activity_score), ]
  mech <- function(d, n = 3) { d <- utils::head(d, n); if (!nrow(d)) return("none")
    paste0(d$compound, ifelse(d$mechanism %in% c("", "unknown") | grepl("^legacy", d$mechanism), "", paste0(" (", d$mechanism, ")")), collapse = "; ") }
  star <- any(!utils::head(act$target_measured, 3)) || any(!utils::head(inh$target_measured, 3))
  sig <- o$pw[o$pw$padj < 0.05, , drop = FALSE]
  tx <- list()
  tx$input <- if (anyNA(o$coverage)) sprintf("%s genes were ranked by %s.", .fmt_n(nrow(p$rank_input)), p$stat_method) else
    sprintf("%s genes were ranked by %s. Your table covers %.0f%% of the genes in the knockout signatures and %.0f%% of those in the compound signatures.",
    .fmt_n(nrow(p$rank_input)), p$stat_method, 100 * o$coverage[["crispr"]], 100 * o$coverage[["chem"]])
  tx$crispr <- sprintf("%s of %s gene knockouts gave a call: %s Activated and %s Inhibited (KO-like). %s (%.0f%%) were Discordant.",
    .fmt_n(nrow(act) + nrow(inh)), .fmt_n(nrow(ct)), .fmt_n(nrow(act)), .fmt_n(nrow(inh)), .fmt_n(sum(ct$call == "Discordant")),
    100 * mean(ct$call == "Discordant"))
  tx$crispr_top <- sprintf("Strongest KO-like: %s. Strongest Activated: %s.%s", .rp_top_names(inh, "target"), .rp_top_names(act, "target"),
    if (star) " (* gene not in your table)" else "")
  tx$pathway <- if (o$pathway_mode == "ora")
    sprintf("%s pathways are over-represented among the called genes, led by %s.", .fmt_n(nrow(sig)), paste(utils::head(sig$term, 2), collapse = "; "))
    else if (o$pathway_mode == "none") sprintf("The pathway analysis did not run (%s).", o$pathway_error)
    else if (nrow(sig)) sprintf("No pathway is over-represented among the called genes. Ranking all %s knockout targets by activity instead, %s pathways sit at one end, led by %s (%s).",
      .fmt_n(nrow(ct)), .fmt_n(nrow(sig)), sig$term[1], sig$end_label[1])
    else "No pathway is over-represented among the called genes, and none sits at either end of the activity ranking."
  tx$chem <- sprintf("%s of %s compounds gave a call: %s Concordant (similar profile) and %s Opposing (reverse profile).",
    .fmt_n(nrow(conc) + nrow(opp)), .fmt_n(nrow(cc)), .fmt_n(nrow(conc)), .fmt_n(nrow(opp)))
  tx$chem_top <- sprintf("Top Concordant: %s. Top Opposing: %s.", mech(conc), mech(opp))
  att <- o$checks[o$checks$status != "OK", ]
  tx$attention <- if (nrow(att)) sprintf("%d of %d quality checks need a look: %s. Page %d explains each one.", nrow(att), nrow(o$checks),
    paste(tolower(att$check), collapse = ", "), .rp_pn(o, "cautions")) else "All quality checks passed."
  tx
}

# =========================================================================================================
# 5. Page furniture (pure grid; sizes in mm)
# =========================================================================================================
.rp_L <- list(page_w = 210, page_h = 297, m = 12, W = 186, H = 273, foot = 7)
.rp_status_col <- c(OK = "#2F7D4F", Check = "#A86A12", Caution = "#B5443B")
.rp_status_tint <- c(OK = "#E3F0E7", Check = "#F7ECD9", Caution = "#F4DEDB")
.rp_lh <- function(size, lh = 1.3) size * 0.3528 * lh
.rp_top <- function(mm) grid::unit(1, "npc") - grid::unit(mm, "mm")

.rp_tw <- function(s, size, face = "plain")
  grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = .gp(size, face = face))), "mm", valueOnly = TRUE)
.rp_fit <- function(s, width, size, face = "plain") {
  s <- ifelse(is.na(s), "", as.character(s))
  vapply(s, function(z) {
    if (.rp_tw(z, size, face) <= width) return(z)
    lo <- 1; hi <- nchar(z)
    while (lo < hi) { mid <- (lo + hi + 1) %/% 2; if (.rp_tw(paste0(substr(z, 1, mid), "\u2026"), size, face) <= width) lo <- mid else hi <- mid - 1 }
    paste0(trimws(substr(z, 1, lo)), "\u2026")
  }, "", USE.NAMES = FALSE)
}
.rp_wrap <- function(text, width, size, face = "plain", first = width) {
  out <- character(); lim <- first
  for (para in strsplit(as.character(text), "\n", fixed = TRUE)[[1]]) {
    words <- strsplit(para, " ", fixed = TRUE)[[1]]; line <- ""
    for (w in words) {
      cand <- if (nzchar(line)) paste(line, w) else w
      if (nzchar(line) && .rp_tw(cand, size, face) > lim) { out <- c(out, line); line <- w; lim <- width } else line <- cand
    }
    out <- c(out, line); lim <- width
  }
  if (!length(out)) "" else out
}
.rp_blk <- function(grob, h, w = .rp_L$W) list(grob = grob, h = h, w = w)
.rp_space <- function(h) .rp_blk(grid::nullGrob(), h)

.rp_text_blk <- function(text, size = 7.4, col = .tok$ink2, face = "plain", width = .rp_L$W, lh = 1.32, after = 1.6) {
  lines <- unlist(lapply(text, .rp_wrap, width = width, size = size, face = face))
  g <- grid::textGrob(paste(lines, collapse = "\n"), x = 0, y = 1, just = c("left", "top"), gp = .gp(size, col, face, lineheight = lh / 1.2))
  .rp_blk(g, length(lines) * .rp_lh(size, lh) + after, width)
}
# Bullets; an item is a string or c(lead, body) with a bold lead.
.rp_bullets <- function(items, size = 7.2, width = .rp_L$W, col = .tok$ink2, gap = 1.1, lh = 1.3, bullet = "\u2022") {
  g <- list(); y <- 0; ind <- 3.2; step <- .rp_lh(size, lh)
  for (it in items) {
    lead <- if (length(it) == 2) it[[1]] else NULL; body <- if (length(it) == 2) it[[2]] else it[[1]]
    lw <- if (is.null(lead)) 0 else .rp_tw(paste0(lead, " "), size, "bold")
    lines <- .rp_wrap(body, width - ind, size, first = width - ind - lw)
    g[[length(g) + 1]] <- grid::textGrob(bullet, x = grid::unit(0.4, "mm"), y = .rp_top(y), just = c("left", "top"), gp = .gp(size, .tok$muted))
    if (!is.null(lead)) g[[length(g) + 1]] <- grid::textGrob(lead, x = grid::unit(ind, "mm"), y = .rp_top(y), just = c("left", "top"), gp = .gp(size, .tok$ink, "bold"))
    for (k in seq_along(lines))
      g[[length(g) + 1]] <- grid::textGrob(lines[k], x = grid::unit(ind + if (k == 1) lw else 0, "mm"), y = .rp_top(y + (k - 1) * step),
        just = c("left", "top"), gp = .gp(size, col))
    y <- y + length(lines) * step + gap
  }
  .rp_blk(do.call(grid::grobTree, g), y + 0.6, width)
}
.rp_section <- function(title, note = NULL, width = .rp_L$W, size = 8.6) {
  g <- list(grid::textGrob(title, x = 0, y = .rp_top(1.2), just = c("left", "top"), gp = .gp(size, .tok$ink, "bold")),
            grid::segmentsGrob(0, .rp_top(5.6), 1, .rp_top(5.6), gp = grid::gpar(col = .tok$rule, lwd = 0.6)))
  if (!is.null(note)) g[[3]] <- grid::textGrob(note, x = 1, y = .rp_top(1.7), just = c("right", "top"), gp = .gp(6.3, .tok$muted))
  .rp_blk(do.call(grid::grobTree, g), 7.4, width)
}
# Table: styles = list(column = function(value) list(col =, face =, dot =, fill =)); wrap = columns allowed to wrap.
.rp_tbl <- function(df, widths, size = 6.4, align = NULL, wrap = character(), styles = list(), lh = 1.22, pad = 0.9,
                 header = names(df), max_lines = 3, zebra = TRUE) {
  k <- ncol(df); align <- align %||% ifelse(vapply(df, is.numeric, TRUE), "right", "left")
  x0 <- cumsum(c(0, widths))[seq_len(k)]; step <- .rp_lh(size, lh)
  cell_lines <- function(v, j, face = "plain") {
    v <- ifelse(is.na(v), "", as.character(v))
    if (names(df)[j] %in% wrap) utils::head(.rp_wrap(v, widths[j] - 2.4, size, face), max_lines) else .rp_fit(v, widths[j] - 2.4, size, face)
  }
  rows <- c(list(lapply(seq_len(k), function(j) cell_lines(header[j], j, "bold"))),
            lapply(seq_len(nrow(df)), function(i) lapply(seq_len(k), function(j) cell_lines(df[[j]][i], j))))
  hts <- vapply(rows, function(r) max(lengths(r)) * step + 2 * pad, 0)
  g <- list(); y <- 0
  for (i in seq_along(rows)) {
    fill <- if (i == 1) "#ECEBE7" else if (zebra && i %% 2 == 1) "#F7F6F3" else NA
    if (!is.na(fill)) g[[length(g) + 1]] <- grid::rectGrob(x = 0, y = .rp_top(y), width = grid::unit(sum(widths), "mm"),
      height = grid::unit(hts[i], "mm"), just = c("left", "top"), gp = grid::gpar(fill = fill, col = NA))
    for (j in seq_len(k)) {
      st <- if (i > 1 && names(df)[j] %in% names(styles)) styles[[names(df)[j]]](df[[j]][i - 1]) else list()
      txt <- rows[[i]][[j]]; xx <- x0[j] + 1.2
      if (!is.null(st$fill)) {
        fw <- min(widths[j] - 1.2, .rp_tw(txt[1], size, st$face %||% "plain") + 2.4)
        g[[length(g) + 1]] <- grid::roundrectGrob(x = grid::unit(x0[j] + 0.6, "mm"), y = .rp_top(y + pad - 0.35),
          width = grid::unit(fw, "mm"), height = grid::unit(step + 0.7, "mm"), just = c("left", "top"),
          r = grid::unit(0.6, "mm"), gp = grid::gpar(fill = st$fill, col = NA))
      }
      if (!is.null(st$dot)) {
        g[[length(g) + 1]] <- grid::circleGrob(x = grid::unit(x0[j] + 2.3, "mm"), y = .rp_top(y + pad + step * 0.45), r = grid::unit(1.1, "mm"),
          gp = grid::gpar(fill = st$dot, col = NA)); xx <- xx + 2.7
      }
      jx <- if (align[j] == "right") grid::unit(x0[j] + widths[j] - 1.2, "mm") else grid::unit(xx, "mm")
      g[[length(g) + 1]] <- grid::textGrob(paste(txt, collapse = "\n"), x = jx, y = .rp_top(y + pad), just = c(align[j], "top"),
        gp = .gp(size, st$col %||% (if (i == 1) .tok$ink2 else .tok$ink), st$face %||% (if (i == 1) "bold" else "plain"), lineheight = lh / 1.2))
    }
    y <- y + hts[i]
  }
  g[[length(g) + 1]] <- grid::segmentsGrob(0, .rp_top(hts[1]), grid::unit(sum(widths), "mm"), .rp_top(hts[1]), gp = grid::gpar(col = "#BDBBB5", lwd = 0.6))
  .rp_blk(do.call(grid::grobTree, g), y + 0.8, sum(widths))
}
# A plot as a block of its recommended size (attr size_mm) or the size given.
.rp_fig <- function(p, w = NULL, h = NULL) {
  s <- attr(p, "size_mm"); w <- w %||% s[["width"]]; h <- h %||% s[["height"]]
  g <- if (inherits(p, "patchwork")) patchwork::patchworkGrob(p) else if (inherits(p, "ggplot")) ggplot2::ggplotGrob(p) else p
  .rp_blk(g, h, w)
}
# Blocks side by side; each drawn at its own size from the top left of its slot.
.rp_row <- function(items, gap = 2, after = 1.5) {
  x <- 0; kids <- list()
  for (it in items) {
    vp <- grid::viewport(x = grid::unit(x, "mm"), y = grid::unit(1, "npc"), width = grid::unit(it$w, "mm"), height = grid::unit(it$h, "mm"),
                         just = c("left", "top"))
    kids[[length(kids) + 1]] <- grid::gTree(children = grid::gList(it$grob), vp = vp)
    x <- x + it$w + gap
  }
  .rp_blk(grid::gTree(children = do.call(grid::gList, kids)), max(vapply(items, `[[`, 0, "h")) + after, x - gap)
}
.rp_col <- function(items, width) {     # blocks stacked inside one column
  y <- 0; kids <- list()
  for (it in items) {
    vp <- grid::viewport(x = 0, y = .rp_top(y), width = grid::unit(it$w, "mm"), height = grid::unit(it$h, "mm"), just = c("left", "top"))
    kids[[length(kids) + 1]] <- grid::gTree(children = grid::gList(it$grob), vp = vp); y <- y + it$h
  }
  .rp_blk(grid::gTree(children = do.call(grid::gList, kids)), y, width)
}
# Rounded box holding blocks (build the inner blocks at width - 2 * pad, minus 1.3 with a bar).
.rp_box <- function(items, width = .rp_L$W, fill = .tok$pill, pad = 2.6, bar = NULL) {
  off <- if (is.null(bar)) 0 else 1.3
  inner <- .rp_col(items, width - 2 * pad - off)
  h <- inner$h + 2 * pad - 1
  vp <- grid::viewport(x = grid::unit(pad + off, "mm"), y = .rp_top(pad), width = grid::unit(inner$w, "mm"),
                       height = grid::unit(inner$h, "mm"), just = c("left", "top"))
  g <- list(grid::roundrectGrob(x = 0, y = 1, just = c("left", "top"), width = grid::unit(width, "mm"), height = grid::unit(h, "mm"),
                                r = grid::unit(1.2, "mm"), gp = grid::gpar(fill = fill, col = NA)))
  if (!is.null(bar)) g[[2]] <- grid::rectGrob(x = 0, y = 1, just = c("left", "top"), width = grid::unit(1.3, "mm"), height = grid::unit(h, "mm"),
                                               gp = grid::gpar(fill = bar, col = NA))
  g[[length(g) + 1]] <- grid::gTree(children = grid::gList(inner$grob), vp = vp)
  .rp_blk(do.call(grid::grobTree, g), h + 2, width)
}
.rp_inner <- function(width, bar = FALSE) width - 5.2 - if (bar) 1.3 else 0
.rp_box_title <- function(t, width, size = 6.3, col = .tok$muted)
  .rp_blk(grid::textGrob(toupper(t), x = 0, y = .rp_top(0.6), just = c("left", "top"), gp = .gp(size, col, "bold")), 4.4, width)

# Stat tiles in labelled groups: groups = list(list(title, tiles = list(list(label, value, cls)))).
.rp_tiles <- function(groups, width = .rp_L$W, gap = 2, h = 15) {
  n <- sum(vapply(groups, function(g) length(g$tiles), 0)); tw <- (width - gap * (n - 1)) / n
  g <- list(); x <- 0
  for (grp in groups) {
    g[[length(g) + 1]] <- grid::textGrob(toupper(grp$title), x = grid::unit(x, "mm"), y = .rp_top(0.5), just = c("left", "top"), gp = .gp(6, .tok$muted, "bold"))
    for (t in grp$tiles) {
      cls <- t$cls
      fill <- if (is.null(cls)) .tok$pill else .pal$tint[[cls]]; tc <- if (is.null(cls)) .tok$ink else .pal$text[[cls]]
      g[[length(g) + 1]] <- grid::roundrectGrob(x = grid::unit(x, "mm"), y = .rp_top(4.2), width = grid::unit(tw, "mm"), height = grid::unit(h, "mm"),
        just = c("left", "top"), r = grid::unit(1, "mm"), gp = grid::gpar(fill = fill, col = NA))
      g[[length(g) + 1]] <- grid::textGrob(t$label, x = grid::unit(x + 2.3, "mm"), y = .rp_top(7.3), just = c("left", "centre"),
        gp = .gp(6.3, if (is.null(cls)) .tok$ink2 else tc, if (is.null(cls)) "plain" else "bold"))
      g[[length(g) + 1]] <- grid::textGrob(t$value, x = grid::unit(x + 2.3, "mm"), y = .rp_top(4.2 + h - 2.4), just = c("left", "bottom"),
        gp = .gp(14, tc, "bold"))
      x <- x + tw + gap
    }
  }
  .rp_blk(do.call(grid::grobTree, g), h + 6.5, width)
}
.rp_status_style <- function(v) list(col = .rp_status_col[[v]], face = "bold", dot = .rp_status_col[[v]])
.rp_call_key <- c(Activated = "Activated", `Inhibited (KO-like)` = "Inhibited", Concordant = "Activated", Opposing = "Inhibited", Discordant = "Discordant",
               `Activated end` = "Activated", `KO-like end` = "Inhibited", Inhibited = "Inhibited")
.rp_call_style <- function(v) {
  if (!v %in% names(.rp_call_key)) return(list(col = .tok$muted))
  list(col = .pal$text[[.rp_call_key[[v]]]], face = "bold", fill = .pal$tint[[.rp_call_key[[v]]]])
}

# Running head, page title (a question) and its data-specific answer.
.rp_page_head <- function(o, key, answer) {
  num <- .rp_pn(o, key); question <- .rp_specs[[key]]$q
  run <- .rp_blk(grid::grobTree(
    grid::textGrob(sprintf("L1000PAI report  \u00b7  %s", .rp_fit(o$meta$contrast, 150, 6.4)), x = 0, y = .rp_top(0.4), just = c("left", "top"), gp = .gp(6.4, .tok$muted)),
    grid::textGrob(sprintf("%d / %d", num, o$meta$n_pages), x = 1, y = .rp_top(0.4), just = c("right", "top"), gp = .gp(6.4, .tok$muted)),
    grid::segmentsGrob(0, .rp_top(4), 1, .rp_top(4), gp = grid::gpar(col = .tok$rule, lwd = 0.6))), 6.5)
  ttl <- .rp_blk(grid::grobTree(
    grid::textGrob(as.character(num), x = 0, y = .rp_top(0.5), just = c("left", "top"), gp = .gp(13, .pal$marks[["Inhibited"]], "bold")),
    grid::textGrob(question, x = grid::unit(6.5, "mm"), y = .rp_top(0.5), just = c("left", "top"), gp = .gp(13, .tok$ink, "bold"))), 8.4)
  list(run, ttl, .rp_text_blk(answer, size = 7.9, col = .tok$ink2, after = 2.6))
}
# Draw a built page on the current (already started) page.
.rp_draw_page <- function(blocks, o, num) {
  grid::pushViewport(grid::viewport(x = grid::unit(.rp_L$m, "mm"), y = grid::unit(.rp_L$m, "mm"), width = grid::unit(.rp_L$W, "mm"),
                                    height = grid::unit(.rp_L$H, "mm"), just = c("left", "bottom")))
  y <- 0
  for (b in blocks) {
    grid::pushViewport(grid::viewport(x = 0, y = .rp_top(y), width = grid::unit(b$w, "mm"), height = grid::unit(b$h, "mm"), just = c("left", "top")))
    grid::grid.draw(b$grob); grid::popViewport(); y <- y + b$h
  }
  grid::grid.segments(0, grid::unit(4.2, "mm"), 1, grid::unit(4.2, "mm"), gp = grid::gpar(col = .tok$rule, lwd = 0.6))
  grid::grid.text(o$meta$footer, x = 0, y = grid::unit(1.6, "mm"), just = c("left", "centre"), gp = .gp(5.7, .tok$muted))
  grid::popViewport()
  invisible(y)
}

# =========================================================================================================
# 6. Plots drawn here (input histogram, proliferation bars, mechanism classes)
# =========================================================================================================
.rp_plot_stat_hist <- function(prep, w = 92, h = 52) {
  d <- prep$rank_input
  lim <- stats::quantile(abs(d$stat), 0.999, na.rm = TRUE)
  d$side <- ifelse(d$stat > 0, "up", "down"); d$v <- pmax(pmin(d$stat, lim), -lim)
  p <- ggplot2::ggplot(d, ggplot2::aes(v, fill = side)) +
    ggplot2::geom_histogram(ggplot2::aes(colour = side), bins = 60, boundary = 0, linewidth = 0.15) +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c(up = .pal$marks[["Activated"]], down = .pal$marks[["Inhibited"]]), guide = "none") +
    ggplot2::scale_colour_manual(values = c(up = .pal$marks[["Activated"]], down = .pal$marks[["Inhibited"]]), guide = "none") +
    ggplot2::scale_x_continuous(labels = .lab_minus) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.05))) +
    ggplot2::labs(x = sprintf("Ranking statistic (%s)", prep$stat_col), y = "Genes", title = "Ranking statistic",
      subtitle = sprintf("%s genes up (> 0), %s down (< 0)", .fmt_n(sum(d$stat > 0)), .fmt_n(sum(d$stat < 0)))) +
    theme_pai(grid = "y")
  .with_size(p, w, h)
}
.rp_plot_prolif <- function(pr, w = 92, h = 36) {
  d <- pr; d$label <- .rp_hall_name(d$gene_set)
  d$txt <- sprintf("NES %s,  adj. P %s", .fmt_num(d$nes, 2), .rp_fmt_padj(d$padj))
  lim <- max(3, ceiling(max(abs(d$nes))))
  p <- ggplot2::ggplot(d, ggplot2::aes(nes, label)) +
    ggplot2::annotate("rect", xmin = -2, xmax = 2, ymin = -Inf, ymax = Inf, fill = "#F1F0EC") +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::geom_col(ggplot2::aes(fill = nes > 0), width = 0.55) +
    ggplot2::geom_text(ggplot2::aes(x = ifelse(nes > 0, nes + 0.15, nes - 0.15), label = txt, hjust = ifelse(nes > 0, 0, 1)),
      size = 2.1, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = .pal$marks[["Activated"]], `FALSE` = .pal$marks[["Inhibited"]]), guide = "none") +
    ggplot2::scale_x_continuous(limits = c(-lim - 3.4, lim + 3.4), breaks = seq(-lim, lim, 1), labels = .lab_minus) +
    ggplot2::labs(x = "NES in your ranking (grey band: |NES| < 2)", y = NULL, title = "Proliferation signal",
                  subtitle = "MSigDB Hallmark sets; flagged at adj. P < 0.05 and |NES| \u2265 2") +
    theme_pai(grid = "none")
  .with_size(p, w, h)
}
.rp_plot_classes <- function(kt, n = 10, w = 186) {
  d <- kt[kt$mechanism_class != .rp_no_class & kt$n >= 5 & (kt$concordant + kt$opposing) > 0, ]
  d <- utils::head(d, n)
  if (!nrow(d)) return(NULL)
  d$lab <- sprintf("%s  (%s)", d$mechanism_class, .fmt_n(d$n))
  d$lab <- factor(d$lab, levels = rev(d$lab))
  long <- rbind(data.frame(lab = d$lab, side = "Opposing", k = -d$opposing, e = -d$expected_opposing, q = d$padj_opposing),
                data.frame(lab = d$lab, side = "Concordant", k = d$concordant, e = d$expected_concordant, q = d$padj_concordant))
  long$txt <- ifelse(long$k == 0, "", paste0(abs(long$k), ifelse(!is.na(long$q) & long$q < 0.05, " *", "")))
  m <- max(abs(c(long$k, long$e))) * 1.25
  p <- ggplot2::ggplot(long, ggplot2::aes(k, lab)) +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$ink2, linewidth = 0.3) +
    ggplot2::geom_col(ggplot2::aes(fill = side), width = 0.62) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = e, xmax = e), width = 0.8, colour = .tok$ink, linewidth = 0.45) +
    ggplot2::geom_text(ggplot2::aes(x = ifelse(k >= 0, pmax(k, e) + m * 0.02, pmin(k, e) - m * 0.02), label = txt,
      hjust = ifelse(k >= 0, 0, 1)), size = 2.2, colour = .tok$ink, family = .pai_font()) +
    ggplot2::scale_fill_manual(values = c(Opposing = .pal$marks[["Inhibited"]], Concordant = .pal$marks[["Activated"]]),
      breaks = c("Opposing", "Concordant"), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(-m, m), labels = function(v) .lab_minus(abs(v))) +
    ggplot2::labs(x = "Called compounds   \u2190 Opposing  |  Concordant \u2192", y = NULL,
      title = "Mechanism classes among the called compounds",
      subtitle = "Bars: called compounds of the class (class size in brackets). Black tick: number expected by chance. * adj. P < 0.05 (one-sided Fisher, BH)") +
    theme_pai(grid = "x") + ggplot2::theme(legend.position = "top", legend.justification = "left")
  .with_size(p, w, 27 + 5 * nrow(d))
}

# =========================================================================================================
# 7. Pages
# =========================================================================================================
# Pages of the report: question (page title), what the page holds (contents table), page function
# and the number of table rows it starts with (rows are dropped until the page fits).
.rp_specs <- list(
  overview = list(q = "Overview", what = "Headline numbers, key findings, quality checks, how to read the calls",
                  fun = ".rp_page_overview"),
  input = list(q = "Was your table read correctly?", what = "Input summary, statistic distribution, proliferation check, top genes",
               fun = ".rp_page_input"),
  crispr = list(q = "Which gene knockouts does your contrast resemble?",
                what = "All knockouts, the strongest calls, and the gene's own change in your data", fun = ".rp_page_crispr", n = 10),
  evidence = list(q = "How strong is the evidence for the top gene calls?",
                  what = "Both halves of each top knockout signature; the genes that drive them", fun = ".rp_page_evidence"),
  pathways = list(q = "Which pathways do the gene calls point to?", what = "Pathways of the gene calls, or along the activity ranking",
                  fun = ".rp_page_pathways", n = 12),
  enrich_overview = list(q = "Which collections and directions are enriched?",
                         what = "Strongest terms of each collection for Activated and Inhibited genes", fun = ".rp_page_enrich_overview"),
  enrich_terms = list(q = "Which terms stand out, and which genes drive them?", what = "Every tested term by direction; top terms with their genes",
                      fun = ".rp_page_enrich_terms", n = 12),
  enrich_circle = list(q = "How do the top terms and their genes fit together?", what = "The ten top terms with all their member genes",
                       fun = ".rp_page_enrich_circle"),
  enrich_chord = list(q = "Which genes do the top terms share?", what = "Genes shared by the top terms: chord and Venn diagrams",
                      fun = ".rp_page_enrich_chord"),
  enrich_bars = list(q = "How strong is each top term, and which of its genes belong together?",
                     what = "Top terms as bars; their genes clustered by shared terms", fun = ".rp_page_enrich_bars", n = 15, n_min = 8),
  enrich_heat = list(q = "How active is each gene of the top terms?", what = "Heat map of the genes of the top terms by activity score",
                     fun = ".rp_page_enrich_heat", n = 30, n_min = 16),
  enrich_armmap = list(q = "Where do the genes of each top term sit on the two halves?", what = "Both halves of every gene of the top terms, term by term",
                       fun = ".rp_page_enrich_armmap"),
  enrich_hubs = list(q = "Which genes link the enriched pathways?", what = "Pathway modules and the hub genes found in many of them",
                     fun = ".rp_page_enrich_hubs", n = 15),
  enrich_hubnet = list(q = "How are the hub genes wired to the modules?", what = "Hub genes and modules as a chord diagram and a network",
                       fun = ".rp_page_enrich_hubnet"),
  chem = list(q = "Which compounds produce a similar or opposite profile?", what = "All compounds; the strongest with mechanism, evidence, common hits",
              fun = ".rp_page_chem", n = 8),
  chem_evidence = list(q = "How strong is the evidence for the top compound calls?",
                       what = "Both halves of each top compound signature; the genes that drive them", fun = ".rp_page_chem_evidence"),
  classes = list(q = "Do compounds with the same mechanism agree?", what = "Curated drug families and mechanism classes among the calls",
                 fun = ".rp_page_classes"),
  cautions = list(q = "What should you check before drawing conclusions?", what = "Flags for this contrast, notes for every report, next steps",
                  fun = ".rp_page_cautions"),
  methods = list(q = "How was this report made?", what = "Methods, settings, versions, files and references", fun = ".rp_page_methods"))
# The pathway-enrichment pages depend on what the pathway analysis found.
.rp_page_keys <- function(mode, enrichment = "full") {
  enr <- if (enrichment == "none") character() else switch(mode,
    ora = if (enrichment == "core") c("enrich_overview", "enrich_terms", "enrich_hubs") else
      c("enrich_overview", "enrich_terms", "enrich_circle", "enrich_chord", "enrich_bars", "enrich_heat", "enrich_armmap",
        "enrich_hubs", "enrich_hubnet"),
    gsea = if (enrichment == "core") "enrich_hubs" else c("enrich_hubs", "enrich_hubnet"), character())
  c("overview", "input", "crispr", "evidence", "pathways", enr, "chem", "chem_evidence", "classes", "cautions", "methods")
}
.rp_pn <- function(o, key) match(key, o$page_keys)
.rp_pct <- function(v, d = 1) ifelse(is.na(v), "not available", sprintf(paste0("%.", d, "f%%"), 100 * v))
.rp_fmt_padj <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, .fmt_p(p), sprintf("%.3f", p)))
.rp_sgn <- function(v, d = 2) ifelse(is.na(v), "", .fmt_num(v, d))

.rp_page_overview <- function(o) {
  ct <- o$tabs$crispr; cc <- o$tabs$chem; tx <- o$text
  n <- function(d, k) .fmt_n(sum(d$call == k))
  list(
    .rp_title_blk(o),
    .rp_tiles(list(
      list(title = "Input", tiles = list(list(label = "Genes ranked", value = .fmt_n(nrow(o$prep$rank_input))))),
      list(title = "Gene knockouts (CRISPR)", tiles = list(
        list(label = "Activated", value = n(ct, "Activated"), cls = "Activated"),
        list(label = "Inhibited (KO-like)", value = n(ct, "Inhibited (KO-like)"), cls = "Inhibited"),
        list(label = "Discordant", value = n(ct, "Discordant"), cls = "Discordant"))),
      list(title = "Compounds (Chem)", tiles = list(
        list(label = "Concordant", value = n(cc, "Concordant"), cls = "Activated"),
        list(label = "Opposing", value = n(cc, "Opposing"), cls = "Inhibited"),
        list(label = "Discordant", value = n(cc, "Discordant"), cls = "Discordant"))))),
    .rp_section("Key findings"),
    .rp_bullets(list(c("Input.", tx$input), c("Gene knockouts.", paste(tx$crispr, tx$crispr_top)), c("Pathways.", tx$pathway),
                  c("Compounds.", paste(tx$chem, tx$chem_top)), c("Attention.", tx$attention)), size = 7.3),
    .rp_space(1),
    .rp_section("Quality checks", sprintf("OK  \u00b7  Check: read page %d  \u00b7  Caution: interpret with care", .rp_pn(o, "cautions"))),
    .rp_tbl(data.frame(Check = o$checks$check, Status = o$checks$status, Result = o$checks$result, `What it means` = o$checks$meaning,
                    check.names = FALSE), c(28, 17, 54, 87), size = 6.5, wrap = c("Result", "What it means"),
         styles = list(Status = .rp_status_style)),
    .rp_space(2),
    .rp_section("How to read the calls"),
    .rp_glossary(o),
    .rp_space(1),
    .rp_section("In this report", if ("xlsx" %in% names(o$meta$files)) sprintf("companion workbook: %s", basename(o$meta$files[["xlsx"]]))),
    .rp_tbl(.rp_contents_table(o), c(12, 76, 98), size = 6.4, align = c("left", "left", "left"), wrap = "What you will find"))
}
# Contents of the report; the pathway-enrichment pages share one row when there are more than one.
.rp_contents_table <- function(o) {
  k <- o$page_keys; pg <- seq_along(k); enr <- grepl("^enrich_", k)
  rows <- data.frame(Page = as.character(pg), Question = vapply(k, function(z) .rp_specs[[z]]$q, ""),
                     `What you will find` = vapply(k, function(z) .rp_specs[[z]]$what, ""), check.names = FALSE, stringsAsFactors = FALSE)
  if (sum(enr) > 1) {
    i <- which(enr)
    what <- c(enrich_overview = "collections", enrich_terms = "top terms and genes", enrich_circle = "term circle",
              enrich_chord = "chord and Venn", enrich_bars = "bars and gene clusters", enrich_heat = "heat map", enrich_armmap = "arm map",
              enrich_hubs = "hub genes and modules", enrich_hubnet = "hub network")
    rows <- rbind(rows[seq_len(min(i) - 1), ],
                  data.frame(Page = sprintf("%d\u2013%d", min(i), max(i)), Question = "Pathway enrichment in depth",
                             `What you will find` = paste(what[k[i]], collapse = ", "), check.names = FALSE, stringsAsFactors = FALSE),
                  rows[(max(i) + 1):nrow(rows), ])
  }
  rows[-1, , drop = FALSE]
}
.rp_title_blk <- function(o) {
  cl <- .rp_wrap(o$meta$contrast, .rp_L$W, 10.5); y2 <- 9.6 + length(cl) * .rp_lh(10.5, 1.25) + 0.8
  .rp_blk(grid::grobTree(
    grid::textGrob("L1000PAI report", x = 0, y = .rp_top(1), just = c("left", "top"), gp = .gp(18, .tok$ink, "bold")),
    grid::textGrob(paste(cl, collapse = "\n"), x = 0, y = .rp_top(9.6), just = c("left", "top"), gp = .gp(10.5, .tok$ink2, lineheight = 1.25 / 1.2)),
    grid::textGrob(o$meta$byline, x = 0, y = .rp_top(y2), just = c("left", "top"), gp = .gp(6.6, .tok$muted))), y2 + 5.8)
}
.rp_glossary <- function(o, size = 6.8) {
  sc <- .pai_semantics("crispr", .pai_score_direction(o$x, "crispr")); sx <- .pai_semantics("chem", .pai_score_direction(o$xc, "chem"))
  rows <- list(
    list("Activated", "Activated", sprintf("Gene knockout: %s. The gene behaves as if more active.", sc$explain[["Activated"]])),
    list("Inhibited", "Inhibited (KO-like)", sprintf("Gene knockout: %s. The gene behaves as if lost or inhibited.", sc$explain[["Inhibited"]])),
    list("Activated", "Concordant", sprintf("Compound: %s. The compound produces a similar profile in cell lines.", sx$explain[["Activated"]])),
    list("Inhibited", "Opposing", sprintf("Compound: %s. The compound produces the reverse profile; a candidate to test, not a treatment.", sx$explain[["Inhibited"]])),
    list("Discordant", "Discordant", "The up-genes and the down-genes of the signature both moved the same way, so no direction can be read."),
    list("Not_significant", "Not called", "Below the call thresholds: |activity score| \u2265 1, coherent adj. P < 0.05 and both halves aligned."))
  g <- list(); y <- 0; step <- .rp_lh(size, 1.3)
  for (r in rows) {
    lines <- .rp_wrap(r[[3]], .rp_L$W - 32, size)
    g[[length(g) + 1]] <- grid::circleGrob(x = grid::unit(1.3, "mm"), y = .rp_top(y + step * 0.45), r = grid::unit(1.2, "mm"),
      gp = grid::gpar(fill = .pal$marks[[r[[1]]]], col = NA))
    g[[length(g) + 1]] <- grid::textGrob(r[[2]], x = grid::unit(4, "mm"), y = .rp_top(y), just = c("left", "top"), gp = .gp(size, .tok$ink, "bold"))
    g[[length(g) + 1]] <- grid::textGrob(paste(lines, collapse = "\n"), x = grid::unit(32, "mm"), y = .rp_top(y), just = c("left", "top"),
      gp = .gp(size, .tok$ink2, lineheight = 1.3 / 1.2))
    y <- y + length(lines) * step + 1.1
  }
  .rp_blk(do.call(grid::grobTree, g), y + 1)
}

.rp_page_input <- function(o) {
  p <- o$prep; c_ <- p$counts
  kv <- data.frame(Item = c("Rows in your table", "Gene column", "Identifier type", "Ranking statistic", "Rows without a gene symbol",
                            "Rows without a statistic", "Duplicate genes collapsed", "Genes used", "Coverage, knockout signatures",
                            "Coverage, compound signatures"),
                   Value = c(.fmt_n(c_[["rows"]]), p$gene_col, p$id_type, p$stat_method, .fmt_n(c_[["no_symbol"]]), .fmt_n(c_[["no_stat"]]),
                             sprintf("%s (largest |statistic| kept)", .fmt_n(c_[["duplicates"]])), .fmt_n(c_[["used"]]),
                             .rp_pct(o$coverage[["crispr"]]), .rp_pct(o$coverage[["chem"]])),
                   stringsAsFactors = FALSE)
  top_tab <- function(dir) {
    q <- p$qc[p$qc$status == "used", , drop = FALSE]; q <- q[order(if (dir == "up") -q$statistic else q$statistic), , drop = FALSE]
    q <- utils::head(q, 12)
    d <- data.frame(Gene = q$gene, Statistic = .rp_sgn(q$statistic), stringsAsFactors = FALSE)
    if (!is.null(q$log_fold_change)) d$`Fold change (log)` <- .rp_sgn(q$log_fold_change, 3)
    if (!is.null(q$adj_p)) d$`adj. P` <- .rp_fmt_padj(q$adj_p)
    d
  }
  wd <- function(d) { w <- c(24, 18, 26, 24)[seq_len(ncol(d))]; w * 91 / sum(w) }
  up <- top_tab("up"); dn <- top_tab("down")
  pr <- o$prolif; strong <- pr[pr$strong, , drop = FALSE]
  ptxt <- if (nrow(strong)) sprintf(paste("Your ranking has a strong proliferation component (%s). Many L1000 signatures share this axis,",
    "so knockouts and compounds that change cell growth can rank high for this reason alone. Treat top hits that are cell-cycle genes",
    "or antiproliferative compounds with care (pages %d and %d)."),
    paste(sprintf("%s NES %s", .rp_hall_name(strong$gene_set), .fmt_num(strong$nes, 1)), collapse = ", "), .rp_pn(o, "crispr"), .rp_pn(o, "chem")) else
    "No strong proliferation component: the E2F targets and G2M checkpoint gene sets are not concentrated at either end of your ranking."
  answer <- sprintf("%s of %s rows were used, ranked by %s. %s", .fmt_n(c_[["used"]]), .fmt_n(c_[["rows"]]), p$stat_method,
    if (nrow(strong)) "The ranking has a strong proliferation component; see below." else "No strong proliferation component was found.")
  bar <- nrow(strong) > 0
  c(.rp_page_head(o, "input", answer),
    list(.rp_row(list(.rp_col(list(.rp_section("What was read", width = 92), .rp_tbl(kv, c(44, 48), size = 6.5, align = c("left", "left"))), 92),
                   .rp_fig(.rp_plot_stat_hist(p))), gap = 2),
         .rp_space(1),
         .rp_row(list(if (all(is.na(pr$nes))) .rp_text_blk("Proliferation signal: not tested, too few genes of the gene sets in your table.", size = 7, width = 92)
                      else .rp_fig(.rp_plot_prolif(pr)),
                   .rp_col(list(.rp_space(3), .rp_box(list(.rp_text_blk(ptxt, size = 6.9, width = .rp_inner(92, bar), after = 0.4)), width = 92,
                     fill = if (bar) .rp_status_tint[["Check"]] else .tok$pill, bar = if (bar) .rp_status_col[["Check"]] else NULL)), 92)), gap = 2),
         .rp_space(1),
         .rp_row(list(.rp_col(list(.rp_section("Most increased genes", width = 91), .rp_tbl(up, wd(up), size = 6.4)), 91),
                   .rp_col(list(.rp_section("Most decreased genes", width = 91), .rp_tbl(dn, wd(dn), size = 6.4)), 91)), gap = 4)))
}

.rp_target_rows <- function(d) {
  own <- ifelse(!d$target_measured, "not measured",
    sprintf("%s %s", ifelse(d$target_statistic > 0, "\u2191", "\u2193"), .fmt_num(d$target_statistic, 2)))
  data.frame(Gene = d$target, Score = .rp_sgn(d$activity_score), `adj. P` = .rp_fmt_padj(d$padj_coherent), `Gene in your data` = own,
             check.names = FALSE, stringsAsFactors = FALSE)
}
.rp_page_crispr <- function(o, n = 10) {
  ct <- o$tabs$crispr; tx <- o$text
  act <- ct[ct$call == "Activated", ]; inh <- ct[ct$call == "Inhibited (KO-like)", ]; inh <- inh[order(inh$activity_score), ]
  unm <- function(v) if (identical(v, "not measured")) list(col = .tok$muted) else list()
  tw <- c(24, 15, 18, 34) * 91 / 91
  c(.rp_page_head(o, "crispr", paste(tx$crispr, tx$crispr_top)),
    list(.rp_row(list(.rp_cached_fig(o, "arm_crispr", pai_armplot(o$x, n_label = 6)), .rp_cached_fig(o, "dot_crispr", pai_dotplot(o$x, n = 10))), gap = 2),
         .rp_row(list(.rp_col(list(.rp_section("Strongest Activated", sprintf("%s in total", .fmt_n(nrow(act))), width = 91),
                             .rp_tbl(.rp_target_rows(utils::head(act, n)), tw, styles = list(`Gene in your data` = unm))), 91),
                   .rp_col(list(.rp_section("Strongest Inhibited (KO-like)", sprintf("%s in total", .fmt_n(nrow(inh))), width = 91),
                             .rp_tbl(.rp_target_rows(utils::head(inh, n)), tw, styles = list(`Gene in your data` = unm))), 91)), gap = 4),
         .rp_text_blk(paste("Score: activity score (positive = Activated). adj. P: coherent test, BH-adjusted. Gene in your data: the knocked-out",
           "gene's own statistic in your table (\u2191 up, \u2193 down). A KO-like call whose gene is itself down, or an Activated call whose gene is up,",
           "is more consistent; 'not measured' means the call rests on the signature alone. Every target is in the Excel sheet CRISPR_targets."),
           size = 6.4, col = .tok$muted)))
}
# An evidence card at most the content width (the card widens its name column for long names).
.rp_card <- function(o, x, key) .rp_cached_fig(o, key, {
  p <- pai_evidence_card(x, n = 20, width = .rp_L$W)
  w <- attr(p, "size_mm")[["width"]]
  if (w > .rp_L$W) pai_evidence_card(x, n = 20, width = .rp_L$W - (w - .rp_L$W)) else p
})
# Card, then the leading-edge network (pai_cnetplot) beside a column of reading notes.
.rp_evidence_blocks <- function(o, x, lib, key) {
  chem <- lib == "chem"
  dirs <- if (chem) c("concordance", "opposition") else c("activation", "inhibition")
  card_notes <- list(
    c("arms.", sprintf(paste("Upper bar: the up half of the signature, lower bar: the down half, both flipped so that right points toward",
      "%s. Two bars on the same side make a coherent call."), dirs[1])),
    c("NES Up, NES Down.", "Enrichment of each half in your ranking (+: among your increased genes, \u2212: decreased)."),
    c("score.", if (chem) "NES Up \u2212 NES Down: + Concordant (similar profile), \u2212 Opposing (reverse profile)." else
      "NES Down \u2212 NES Up: + Activated (reverse of the knockout), \u2212 Inhibited (KO-like)."),
    c("p act, p inh, q coh.", sprintf("One-sided P for %s and for %s; BH-adjusted coherent P that decides the call (< 0.05).", dirs[1], dirs[2])))
  net_notes <- list(
    c("Squares:", sprintf("the two strongest %s on each side.", if (chem) "compounds" else "knockouts")),
    c("Dots:", "up to six leading-edge genes per square, coloured by their own change in your table (red up, blue down), sized by |statistic|; bold genes are shared."),
    c("Lines:", "solid through the up half of the signature, dashed through the down half."))
  net <- .rp_cached_fig(o, paste0("cnet_", lib), pai_cnetplot(x, n = 4))
  sw <- .rp_L$W - net$w - 2; iw <- .rp_inner(sw)
  side <- .rp_col(list(.rp_space(1),
    .rp_box(list(.rp_box_title("Reading the card", iw), .rp_bullets(card_notes, size = 6.5, width = iw, gap = 0.8)), width = sw), .rp_space(1),
    .rp_box(list(.rp_box_title("Reading the network", iw), .rp_bullets(net_notes, size = 6.5, width = iw, gap = 0.8)), width = sw)), sw)
  list(.rp_card(o, x, paste0("card_", lib)), .rp_space(2),
       .rp_section(sprintf("Which of your genes carry the strongest %s calls?", if (chem) "compound" else "knockout")),
       .rp_row(list(net, side), gap = 2))
}
.rp_page_evidence <- function(o) {
  called <- .rp_called(o$tabs$crispr)
  answer <- sprintf(paste("Each row shows the two halves of a knockout signature: the genes that go up and the genes that go down after the",
    "knockout. A call is strongest when both halves point the same way with small P values. The %d strongest of the %s calls are shown;",
    "below, the genes of your table that drive the four strongest."), min(20, nrow(called)), .fmt_n(nrow(called)))
  c(.rp_page_head(o, "evidence", answer), .rp_evidence_blocks(o, o$x, "crispr", "card"))
}
.rp_page_chem_evidence <- function(o) {
  called <- .rp_called(o$tabs$chem)
  answer <- sprintf(paste("Each row shows the two halves of a compound signature: the genes that go up and the genes that go down after",
    "treatment. A call is strongest when both halves point the same way with small P values. The %d strongest of the %s calls are",
    "shown; below, the genes of your table that drive the four strongest."), min(20, nrow(called)), .fmt_n(nrow(called)))
  c(.rp_page_head(o, "chem_evidence", answer), .rp_evidence_blocks(o, o$xc, "chem", "card_chem"))
}
.rp_page_pathways <- function(o, n = 12) {
  pw <- o$pw; tx <- o$text
  if (o$pathway_mode == "none")
    return(c(.rp_page_head(o, "pathways", tx$pathway),
             list(.rp_box(list(.rp_text_blk(paste("Pathway gene sets could not be read or built, so this page is empty. The other pages are",
               "not affected. Run pai_enrich() on the CRISPR result in the R objects file to see the error, for example a missing",
               "network connection on first use."), size = 7, width = .rp_inner(.rp_L$W), after = 0.4))))))
  sig <- pw[pw$padj < 0.05, , drop = FALSE]
  tab <- utils::head(sig, n)
  if (o$pathway_mode == "ora") {
    fig <- .rp_cached_fig(o, "ladder", pai_ladder(o$e, n = 10))
    td <- data.frame(Collection = tab$collection, Pathway = tab$term,
                     Direction = ifelse(tab$direction == "Inhibited", "Inhibited (KO-like)", tab$direction), Genes = tab$count,
                     `adj. P` = .rp_fmt_padj(tab$padj), check.names = FALSE)
    tw <- c(20, 100, 30, 14, 22)
    why <- paste("The Activated and the Inhibited (KO-like) genes were tested separately for over-represented pathways; the",
                 "background is every gene the knockout library scored.")
    read <- paste("Each row is a pathway over-represented among the calls. The counts give its Activated and Inhibited members; the",
                  "points show where each member's activity score falls. Pathways summarise the gene calls; they are not separate evidence.")
  } else {
    fig <- .rp_cached_fig(o, "pathwayplot", pai_pathwayplot(o$g, ora = o$e, n = 4))
    td <- data.frame(Collection = tab$collection, Pathway = tab$term, `Gene end` = tab$end_label, NES = .rp_sgn(tab$nes),
                     Genes = tab$size, `adj. P` = .rp_fmt_padj(tab$padj), check.names = FALSE)
    tw <- c(20, 98, 24, 12, 12, 20)
    why <- sprintf(paste("No pathway was over-represented among the %s called genes (smallest adj. P %.2f). The ranked test uses every scored",
      "gene instead: all %s knockout targets are ordered by activity score, and each pathway is tested for concentration at the",
      "Activated end (NES > 0) or the KO-like end (NES < 0)."), .fmt_n(nrow(.rp_called(o$tabs$crispr))),
      min(.as_tt(o$e)$table$p.adjust), .fmt_n(nrow(o$tabs$crispr)))
    read <- paste("Each pathway is a group of knockout targets. A pathway whose genes cluster at one end suggests that the process",
                  "as a whole looks activated or lost in your contrast. Pathways summarise the gene calls; they are not separate evidence.")
  }
  boxes <- function(w) list(
    .rp_box(list(.rp_box_title("Why this view", .rp_inner(w)), .rp_text_blk(why, size = 6.8, width = .rp_inner(w), after = 0.4)), width = w),
    .rp_box(list(.rp_box_title("Reading the plot", .rp_inner(w)), .rp_text_blk(read, size = 6.8, width = .rp_inner(w), after = 0.4)), width = w))
  sw <- .rp_L$W - fig$w - 2
  top <- if (sw >= 50) { b <- boxes(sw); list(.rp_row(list(fig, .rp_col(list(.rp_space(1), b[[1]], .rp_space(1), b[[2]]), sw)), gap = 2)) } else
    list(fig, .rp_row(boxes((.rp_L$W - 3) / 2), gap = 3))
  c(.rp_page_head(o, "pathways", tx$pathway), top,
    list(.rp_space(1),
         .rp_section("Top pathways", sprintf("%s with adj. P < 0.05; all in the Excel sheet Pathways", .fmt_n(nrow(sig)))),
         .rp_tbl(td, tw, wrap = "Pathway", max_lines = 2, styles = list(`Gene end` = .rp_call_style, Direction = .rp_call_style))))
}
# ---- pathway enrichment in depth (pai_enrich result; hubs also from pai_gsea_pathways) -------------------
.rp_hub_label <- function(o, dir) if (identical(o$pathway_mode, "gsea"))
  c(Inhibited = "KO-like end of the ranking", Activated = "Activated end of the ranking")[[dir]] else
  c(Inhibited = "Inhibited (KO-like) genes", Activated = "Activated genes")[[dir]]
.rp_page_enrich_overview <- function(o) {
  t <- o$e$table; sig <- t[t$p.adjust < 0.05, , drop = FALSE]
  cols <- unique(t$collection)
  per <- data.frame(Collection = cols,
    `Terms tested` = vapply(cols, function(k) length(unique(t$ID[t$collection == k])), 0L),
    `Inhibited (KO-like)` = vapply(cols, function(k) sum(sig$collection == k & sig$direction == "Inhibited"), 0L),
    Activated = vapply(cols, function(k) sum(sig$collection == k & sig$direction == "Activated"), 0L),
    `Smallest adj. P` = .rp_fmt_padj(vapply(cols, function(k) min(t$p.adjust[t$collection == k]), 0)),
    check.names = FALSE, stringsAsFactors = FALSE)
  answer <- sprintf(paste("%s term-direction pairs are over-represented: %s among the Inhibited (KO-like) genes and %s among the Activated",
    "genes, over %d collections. The plot shows the strongest terms of each collection, both directions side by side."),
    .fmt_n(nrow(sig)), .fmt_n(sum(sig$direction == "Inhibited")), .fmt_n(sum(sig$direction == "Activated")), length(cols))
  fig <- .rp_cached_fig(o, "enrichplot", pai_enrichplot(o$e, n = 4))
  sw <- .rp_L$W - fig$w - 2; iw <- .rp_inner(sw)
  side <- .rp_col(list(.rp_space(1), .rp_box(list(.rp_box_title("Reading the plot", iw), .rp_text_blk(paste(
    "Each row is a pathway. The dot on the right belongs to the Inhibited (KO-like) genes, the dot on the left to the Activated genes;",
    "the farther from the centre, the smaller the adjusted P. The number beside a dot counts the called genes in the pathway; open or",
    "shaded dots are not significant. N is the number of scored knockout genes the collection annotates, the background of the test."),
    size = 6.8, width = iw, after = 0.4)), width = sw)), sw)
  bw <- (.rp_L$W - 4) / 3; ib <- .rp_inner(bw)
  guide <- function(title, body) .rp_box(list(.rp_box_title(title, ib), .rp_text_blk(body, size = 6.6, width = ib, after = 0.4)), width = bw)
  c(.rp_page_head(o, "enrich_overview", answer),
    list(.rp_row(list(fig, side), gap = 2), .rp_space(1),
         .rp_section("Significant terms per collection", "adj. P < 0.05, BH within each collection and direction"),
         .rp_tbl(per, c(46, 30, 36, 30, 44)), .rp_space(2),
         .rp_section("How to read the next pages"),
         .rp_row(list(
           guide(sprintf("Bubble plot (page %d)", .rp_pn(o, "enrich_terms")), paste("Every tested term is a bubble: higher means more significant,",
             "left means its called genes are mostly Inhibited (KO-like), right mostly Activated. The ten strongest terms that do not repeat each",
             "other are numbered T1-T10 and listed beside the plot.")),
           guide(sprintf("Term-gene network (page %d)", .rp_pn(o, "enrich_terms")), paste("The five strongest terms in the middle, with the",
             "called genes they contain: Inhibited genes on the left, Activated on the right, placed by activity score. A gene linked to several",
             "terms is shared by them.")),
           guide(sprintf("Hub genes and modules (page %d)", .rp_pn(o, "enrich_hubs")), paste("Terms that share most of their genes form a",
             "module (one colour). A hub gene sits in many modules at once, so it links otherwise separate processes; a dagger marks genes",
             "that are in many terms of any kind."))), gap = 2)))
}
.rp_page_enrich_terms <- function(o, n = 12) {
  tt <- .as_tt(o$e); top <- .pick_terms(tt, n = 10)
  lead <- if (nrow(top)) sprintf("T1 is %s (%s, adj. P %s).", top$Description[1], if (top$direction[1] == "Inhibited") "Inhibited (KO-like) genes" else
    "Activated genes", .rp_fmt_padj(top$p.adjust[1])) else ""
  answer <- paste(sprintf("%s terms were tested; %s are significant in at least one direction.", .fmt_n(nrow(tt$table)),
    .fmt_n(sum(tt$table$p.adjust < 0.05))), lead)
  c(.rp_page_head(o, "enrich_terms", answer),
    list(.rp_cached_fig(o, "gobubble", pai_gobubble(o$e, width = .rp_L$W)), .rp_space(1.5),
         .rp_cached_fig(o, paste0("term_cnet_", n), pai_term_cnetplot(o$e, n = 5, n_genes = n))))
}
.rp_page_enrich_hubs <- function(o, n = 15) {
  h <- o$hubs
  say_dir <- function(k) if (inherits(h[[k]], "pai_hubs")) sprintf("%s: %s", .rp_hub_label(o, k),
    paste(utils::head(h[[k]]$genes$gene, 3), collapse = ", ")) else sprintf("%s: too few significant terms", .rp_hub_label(o, k))
  answer <- paste(if (identical(o$pathway_mode, "gsea"))
    "Pathways at each end of the activity ranking are grouped into modules of terms that share most of their genes." else
    "Significant terms of each direction are grouped into modules of terms that share most of their genes.",
    sprintf("Hub genes, found in the most modules: %s; %s.", say_dir("Inhibited"), say_dir("Activated")))
  rows <- list()
  one <- sum(vapply(h, inherits, TRUE, "pai_hubs")) == 1      # room for the module table of the only direction drawn
  for (k in c("Inhibited", "Activated")) {
    hk <- h[[k]]
    if (inherits(hk, "pai_hubs")) {
      rows <- c(rows, list(
        .rp_section(.rp_hub_label(o, k), sprintf("%d modules from %s terms", nrow(hk$modules), .fmt_n(nrow(hk$terms)))),
        .rp_row(list(.rp_cached_fig(o, sprintf("hub_%s_%d", k, n), pai_hubplot(hk, n_genes = n)),
                     .rp_cached_fig(o, sprintf("emap_%s", k), pai_emapplot(hk))), gap = 2)))
      if (one) {
        m <- utils::head(hk$modules, 8)
        rows <- c(rows, list(.rp_section("Modules", "all modules and their genes in the Excel sheets Hub_modules and Hub_genes"),
          .rp_tbl(data.frame(Module = m$module, `Most significant term` = m$label, Collection = m$label_source, Terms = m$n_terms,
                             Genes = m$n_genes, `adj. P` = .rp_fmt_padj(m$best_padj), check.names = FALSE),
                  c(16, 96, 22, 16, 16, 20), wrap = "Most significant term", max_lines = 2), .rp_space(2)))
      }
    } else {
      rows <- c(rows, list(.rp_section(.rp_hub_label(o, k)),
        .rp_text_blk(sprintf("Not drawn: %s", sub(" A larger cutoff.*$", "", as.character(hk))), size = 7, after = 3)))
    }
  }
  c(.rp_page_head(o, "enrich_hubs", answer), rows)
}
# A reading box beside (width w) or below a plot.
.rp_read <- function(text, w = .rp_L$W, title = "Reading the plot")
  .rp_box(list(.rp_box_title(title, .rp_inner(w)), .rp_text_blk(text, size = 6.8, width = .rp_inner(w), after = 0.4)), width = w)
# A plot of width <= 150 with its reading box in the column beside it.
.rp_fig_side <- function(fig, text) {
  sw <- .rp_L$W - fig$w - 2
  .rp_row(list(fig, .rp_col(list(.rp_space(1), .rp_read(text, sw)), sw)), gap = 2)
}
.rp_tnote <- function(o) sprintf("T1, T2, ... are the same terms on every page (listed on page %d).", .rp_pn(o, "enrich_terms"))
.rp_page_enrich_circle <- function(o) {
  c(.rp_page_head(o, "enrich_circle", paste("The ten top terms with all their scored member genes; the genes with the strongest",
      "activity in each term are named.", .rp_tnote(o))),
    list(.rp_cached_fig(o, "gocircle", pai_gocircle(o$e, n = 10, members = "all", label_genes = 1)), .rp_space(1.5),
         .rp_read(paste("Rings from the inside: the bar height is the term's \u2212log10 adjusted P (colour: balanced z, blue = mostly",
           "Inhibited genes, red = mostly Activated; dashed circle: adj. P = 0.05); the next ring splits the term's called genes into",
           "Inhibited | Activated; the outer ring places every scored member gene by activity score (solid line 0, dashed \u00b11), grey",
           "when the gene is not called. The named genes are the most extreme members of each term."))))
}
.rp_page_enrich_chord <- function(o) {
  c(.rp_page_head(o, "enrich_chord", paste("The chord diagram joins the five strongest terms to the genes they share; the Venn diagram",
      "counts the genes shared by the three strongest.", .rp_tnote(o))),
    list(.rp_cached_fig(o, "gochord", pai_gochord(o$e)), .rp_space(1.5),
         .rp_fig_side(.rp_cached_fig(o, "govenn", pai_govenn(o$e, width = 130)), paste("Chord: a ribbon joins a gene (left, coloured by",
           "activity score) to each term that contains it, so genes with many ribbons are shared. Venn: each region gives the number of",
           "called genes it holds, split into Inhibited (\u2193) and Activated (\u2191). Large shared regions mean the terms describe",
           "the same genes and should be read as one finding."))))
}
.rp_page_enrich_bars <- function(o, n = 15) {
  c(.rp_page_head(o, "enrich_bars", paste("Bars rank the top terms by significance and colour them by direction; the circle groups the",
      "genes of the five strongest terms by the terms they share.", .rp_tnote(o))),
    list(.rp_cached_fig(o, sprintf("gobar_%d", n), pai_gobar(o$e, n = n)), .rp_space(1.5),
         .rp_fig_side(.rp_cached_fig(o, "gocluster", pai_gocluster(o$e, width = 150)), paste("Bars: length = \u2212log10 adjusted P,",
           "colour = z-score (blue: mostly Inhibited genes, red: mostly Activated), numbers = Inhibited \u2193 and Activated \u2191 genes.",
           "Circle: the tree clusters genes by the terms that contain them; rings from the inside give the gene's activity score, then",
           "membership in T1-T5. Neighbouring genes share terms."))))
}
.rp_page_enrich_heat <- function(o, n = 30) {
  c(.rp_page_head(o, "enrich_heat", paste("Each column is one of the ten top terms and each row a called gene found in many of them;",
      "the colour is the gene's activity score.", .rp_tnote(o))),
    list(.rp_fig_side(.rp_cached_fig(o, sprintf("goheat_%d", n), pai_goheat(o$e, n_genes = n)), paste("Rows: the called genes found",
      "in the most top terms. Columns: the terms. A coloured cell means the gene belongs to the term (blue Inhibited, red Activated,",
      "darker = stronger); grey means it does not. The bars above each column give the term's Inhibited/Activated share and its",
      "balanced z. Rows coloured across many columns are genes shared by several terms."))))
}
.rp_page_enrich_armmap <- function(o) {
  c(.rp_page_head(o, "enrich_armmap", paste("Every scored gene of each top term is drawn on the two halves of its knockout signature.",
      .rp_tnote(o))),
    list(.rp_cached_fig(o, "armmap", pai_armmap(o$e)), .rp_space(1.5),
         .rp_read(paste("x: the up half of the signature, y: the down half, both oriented so that the warm top-right quadrant means",
           "Activated and the cool bottom-left means Inhibited (KO-like); the other two quadrants are Discordant. The framed first panel",
           "shows all scored genes as the reference; the percentages give the share of the panel's genes in each quadrant. A term whose",
           "genes crowd one coloured quadrant more than the reference points in that direction."))))
}
.rp_page_enrich_hubnet <- function(o) {
  h <- o$hubs; rows <- list()
  for (k in c("Inhibited", "Activated")) {
    hk <- h[[k]]
    rows <- c(rows, if (inherits(hk, "pai_hubs")) list(
      .rp_section(.rp_hub_label(o, k), sprintf("%d modules from %s terms", nrow(hk$modules), .fmt_n(nrow(hk$terms)))),
      .rp_row(list(.rp_cached_fig(o, sprintf("chord_%s", k), pai_chordplot(hk)),
                   .rp_cached_fig(o, sprintf("mcnet_%s", k), pai_enrich_cnetplot(hk))), gap = 2)) else
      list(.rp_section(.rp_hub_label(o, k)), .rp_text_blk(sprintf("Not drawn: %s", sub(" A larger cutoff.*$", "", as.character(hk))),
                                                          size = 7, after = 3)))
  }
  c(.rp_page_head(o, "enrich_hubnet", sprintf(paste("Left: hub genes (number = modules) joined to the modules they belong to.",
      "Right: the five largest modules (squares) with the genes in most modules, coloured by activity score; bold names are hubs.",
      "Modules as on page %d."), .rp_pn(o, "enrich_hubs"))), rows)
}
# A plot built once per report (pages are rebuilt while fitting their tables).
.rp_cached_fig <- function(o, key, expr) {
  if (is.null(o$figs[[key]])) o$figs[[key]] <- .rp_fig(suppressMessages(expr))
  o$figs[[key]]
}
.rp_compound_rows <- function(d) {
  ev <- ifelse(is.na(d$evidence_tier), "", paste0(d$evidence_tier, ifelse(d$evidence_confidence == "", "",
    paste0(", ", ifelse(d$evidence_confidence %in% names(.rp_conf_text), .rp_conf_text[d$evidence_confidence], d$evidence_confidence)))))
  pos <- ifelse(is.na(d$common_reference_percentile), "", sprintf("%.0f%%", 100 * d$common_reference_percentile))
  data.frame(Compound = d$compound, Score = .rp_sgn(d$activity_score), `adj. P` = .rp_fmt_padj(d$padj_coherent),
             `Mechanism (target)` = d$mechanism, Evidence = ev, `Usual position` = pos,
             Note = ifelse(d$common_hit, "common hit", ""), check.names = FALSE, stringsAsFactors = FALSE)
}
.rp_page_chem <- function(o, n = 8) {
  cc <- o$tabs$chem; tx <- o$text
  conc <- cc[cc$call == "Concordant", ]; opp <- cc[cc$call == "Opposing", ]; opp <- opp[order(opp$activity_score), ]
  tw <- c(29, 12, 16, 62, 25, 23, 19)
  st <- list(Note = function(v) if (v == "common hit") list(col = .rp_status_col[["Check"]], face = "bold") else list(),
             `Mechanism (target)` = function(v) if (v %in% c("unknown", "not a small molecule") || grepl("^legacy", v)) list(col = .tok$muted) else list())
  c(.rp_page_head(o, "chem", paste(tx$chem, tx$chem_top)),
    list(.rp_row(list(.rp_cached_fig(o, "arm_chem", pai_armplot(o$xc, n_label = 6)), .rp_cached_fig(o, "dot_chem", pai_dotplot(o$xc, n = 10))), gap = 2),
         .rp_section("Strongest Concordant (similar profile)", sprintf("%s in total", .fmt_n(nrow(conc)))),
         .rp_tbl(.rp_compound_rows(utils::head(conc, n)), tw, styles = st),
         .rp_section("Strongest Opposing (reverse profile)", sprintf("%s in total", .fmt_n(nrow(opp)))),
         .rp_tbl(.rp_compound_rows(utils::head(opp, n)), tw, styles = st),
         .rp_text_blk(paste("Mechanism: Chem taxonomy 3.0.1. Evidence: A curated mechanism (ChEMBL, Guide to PHARMACOLOGY, LINCS 2020 or an author",
           "decision), B bioactivity only, C/D legacy label only, U unknown, X not a small molecule; then how many sources agree. Usual position:",
           "the compound's median percentile in 12 unrelated public perturbation experiments, in 9 groups (100% = always the most similar). A common hit sits at the",
           "same end in most of them, so its match here is less specific."), size = 6.3, col = .tok$muted)))
}
.rp_page_classes <- function(o) {
  kt <- o$classes
  sigc <- kt[!is.na(kt$padj_concordant) & kt$padj_concordant < 0.05, ]; sigo <- kt[!is.na(kt$padj_opposing) & kt$padj_opposing < 0.05, ]
  say <- function(s) if (nrow(s)) sprintf("%d (%s)", nrow(s), paste(utils::head(s$mechanism_class, 2), collapse = "; ")) else "none"
  answer <- sprintf(paste("If a match reflects a mechanism, compounds that share it should agree. Mechanism classes over-represented",
    "among the calls: Concordant, %s; Opposing, %s."), say(sigc), say(sigo))
  fp <- tryCatch(suppressMessages(pai_familyplot(o$xc)), error = function(e) NULL)
  fw <- if (is.null(fp)) 92 else attr(fp, "size_mm")[["width"]]
  sw <- .rp_L$W - fw - 2; iw <- .rp_inner(sw)
  side <- .rp_col(list(.rp_space(1), .rp_box(list(.rp_box_title("Curated drug families", iw),
    .rp_text_blk(paste("Each bar is a hand-checked family of compounds with the same direct target (taxonomy 3.0.1, evidence tiers A and B).",
      "It shows the share of the family's compounds called on each side, with the count. One member agreeing is weak support; most",
      "members agreeing is strong support."), size = 6.8, width = iw, after = 0.4)), width = sw)), sw)
  cp <- .rp_plot_classes(o$classes)
  tab <- utils::head(kt[kt$mechanism_class != .rp_no_class, ], 8)
  td <- data.frame(`Mechanism class` = tab$mechanism_class, Compounds = tab$n, Concordant = tab$concordant,
                   `adj. P (Conc.)` = .rp_fmt_padj(tab$padj_concordant), Opposing = tab$opposing, `adj. P (Opp.)` = .rp_fmt_padj(tab$padj_opposing),
                   check.names = FALSE)
  c(.rp_page_head(o, "classes", answer),
    list(.rp_row(list(if (is.null(fp)) .rp_text_blk("No curated drug family has members among the scored compounds.", size = 7, width = fw)
                      else .rp_fig(fp), side), gap = 2),
         if (!is.null(cp)) .rp_fig(cp) else .rp_text_blk("No mechanism class has called compounds."),
         .rp_space(1),
         .rp_section("Classes with the most calls", "all classes in the Excel sheet Chem_classes"),
         .rp_tbl(td, c(80, 20, 20, 24, 18, 24), align = c("left", "right", "right", "right", "right", "right"))))
}
.rp_page_cautions <- function(o) {
  ch <- o$checks; ct <- o$tabs$crispr; cc <- o$tabs$chem
  cards <- list(); iw <- .rp_inner(.rp_L$W, bar = TRUE)
  add <- function(status, title, body) cards[[length(cards) + 1]] <<- .rp_box(list(
    .rp_blk(grid::grobTree(
      grid::textGrob(toupper(status), x = 0, y = .rp_top(0.7), just = c("left", "top"), gp = .gp(6.2, .rp_status_col[[status]], "bold")),
      grid::textGrob(title, x = grid::unit(14, "mm"), y = .rp_top(0.4), just = c("left", "top"), gp = .gp(7.8, .tok$ink, "bold"))), 5, iw),
    .rp_text_blk(body, size = 6.9, width = iw, after = 0.3)), fill = .rp_status_tint[[status]], bar = .rp_status_col[[status]])
  s <- stats::setNames(ch$status, ch$check)
  if (s[["Proliferation signal"]] != "OK") {
    sp <- o$prolif[o$prolif$strong, ]
    add(s[["Proliferation signal"]], "Strong proliferation signal in your ranking", sprintf(paste(
      "%s. Many L1000 signatures slow or speed cell growth, so a contrast with a cell-cycle component matches them for that reason alone.",
      "Tissue contrasts that differ in dividing cells (for example immune or glial cells) often show this. Before reading a top hit as",
      "specific, check whether it is a cell-cycle gene or an antiproliferative compound, and whether its call holds in a ranking without",
      "the cell-cycle genes."), paste(sprintf("%s: NES %s, adj. P %s", .rp_hall_name(sp$gene_set), .fmt_num(sp$nes, 2),
      .rp_fmt_padj(sp$padj)), collapse = "; ")))
  }
  if (s[["Discordant share"]] != "OK")
    add(s[["Discordant share"]], "Many signatures moved in both directions at once", sprintf(paste(
      "%.0f%% of the gene knockouts and %.0f%% of the compounds are Discordant: the up-genes and the down-genes of their signature both",
      "moved the same way in your contrast. This happens when large groups of genes change together, as with a shift in cell-type",
      "composition, and it lowers the number of clear calls. A deconvolution or a cell-type-adjusted contrast can help."),
      100 * mean(ct$call == "Discordant"), 100 * mean(cc$call == "Discordant")))
  tc <- .rp_called(cc)[.rp_called(cc)$common_hit, ]; tc <- tc[order(-abs(tc$activity_score)), ]
  if (nrow(tc))
    add(if (s[["Common hits"]] == "OK") "Check" else s[["Common hits"]], "Some compound matches are common hits", sprintf(paste(
      "%s of the %s called compounds usually sit at the same end in unrelated perturbation experiments (12 public GEO contrasts in 9",
      "groups), so their match is partly generic: %s. Their rank says less about your contrast than that of a compound whose usual",
      "position is near the middle. They are marked 'common hit' on page %d and in the Excel sheet Chem_compounds."),
      .fmt_n(nrow(tc)), .fmt_n(nrow(.rp_called(cc))), paste(utils::head(tc$compound, 6), collapse = ", "), .rp_pn(o, "chem")))
  um <- .rp_called(ct)[!.rp_called(ct)$target_measured, ]; um <- um[order(-abs(um$activity_score)), ]
  if (nrow(um))
    add(if (s[["Unmeasured targets"]] == "OK") "Check" else s[["Unmeasured targets"]], "Some called genes are not in your table", sprintf(paste(
      "%s of the %s gene calls are for genes your table does not measure (for example %s). Their call rests on the genes that change",
      "after the knockout, not on the gene itself, so it cannot be checked against the gene's own expression."),
      .fmt_n(nrow(um)), .fmt_n(nrow(.rp_called(ct))), paste(utils::head(um$target, 5), collapse = ", ")))
  if (s[["Pathways"]] != "OK")
    add("Check", "No pathway among the called genes", sprintf(paste("The called genes do not share a pathway beyond chance. Page %d shows the ranked",
      "test over all genes instead; it describes where pathways sit along the activity ranking and is weaker evidence than",
      "over-representation among the calls."), .rp_pn(o, "pathways")))
  if (s[["Gene coverage"]] != "OK")
    add(s[["Gene coverage"]], "Low gene coverage", sprintf("Your table measures %.0f%% of the signature genes. Use a table with every tested gene, not only the significant ones.", 100 * o$coverage[["crispr"]]))
  if (s[["Formal calls"]] != "OK")
    add(s[["Formal calls"]], "Few or no calls", "Few targets pass the call thresholds. The rankings are still informative but should be read as descriptive.")
  if (!length(cards)) cards <- list(.rp_text_blk("No check needs attention."))
  general <- list(
    c("Similarity, not proof.", "A call means that the genes changed in your contrast overlap, in a direction-consistent way, with those changed after a knockout or a compound in L1000 experiments. It does not show that the gene or drug target is altered in your samples."),
    c("Cell-line context.", "L1000 signatures come from short-term perturbations of cultured, mostly cancer, cell lines, measured on about 1,000 landmark genes with the rest inferred. Tissue contrasts also reflect cell-type composition."),
    c("Compounds are hypotheses.", "An Opposing compound reverses the profile in cell lines; this is a lead to test, not evidence of benefit or safety."),
    c("Mechanism labels.", "Taxonomy 3.0.1 gives each compound a target and an evidence tier. Rely on tier A; treat tier B with care and lower tiers as unknown."),
    c("Multiple testing.", "Adjusted P values control the false discovery rate within each library; a long list of calls is expected to include some false ones."))
  nxt <- list(
    c("Open the Excel workbook.", "Sheet Summary repeats this report; CRISPR_targets and Chem_compounds hold every target with its halves, P values, mechanism and flags."),
    c("Look for agreement across levels.", "A gene call, its pathway and a compound acting on that gene pointing the same way is stronger than any single call."),
    c("Confirm in your own data or experiments.", "For example the expression or activity of the top genes, or a test of the top compounds in a relevant model."))
  sp_cards <- unlist(lapply(cards, function(cd) list(cd, .rp_space(0.6))), recursive = FALSE)
  c(.rp_page_head(o, "cautions",
      sprintf("%d of the %d quality checks need attention for this contrast. Each flag is explained below, followed by notes that apply to every L1000PAI result.",
              sum(ch$status != "OK"), nrow(ch))),
    list(.rp_section("Flags for this contrast")), sp_cards,
    list(.rp_space(1), .rp_section("Notes for every report"), .rp_bullets(general, size = 6.9), .rp_space(0.5),
         .rp_section("Suggested next steps"), .rp_bullets(nxt, size = 6.9)))
}
.rp_page_methods <- function(o) {
  mx <- pai_metadata(o$x); cp <- mx$classify_params
  sc <- .pai_semantics("crispr", .pai_score_direction(o$x, "crispr")); sx <- .pai_semantics("chem", .pai_score_direction(o$xc, "chem"))
  meth <- c(sprintf(paste("Genes were ranked by %s; for duplicate genes the entry with the largest absolute statistic was kept. Each L1000",
    "consensus signature has an up half and a down half (the genes that go up or down after the perturbation). Both halves were tested",
    "against the ranking with fgsea (multilevel; gene sets of 10\u2013500 genes; eps 1e\u221210; seed %s; observed-tail P values)."), o$prep$stat_method, o$meta$seed),
    sprintf(paste("Activity score: knockouts %s; compounds %s. A call requires |activity score| \u2265 %s, a coherent P value below %s after",
      "Benjamini\u2013Hochberg adjustment over all pairs, and both halves aligned. Pairs whose two halves are both significant with the same",
      "sign are Discordant."), sub("Activity score  ", "", sc$score_label), sub("Activity score  ", "", sx$score_label),
      cp$activity_threshold %||% 1, cp$classify_p_value_threshold %||% 0.05),
    sprintf(paste("Pathways: over-representation of the Activated and of the Inhibited knockout targets, tested separately against all",
      "scored targets, in %s; when no term is significant, pathway GSEA over all targets ranked by activity score. Hub genes: the",
      "significant terms of one direction are grouped into modules by gene overlap (Jaccard similarity 0.25, average linkage) and genes",
      "are ranked by the number of modules they belong to (pai_enrich_hubs). Compound mechanisms:",
      "Chem taxonomy %s (ChEMBL 37, Guide to PHARMACOLOGY 2026.3.1, LINCS 2020; CC BY-SA 4.0). Proliferation: fgsea of the input",
      "ranking on the Hallmark E2F targets and G2M checkpoint sets (%s). Common hits: reference percentiles over 12 public GEO",
      "perturbation contrasts in 9 groups (reference v2, listed in the package file extdata/report/README.md)."), paste(o$meta$collections, collapse = ", "), o$meta$taxonomy, o$meta$msigdb))
  par <- data.frame(Setting = c("Library mode", "Score direction (CRISPR / Chem)", "Gene-set size", "fgsea eps", "Seed", "P value semantics",
                                "Pair P adjustment", "Activity threshold", "Call P threshold", "Coherent test required"),
                    Value = c(mx$library_mode_used %||% "paired",
                              paste(.pai_score_direction(o$x, "crispr"), .pai_score_direction(o$xc, "chem"), sep = " / "),
                              "10\u2013500", "1e\u221210", as.character(o$meta$seed), "observed tail", "BH", as.character(cp$activity_threshold %||% 1),
                              as.character(cp$classify_p_value_threshold %||% 0.05), "yes"), stringsAsFactors = FALSE)
  v <- o$meta$versions
  sw <- data.frame(Component = c("R", "L1000PAI", "L1000PAIdata", "fgsea", "org.Hs.eg.db", "Hallmark sets", "CRISPR library (SHA-256)",
                                 "Chem library (SHA-256)", "Chem taxonomy", "Common-hit reference"),
                   Version = c(paste(R.version$major, R.version$minor, sep = "."), v[["L1000PAI"]], v[["L1000PAIdata"]], v[["fgsea"]],
                               v[["org.Hs.eg.db"]], o$meta$msigdb, paste0(substr(o$meta$sha[["crispr"]], 1, 16), "\u2026"),
                               paste0(substr(o$meta$sha[["chem"]], 1, 16), "\u2026"), o$meta$taxonomy, o$meta$reference_version), stringsAsFactors = FALSE)
  refs <- list(
    "Subramanian A, et al. A next generation connectivity map: L1000 platform and the first 1,000,000 profiles. Cell 2017;171:1437\u20131452.",
    "Evangelista JE, et al. SigCom LINCS: data and metadata search engine for a million gene expression signatures. Nucleic Acids Res 2022;50:W697\u2013W709.",
    "Korotkevich G, et al. Fast gene set enrichment analysis. bioRxiv 2021. doi:10.1101/060012.",
    "Liberzon A, et al. The Molecular Signatures Database Hallmark gene set collection. Cell Syst 2015;1:417\u2013425.",
    "Zdrazil B, et al. The ChEMBL Database in 2023. Nucleic Acids Res 2024;52:D1180\u2013D1192.",
    "Harding SD, et al. The IUPHAR/BPS Guide to PHARMACOLOGY in 2024. Nucleic Acids Res 2024;52:D1438\u2013D1449.",
    sprintf("Chiu Y. L1000PAI: direction-aware perturbation activity inference from ranked transcriptomic signatures. R package version %s.", v[["L1000PAI"]]))
  what <- c(pdf = "this report", xlsx = "every table behind the report, with a column dictionary (sheet README)",
            rds = "the R objects (results, tables, settings) for further analysis; readRDS()")
  files <- data.frame(File = basename(o$meta$files), Contents = unname(what[names(o$meta$files)]), stringsAsFactors = FALSE)
  c(.rp_page_head(o, "methods", "One call of pai_deg_report() produced this report and its Excel workbook from your table. Everything needed to repeat it is listed here."),
    list(.rp_section("Methods"), .rp_text_blk(meth[1], size = 6.9, after = 1.4), .rp_text_blk(meth[2], size = 6.9, after = 1.4),
         .rp_text_blk(meth[3], size = 6.9, after = 1.8),
         .rp_row(list(.rp_col(list(.rp_section("Settings", width = 91), .rp_tbl(par, c(50, 41), size = 6.3, align = c("left", "left"))), 91),
                   .rp_col(list(.rp_section("Software and data", width = 91), .rp_tbl(sw, c(40, 51), size = 6.3, align = c("left", "left"))), 91)), gap = 4),
         .rp_section("Files"), .rp_tbl(files, c(100, 86), size = 6.3, align = c("left", "left"), wrap = c("File", "Contents")),
         .rp_space(1), .rp_section("References"), .rp_bullets(lapply(refs, function(r) r), size = 6.4, bullet = "\u2013", gap = 0.6)))
}

# =========================================================================================================
# 8. Excel workbook
# =========================================================================================================
.deg_dictionary <- function(o) {
  pn <- function(k) .rp_pn(o, k)
  d <- function(sheet, ...) { a <- c(...); data.frame(sheet = sheet, column = names(a), meaning = unname(a), stringsAsFactors = FALSE) }
  rbind(
    d("Summary", item = "What is reported", value = "The number or sentence shown in the PDF"),
    d("Checks", check = "Quality check (PDF page 1)", status = sprintf("OK, Check (read PDF page %d) or Caution (interpret with care)", pn("cautions")),
      result = "What was found", meaning = "Why the check matters", value = "Number behind the status"),
    d("Input_QC", input_row = "Row number in your table", input_id = "Identifier as given", gene = "Gene symbol used",
      statistic = "Ranking statistic", log_fold_change = "Fold change column of your table, if present", adj_p = "Adjusted P value of your table, if present",
      expression = "Average expression column of your table, if present", status = "used; no gene symbol; no statistic; or duplicate gene, smaller |statistic| dropped"),
    d("Proliferation", gene_set = "MSigDB Hallmark gene set", size = "Genes of the set found in your table",
      nes = "Normalised enrichment score in your ranking (positive = among your increased genes)",
      pvalue = "fgsea P value", padj = "BH-adjusted P value", strong = "TRUE when padj < 0.05 and |NES| >= 2 (triggers the proliferation check)",
      leading_edge = "Up to 30 leading-edge genes"),
    d("CRISPR_targets", target = "Knocked-out gene", call = "Activated, Inhibited (KO-like), Discordant or Not called",
      activity_score = "NES Down - NES Up; positive = your contrast is the reverse of the knockout", activity_percentile = "Rank of the score among all targets (100% = most Activated)",
      padj_coherent = "BH-adjusted coherent P value; calls need < 0.05", nes_up_arm = "NES of the genes that go up after the knockout",
      nes_down_arm = "NES of the genes that go down after the knockout", padj_up_arm = "Adjusted P of the up half", padj_down_arm = "Adjusted P of the down half",
      p_activation = "One-sided P for activation", p_inhibition = "One-sided P for inhibition", target_measured = "Is the knocked-out gene in your table",
      target_statistic = "That gene's own statistic in your table", target_percentile = "Its rank in your table (100% = most increased)"),
    d("Pathways", collection = "Gene-set collection", id = "Term identifier", term = "Pathway",
      end = "GSEA: end of the activity ranking where the pathway sits (Activated or Inhibited); ORA: direction of the calls (column direction)",
      nes = "GSEA only: normalised enrichment along the activity ranking", pvalue = "P value", padj = "BH-adjusted P value",
      size = "Genes of the pathway among the scored targets", genes = "Genes driving the result"),
    d("Chem_compounds", compound = "Compound", call = "Concordant (similar profile), Opposing (reverse), Discordant or Not called",
      activity_score = "NES Up - NES Down; positive = resembles the compound", activity_percentile = "Rank among all compounds (100% = most Concordant)",
      padj_coherent = "BH-adjusted coherent P value", nes_up_arm = "NES of the genes that go up after the compound", nes_down_arm = "NES of the genes that go down",
      mechanism = "Target class from Chem taxonomy 3.0.1", target_genes = "Primary target genes", mechanism_family = "Family class",
      mechanism_class = "Broad class", evidence_tier = "A curated mechanism; B bioactivity only; C/D legacy label only; U unknown; X not a small molecule",
      evidence_confidence = "corroborated (2+ sources), single, conflict, author_curated", evidence_sources = "Sources behind the call",
      evidence_records = "Record identifiers (ChEMBL, GtoPdb, LINCS 2020 IDs)", common_reference_percentile = "Median percentile across 12 unrelated public perturbation contrasts (9 groups)",
      percentile_minus_reference = "Percentile here minus the usual percentile (large = more specific to your contrast)",
      common_response = "usually_top (reference >= 80%) or usually_bottom (<= 20%)", common_hit = "Called at the end where the compound usually sits"),
    d("Chem_classes", mechanism_class = "Broad mechanism class (tiers A and B)", n = "Compounds in the class", concordant = "Called Concordant",
      opposing = "Called Opposing", discordant = "Discordant", expected_concordant = "Expected Concordant if calls ignored the class",
      expected_opposing = "Expected Opposing", p_concordant = "One-sided Fisher P, more Concordant than expected", p_opposing = "Same for Opposing",
      padj_concordant = "BH over classes with 5+ compounds", padj_opposing = "BH over classes with 5+ compounds",
      top_concordant = "Up to 5 Concordant members", top_opposing = "Up to 5 Opposing members"),
    d("Chem_families", family = sprintf("Curated drug families (direct mechanism; tiers A and B), as drawn on PDF page %d", pn("classes"))),
    d("Pathway_terms", collection = "Gene-set collection", ID = "Term identifier", Description = "Term",
      padj_activated = "Adjusted P among the Activated genes", padj_inhibited = "Adjusted P among the Inhibited (KO-like) genes",
      p.adjust = "Smaller of the two", direction = "Direction with the smaller adjusted P", count = "Called genes in the term",
      n_activated = "Activated genes in the term", n_inhibited = "Inhibited genes in the term",
      zscore = "(Activated - Inhibited) / sqrt(count)", zscore_balanced = "Activated count against the share of Activated among all calls",
      mean_activity = "Mean activity score of the called genes in the term", geneID = "Called genes in the term"),
    d("Hub_genes", direction = "Inhibited (KO-like) or Activated genes, or the end of the activity ranking", gene = "Gene",
      n_terms = "Significant terms containing the gene", n_modules = "Modules containing the gene", modules = "Those modules",
      activity_score = "The gene's activity score", breadth = "Tested terms containing the gene",
      broad = "In more tested terms than 95% of genes (a generally well-annotated gene)", hub_rank = "Rank by modules, then terms"),
    d("Hub_modules", direction = "As in Hub_genes", module = "Module (M1 = most significant)", label = "Most significant term of the module",
      label_source = "Collection of that term", best_padj = "Its adjusted P", n_terms = "Terms in the module", n_genes = "Genes in the module",
      collections = "Collections contributing terms", cluster = "Internal cluster number"),
    d("Flags", flag = "Kind of flag", item = "Gene, compound or gene set", detail = "Why it is flagged"),
    d("Parameters", setting = "Setting, software or data version", value = "Value"))
}
.deg_flags <- function(o) {
  ct <- o$tabs$crispr; cc <- o$tabs$chem
  f <- list()
  sp <- o$prolif[o$prolif$strong, ]
  if (nrow(sp)) f$p <- data.frame(flag = "proliferation", item = sp$gene_set, detail = sprintf("NES %.2f, adj. P %.2g", sp$nes, sp$padj))
  ch <- cc[cc$common_hit, ]
  if (nrow(ch)) f$c <- data.frame(flag = "common hit", item = ch$compound,
    detail = sprintf("%s here; usual percentile %.0f%% (%s)", ch$call, 100 * ch$common_reference_percentile, ch$common_response))
  um <- .rp_called(ct)[!.rp_called(ct)$target_measured, ]
  if (nrow(um)) f$u <- data.frame(flag = "target not measured", item = um$target, detail = sprintf("%s call; gene absent from your table", um$call))
  low <- .rp_called(cc)[.rp_called(cc)$evidence_tier %in% c("C", "D", "U", "X"), ]
  if (nrow(low)) f$l <- data.frame(flag = "no reliable mechanism", item = low$compound, detail = sprintf("%s call; evidence tier %s", low$call, low$evidence_tier))
  if (!length(f)) return(data.frame(flag = character(), item = character(), detail = character()))
  do.call(rbind, unname(f))
}
# openxlsx 4.2.9 lists drawing parts in the sheet relationships without writing them; strict readers
# (openpyxl, some Excel builds) then refuse the file. Remove relationships whose internal target is absent.
.rp_xlsx_drop_dangling_rels <- function(file) {
  file <- normalizePath(file, mustWork = TRUE)        # zip() works from inside the unpacked folder
  td <- tempfile("xlsx_"); dir.create(td); on.exit(unlink(td, recursive = TRUE), add = TRUE)
  utils::unzip(file, exdir = td)
  changed <- FALSE
  for (r in list.files(file.path(td, "xl", "worksheets", "_rels"), full.names = TRUE)) {
    x <- paste(readLines(r, warn = FALSE, encoding = "UTF-8"), collapse = "")
    items <- regmatches(x, gregexpr("<Relationship [^>]*/>", x))[[1]]
    gone <- vapply(items, function(it) {
      if (grepl('TargetMode="External"', it, fixed = TRUE)) return(FALSE)
      !file.exists(file.path(td, "xl", "worksheets", sub('.*Target="([^"]+)".*', "\\1", it)))
    }, TRUE)
    if (!any(gone)) next
    for (it in items[gone]) x <- sub(it, "", x, fixed = TRUE)
    writeLines(x, r, useBytes = TRUE); changed <- TRUE
  }
  if (!changed) return(invisible(FALSE))
  files <- list.files(td, recursive = TRUE, all.files = TRUE)
  files <- c("[Content_Types].xml", setdiff(files, "[Content_Types].xml"))
  unlink(file); zip::zip(file, files = files, root = td, mode = "mirror")
  invisible(TRUE)
}
.deg_write_xlsx <- function(o, file) {
  wb <- openxlsx::createWorkbook(creator = "L1000PAI")
  openxlsx::modifyBaseFont(wb, fontSize = 10, fontName = "Arial")
  hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#ECEBE7", border = "bottom", borderColour = "#BDBBB5", wrapText = TRUE, valign = "top")
  add <- function(name, df, widths = NULL, num = list(), calls = NULL) {
    openxlsx::addWorksheet(wb, name, gridLines = FALSE)
    openxlsx::writeData(wb, name, df, headerStyle = hs, withFilter = nrow(df) > 0)
    openxlsx::freezePane(wb, name, firstRow = TRUE)
    openxlsx::setColWidths(wb, name, seq_along(df), widths %||% pmin(pmax(nchar(names(df)) + 3, 10), 45))
    for (cl in names(num)) if (cl %in% names(df) && nrow(df)) openxlsx::addStyle(wb, name, openxlsx::createStyle(numFmt = num[[cl]]),
      rows = seq_len(nrow(df)) + 1, cols = match(cl, names(df)), gridExpand = TRUE, stack = TRUE)
    if (!is.null(calls) && calls %in% names(df) && nrow(df)) {
      j <- match(calls, names(df))
      for (v in intersect(names(.rp_call_key), unique(df[[calls]]))) {
        r <- which(df[[calls]] == v)
        openxlsx::addStyle(wb, name, openxlsx::createStyle(fgFill = .pal$tint[[.rp_call_key[[v]]]], fontColour = .pal$text[[.rp_call_key[[v]]]], textDecoration = "bold"),
          rows = r + 1, cols = j, stack = TRUE)
      }
    }
  }
  p3 <- "0.000"; pe <- "0.00E+00"; pc <- "0.0%"
  # README
  openxlsx::addWorksheet(wb, "README", gridLines = FALSE)
  intro <- c("L1000PAI report: companion workbook", o$meta$contrast, o$meta$byline, "",
    "Each sheet holds one table behind the PDF report. Rows can be sorted and filtered; the column dictionary below explains every column.",
    "Calls: Activated / Inhibited (KO-like) for gene knockouts; Concordant (similar profile) / Opposing (reverse profile) for compounds.",
    sprintf("Scores describe similarity to L1000 cell-line signatures, not proof of a mechanism. Read PDF page %d before drawing conclusions.",
            .rp_pn(o, "cautions")), "")
  openxlsx::writeData(wb, "README", data.frame(x = intro), colNames = FALSE)
  openxlsx::addStyle(wb, "README", openxlsx::createStyle(fontSize = 14, textDecoration = "bold"), rows = 1, cols = 1)
  pn <- function(k) .rp_pn(o, k)
  terms_tab <- if (o$pathway_mode == "ora") tryCatch(pai_term_table(o$e)$table, error = function(e) NULL) else NULL
  if (!is.null(terms_tab)) terms_tab$key <- NULL
  hubs <- Filter(function(h) inherits(h, "pai_hubs"), o$hubs %||% list())
  hub_genes <- if (length(hubs)) do.call(rbind, lapply(names(hubs), function(k) data.frame(direction = .rp_hub_label(o, k), hubs[[k]]$genes,
    stringsAsFactors = FALSE, check.names = FALSE))) else NULL
  hub_modules <- if (length(hubs)) do.call(rbind, lapply(names(hubs), function(k) {
    m <- hubs[[k]]$modules; data.frame(direction = .rp_hub_label(o, k), m[, c("module", setdiff(names(m), "module"))], stringsAsFactors = FALSE) })) else NULL
  sheets <- data.frame(sheet = c("Summary", "Checks", "Input_QC", "Proliferation", "CRISPR_targets", "Pathways",
                                 if (!is.null(terms_tab)) "Pathway_terms", if (!is.null(hub_genes)) c("Hub_genes", "Hub_modules"),
                                 "Chem_compounds", "Chem_classes", "Chem_families", "Flags", "Parameters"),
    contents = c("Headline numbers and key findings (PDF page 1)", "Quality checks with status (page 1)",
                 sprintf("Every row of your table and how it was used (page %d)", pn("input")), sprintf("Proliferation check (page %d)", pn("input")),
                 sprintf("All gene knockouts with both halves, P values and own-gene values (pages %d-%d)", pn("crispr"), pn("evidence")),
                 switch(o$pathway_mode, ora = sprintf("Pathways over-represented among the calls, per direction (pages %d-%d)", pn("pathways"), pn("enrich_overview")),
                        gsea = sprintf("Pathway GSEA over all targets ranked by activity (page %d)", pn("pathways")),
                        sprintf("Pathway analysis did not run (page %d)", pn("pathways"))),
                 if (!is.null(terms_tab)) sprintf("One row per tested term with both directions and z-scores (page %d)", pn("enrich_terms")),
                 if (!is.null(hub_genes)) c(sprintf("Hub genes of the pathway modules (page %d)", pn("enrich_hubs")),
                                            sprintf("Pathway modules (page %d)", pn("enrich_hubs"))),
                 sprintf("All compounds with mechanism, evidence and common-hit reference (page %d)", pn("chem")),
                 sprintf("Mechanism classes among the calls (page %d)", pn("classes")), sprintf("Curated drug families among the calls (page %d)", pn("classes")),
                 sprintf("Every flagged item (page %d)", pn("cautions")), sprintf("Settings, versions and data fingerprints (page %d)", pn("methods"))))
  openxlsx::writeData(wb, "README", sheets, startRow = length(intro) + 1, headerStyle = hs)
  dic <- .deg_dictionary(o); dic <- dic[dic$sheet %in% sheets$sheet, , drop = FALSE]
  openxlsx::writeData(wb, "README", dic, startRow = length(intro) + nrow(sheets) + 3, headerStyle = hs)
  openxlsx::setColWidths(wb, "README", 1:3, c(22, 34, 110))
  # Data sheets
  ct <- o$tabs$crispr; cc <- o$tabs$chem
  summ <- data.frame(item = c("Contrast", "Genes ranked", "Ranking statistic", "Coverage, knockout signatures", "Coverage, compound signatures",
                              "Gene knockouts scored", "Activated", "Inhibited (KO-like)", "Discordant (knockouts)",
                              "Compounds scored", "Concordant", "Opposing", "Discordant (compounds)", "Pathway view",
                              paste("Finding:", c("input", "gene knockouts", "top genes", "pathways", "compounds", "top compounds", "attention"))),
    value = c(o$meta$contrast, nrow(o$prep$rank_input), o$prep$stat_method, .rp_pct(o$coverage[["crispr"]]), .rp_pct(o$coverage[["chem"]]), nrow(ct), sum(ct$call == "Activated"), sum(ct$call == "Inhibited (KO-like)"),
              sum(ct$call == "Discordant"), nrow(cc), sum(cc$call == "Concordant"), sum(cc$call == "Opposing"), sum(cc$call == "Discordant"),
              switch(o$pathway_mode, ora = "over-representation among calls", gsea = "GSEA over all targets (no over-represented pathway)",
                     "did not run"),
              unlist(o$text[c("input", "crispr", "crispr_top", "pathway", "chem", "chem_top", "attention")])), stringsAsFactors = FALSE)
  add("Summary", summ, c(34, 120))
  openxlsx::addStyle(wb, "Summary", openxlsx::createStyle(wrapText = TRUE, valign = "top"), rows = 2:(nrow(summ) + 1), cols = 1:2, gridExpand = TRUE, stack = TRUE)
  add("Checks", o$checks[, c("check", "status", "result", "meaning", "value")], c(22, 10, 50, 90, 10))
  for (v in names(.rp_status_col)) { r <- which(o$checks$status == v); if (length(r)) openxlsx::addStyle(wb, "Checks",
    openxlsx::createStyle(fontColour = .rp_status_col[[v]], textDecoration = "bold", fgFill = .rp_status_tint[[v]]), rows = r + 1, cols = 2, stack = TRUE) }
  add("Input_QC", o$prep$qc, NULL, list(statistic = p3, log_fold_change = p3, adj_p = pe, expression = p3))
  add("Proliferation", o$prolif, c(28, 8, 10, 12, 12, 10, 90), list(nes = p3, pvalue = pe, padj = pe))
  add("CRISPR_targets", ct, NULL, list(activity_score = p3, activity_percentile = pc, padj_coherent = pe, nes_up_arm = p3, nes_down_arm = p3,
      padj_up_arm = pe, padj_down_arm = pe, p_activation = pe, p_inhibition = pe, target_statistic = p3, target_percentile = pc), calls = "call")
  pw <- o$pw; pw$end_label <- NULL
  add("Pathways", pw, NULL, list(padj = pe, nes = p3, pvalue = pe))
  if (!is.null(terms_tab)) add("Pathway_terms", terms_tab, c(10, 13, 50, 12, 12, 12, 12, 8, 10, 10, 10, 12, 12, 80),
    list(padj_activated = pe, padj_inhibited = pe, p.adjust = pe, zscore = "0.00", zscore_balanced = "0.00", mean_activity = p3))
  if (!is.null(hub_genes)) {
    add("Hub_genes", hub_genes, c(22, 12, 9, 10, 40, 12, 9, 8, 9), list(activity_score = p3))
    add("Hub_modules", hub_modules, NULL, list(best_padj = pe))
  }
  add("Chem_compounds", cc, c(22, 13, 12, 12, 12, 10, 10, 40, 22, 34, 34, 8, 14, 18, 40, 12, 12, 14, 10),
      list(activity_score = p3, activity_percentile = pc, padj_coherent = pe, nes_up_arm = p3, nes_down_arm = p3,
           common_reference_percentile = pc, percentile_minus_reference = "+0.0%;-0.0%;0.0%"), calls = "call")
  add("Chem_classes", o$classes, c(48, 8, 10, 10, 10, 12, 12, 12, 12, 12, 12, 50, 50),
      list(expected_concordant = "0.0", expected_opposing = "0.0", p_concordant = pe, p_opposing = pe, padj_concordant = pe, padj_opposing = pe))
  add("Chem_families", if (ncol(o$families)) o$families else data.frame(note = "No curated drug family could be drawn."))
  add("Flags", o$flags, c(22, 24, 70))
  add("Parameters", o$params, c(36, 90))
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
  .rp_xlsx_drop_dangling_rels(file)
  invisible(file)
}


# =========================================================================================================
# 9. Main
# =========================================================================================================
.rp_report_dir <- function() system.file("extdata", "report", package = "L1000PAI")
.rp_reference <- function() {
  f <- file.path(.rp_report_dir(), "chem_common_response_reference_v2.csv")
  utils::read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
}
.rp_prolif_sets <- function() {
  f <- file.path(.rp_report_dir(), "hallmark_proliferation_sets.csv")
  d <- utils::read.csv(f, stringsAsFactors = FALSE)
  sets <- split(d$gene_symbol, d$gene_set)
  attr(sets, "db_version") <- d$msigdb_version[1]
  sets
}
# fgsea warns that P values below eps (1e-10) are reported as eps and that some P values may be
# overestimated; both are expected at the package settings and say nothing about the input.
.rp_quiet_fgsea <- function(expr) withCallingHandlers(expr, warning = function(w) {
  if (grepl("P-values were likely overestimated|in reality P-values are less than", conditionMessage(w))) invokeRestart("muffleWarning")
})
# Build a page, dropping table rows until it fits the page.
.rp_fit_page <- function(fun, o, n = NULL, n_min = 4) {
  n_min <- min(n_min %||% 4, n %||% 4)
  limit <- .rp_L$H - .rp_L$foot
  if (is.null(n)) return(fun(o))
  for (k in seq(n, n_min)) {
    b <- fun(o, k)
    if (sum(vapply(b, `[[`, 0, "h")) <= limit) return(b)
  }
  b
}
.rp_error_page <- function(o, key, err) c(.rp_page_head(o, key, "This page could not be drawn."),
  list(.rp_box(list(.rp_text_blk(paste("Error:", conditionMessage(err)), size = 7, width = .rp_inner(.rp_L$W), after = 0.4)))))

#' One-call L1000PAI report from a differential-expression table
#'
#' Scores a differential-expression table against the LINCS L1000 CRISPR knockout and compound
#' libraries and writes three files: a PDF report of up to 19 pages written for readers who are not
#' specialists, an Excel workbook with every table behind the report, and an RDS file with the R
#' objects. Each page of the report answers one question in plain words and states the answer for
#' your data at the top.
#'
#' @details
#' **Input.** One row per gene, as written by limma (`topTable()`), DESeq2 (`results()`) or edgeR
#' (`topTags()`). The gene column is found by its name (`gene`, `symbol`, `SYMBOL`, `gene_name`,
#' `ID`, `ensembl_gene_id`, `entrez`, ...) or taken from the row names, and may hold gene symbols,
#' Ensembl gene IDs or Entrez IDs; the last two are mapped to symbols with org.Hs.eg.db. Genes are
#' ranked by a signed test statistic (`t`, `stat`, ...); without one, by sign(log fold change) times
#' -log10(P value) (`logFC` or `log2FoldChange` with `P.Value`, `pvalue` or `PValue`). When a gene
#' occurs more than once, the row with the largest absolute statistic is kept. Give the full table of
#' tested genes, not only the significant ones.
#'
#' **Scoring.** Both libraries are scored with [run_l1000_pai()] in paired mode (fgsea, gene sets of
#' 10-500 genes, eps 1e-10, observed-tail P values, Benjamini-Hochberg adjustment over pairs; a call
#' needs |activity score| >= 1, a coherent adjusted P value < 0.05 and both arms aligned).
#' Pathways come from [pai_enrich()] (over-representation among the Activated and the Inhibited
#' targets) or, when no term is significant, from [pai_gsea_pathways()]. Compound mechanisms come
#' from the Chem taxonomy 3.0.1 ([get_l1000_annotation_catalog()]) and drug families from
#' [pai_familyplot()].
#'
#' **Pages.** Overview, key findings, quality checks and how to read the calls; input check and
#' proliferation signal; gene knockouts; evidence of the top gene calls ([pai_evidence_card()] and the
#' leading-edge network [pai_cnetplot()]); pathways; pathway enrichment in depth; compounds; evidence
#' of the top compound calls (the same two displays); drug families and mechanism classes; flags and
#' notes; methods, versions and references. The pathway enrichment pages depend on the pathway result: with terms
#' over-represented among the calls, nine pages: [pai_enrichplot()]; [pai_gobubble()] with
#' [pai_term_cnetplot()]; [pai_gocircle()] (top 10, all members, labelled); [pai_gochord()] with
#' [pai_govenn()]; [pai_gobar()] with [pai_gocluster()]; [pai_goheat()]; [pai_armmap()]; and, for each
#' direction, [pai_hubplot()] with [pai_emapplot()], then [pai_chordplot()] with [pai_enrich_cnetplot()]. With only the ranked test,
#' the two hub pages from [pai_gsea_pathways()]; without pathways, none. `enrichment = "core"` keeps
#' the first two and the hub genes.
#'
#' **Checks.** Each check is OK, Check or Caution: gene coverage of the signature genes (>= 90%,
#' >= 70%, lower); number of CRISPR calls (>= 10, 1-9, 0); share of Discordant pairs in either
#' library (<= 10%, <= 30%, higher); a proliferation signal in the input (Hallmark E2F targets or
#' G2M checkpoint with adjusted P < 0.05 and |NES| >= 2 gives Check); common hits among the 20
#' strongest compound calls (<= 3, 4-6, more); called genes absent from the input among the 20
#' strongest gene calls (<= 2, 3-5, more); and whether pathways are over-represented among the calls
#' (OK) or only the ranked test is available (Check).
#'
#' **Common hits.** A compound is a common hit when it is called at the end where it usually sits
#' in unrelated experiments: its median percentile over 12 public GEO perturbation contrasts (9 groups;
#' listed in `system.file("extdata", "report", "README.md", package = "L1000PAI")`) is >= 0.80 and it
#' is Concordant here, or <= 0.20 and it is Opposing.
#'
#' **First run.** The L1000 libraries are downloaded and verified through L1000PAIdata (about 32 MB)
#' and the KEGG and Reactome gene-set libraries are downloaded from Enrichr into the cache; GO gene
#' sets are built from org.Hs.eg.db. Later runs work offline from the caches. A run takes 4-5
#' minutes, most of it scoring, and needs about 6 GB of free memory.
#'
#' The Hallmark E2F targets and G2M checkpoint gene sets shipped with the package are from MSigDB
#' (CC BY 4.0; Liberzon et al., Cell Systems 2015).
#'
#' @param deg A data.frame of differential-expression results, one row per gene.
#' @param contrast A short label of the comparison, shown on every page (for example
#'   `"Treated vs control, liver"`).
#' @param out_dir Folder for the output files; created when needed.
#' @param file_stem Base name of the output files. By default built from `contrast`.
#' @param gene_col,stat_col Column of the gene identifiers and of the ranking statistic, when the
#'   automatic choice should be overridden.
#' @param id_type Type of gene identifier: `"auto"` (default), `"symbol"`, `"ensembl"` or `"entrez"`.
#' @param collections Pathway collections of the pathway page; see [pai_enrich()].
#' @param enrichment Pathway-enrichment pages after the pathway page: `"full"` (default) draws every
#'   enrichment display of the package; `"core"` the overview, the top terms and the hub genes;
#'   `"none"` none.
#' @param excel Write the Excel workbook (default `TRUE`).
#' @param seed Seed of every random step, so that a report can be repeated exactly.
#' @param results Optional `list(crispr = , chem = )` of `L1000PAIResult` objects already computed
#'   from the same table with [run_l1000_pai()] at the settings above, to skip scoring.
#' @param input_label Name of the input shown in the report, for example the file name.
#' @param cache_dir Folder of the gene-set cache; see [get_pai_genesets()].
#' @param download Download (and verify) the L1000 libraries through L1000PAIdata when they are not
#'   yet on this machine. With `FALSE`, a missing library is an error.
#' @param verbose Print progress messages.
#'
#' @return Invisibly, a `pai_deg_report` object: a list with `files` (paths of the PDF, workbook and
#'   RDS), `checks` (the quality checks), `findings` (the sentences of the report), `tables` (input,
#'   CRISPR targets, compounds, pathways, mechanism classes, drug families, proliferation, flags),
#'   `results` (the CRISPR and Chem results and the pathway analyses) and `params`. The same object
#'   is saved in the RDS file.
#' @seealso [run_l1000_pai()], [pai_report()] for a one-page figure of a single library.
#' @family report
#' @export
#' @examples
#' \dontrun{
#' # limma
#' tt <- limma::topTable(fit, coef = "TreatedVsControl", number = Inf)
#' pai_deg_report(tt, contrast = "Treated vs control", out_dir = "L1000PAI_report")
#'
#' # DESeq2 results with Ensembl IDs as row names
#' res <- as.data.frame(DESeq2::results(dds))
#' pai_deg_report(res, contrast = "Knockdown vs control", out_dir = "L1000PAI_report")
#' }
pai_deg_report <- function(deg, contrast = "Contrast", out_dir = ".", file_stem = NULL, gene_col = NULL, stat_col = NULL,
                           id_type = c("auto", "symbol", "ensembl", "entrez"),
                           collections = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome"), enrichment = c("full", "core", "none"),
                           excel = TRUE,
                           seed = 20260813L, results = NULL, input_label = NULL, cache_dir = .pai_cache_dir(),
                           download = TRUE, verbose = TRUE) {
  id_type <- match.arg(id_type); enrichment <- match.arg(enrichment)
  if (!is.data.frame(deg) && !is.matrix(deg)) stop("`deg` must be a data.frame of differential-expression results.", call. = FALSE)
  if (!is.null(results) && (!is.list(results) || !all(c("crispr", "chem") %in% names(results))))
    stop("`results` must be list(crispr = , chem = ) of L1000PAIResult objects.", call. = FALSE)
  if (is.null(results) && !requireNamespace("L1000PAIdata", quietly = TRUE))
    stop("pai_deg_report() scores the L1000 libraries provided by the L1000PAIdata package. Install it with ",
         "remotes::install_github(\"yilin1107/L1000PAIdata\").", call. = FALSE)
  say <- function(...) if (isTRUE(verbose)) message(sprintf(...))
  t0 <- Sys.time()
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  file_stem <- file_stem %||% paste0(sub("_+$", "", substr(gsub("^_+|_+$", "", gsub("[^A-Za-z0-9]+", "_", contrast)), 1, 40)), "_L1000PAI_report")
  files <- stats::setNames(file.path(out_dir, paste0(file_stem, c(".pdf", ".xlsx", "_objects.rds"))), c("pdf", "xlsx", "rds"))
  if (!isTRUE(excel)) files <- files[c("pdf", "rds")]
  op <- options(pai.base_size = 7); on.exit(options(op), add = TRUE)
  genv <- globalenv()   # leave the caller's random-number state as it was
  rng_old <- if (exists(".Random.seed", envir = genv, inherits = FALSE)) get(".Random.seed", envir = genv) else NULL
  on.exit(if (!is.null(rng_old)) assign(".Random.seed", rng_old, envir = genv) else
    if (exists(".Random.seed", envir = genv, inherits = FALSE)) rm(".Random.seed", envir = genv), add = TRUE)
  old_pal <- pai_use_palette("house"); on.exit(pai_use_palette(old_pal %||% "house"), add = TRUE)

  o <- new.env(parent = emptyenv())
  o$prep <- .deg_prepare(deg, gene_col, stat_col, id_type)
  n_used <- nrow(o$prep$rank_input)
  if (n_used < 2000) warning(sprintf("Only %d genes could be used; give the full table of tested genes for reliable scores.", n_used), call. = FALSE)
  say("Read %s genes, ranked by %s.", .fmt_n(n_used), o$prep$stat_method)

  run <- function(lib, label) {
    say("Scoring the %s library (about 1-2 minutes) ...", label)
    .rp_quiet_fgsea(run_l1000_pai(x = o$prep$rank_input, gene_col = "gene", stat_col = "stat", library = lib,
      data_package = "L1000PAIdata", library_mode = "paired", score_direction = "auto", exponent = 1, min_gs_size = 10,
      max_gs_size = 500, eps = 1e-10, pvalue_cutoff = 1, p_adjust_method = "BH", seed = TRUE, by = "fgsea", verbose = FALSE,
      combine_p_method = "max", activity_threshold = 1, classify_p_value_threshold = 0.05, collect_annotation_tables = FALSE,
      rng_seed = seed, keep_all_terms = TRUE, complete_prespecified_terms = TRUE, pvalue_semantics = "observed_tail",
      pair_p_adjust_method = "BH", classify_score_col = "activity_score", classify_require_coherent = TRUE))
  }
  if (is.null(results)) {
    have <- tryCatch(names(L1000PAIdata::list_l1000pai_data_files()), error = function(e) character())
    need <- setdiff(c("crispr", "chem"), have)
    if (length(need)) {
      if (!isTRUE(download)) stop(sprintf("The L1000 %s library is not on this machine; run L1000PAIdata::download_l1000pai_data(c(%s)) or use download = TRUE.",
        paste(need, collapse = " and "), paste0('"', need, '"', collapse = ", ")), call. = FALSE)
      say("Downloading and verifying the L1000 %s (once per machine) ...",
          if (length(need) == 2) "CRISPR and compound libraries" else sprintf("%s library", need))
      L1000PAIdata::download_l1000pai_data(need, quiet = !isTRUE(verbose))
    }
    o$x <- run("crispr", "CRISPR knockout"); o$xc <- run("chem", "compound")
  } else {
    o$x <- results$crispr; o$xc <- results$chem
    for (r in list(o$x, o$xc)) if (!methods::is(r, "L1000PAIResult")) stop("`results` must hold L1000PAIResult objects.", call. = FALSE)
    if (!identical(.pai_library(o$x), "crispr") || !identical(.pai_library(o$xc), "chem"))
      stop("`results$crispr` must come from the crispr library and `results$chem` from the chem library.", call. = FALSE)
    ng <- pai_metadata(o$x)$n_ranked_genes
    if (!is.null(ng) && !identical(as.integer(ng), as.integer(n_used)))
      warning(sprintf("`results` were computed from %s genes but this table gives %s; were they made from the same table?",
                      .fmt_n(ng), .fmt_n(n_used)), call. = FALSE)
  }

  say("Pathways ...")
  o$e <- tryCatch(pai_enrich(o$x, collections = collections, cache_dir = cache_dir), error = function(e) e)
  o$g <- NULL; o$pathway_error <- NULL
  if (inherits(o$e, "error")) {
    o$pathway_error <- conditionMessage(o$e); o$e <- NULL; o$pathway_mode <- "none"
  } else if (any(.as_tt(o$e)$table$p.adjust < 0.05)) {
    o$pathway_mode <- "ora"
  } else {
    o$g <- tryCatch(.rp_quiet_fgsea(pai_gsea_pathways(o$x, collections = collections, seed = seed, cache_dir = cache_dir)),
                    error = function(e) e)
    if (inherits(o$g, "error")) { o$pathway_error <- conditionMessage(o$g); o$g <- NULL; o$pathway_mode <- "none" } else o$pathway_mode <- "gsea"
  }
  if (o$pathway_mode == "ora") {
    t <- o$e$table; t <- t[order(t$p.adjust), ]
    o$pw <- data.frame(collection = t$collection, id = t$ID, term = t$Description, direction = t$direction, padj = t$p.adjust,
                       count = t$k, genes = t$geneID, stringsAsFactors = FALSE)
  } else if (o$pathway_mode == "gsea") {
    t <- o$g$table; t <- t[order(t$p.adjust, -abs(t$NES)), ]
    o$pw <- data.frame(collection = t$collection, id = t$ID, term = t$Description, end = t$end, nes = t$NES, pvalue = t$pvalue,
                       padj = t$p.adjust, size = t$size, genes = t$core_enrichment, stringsAsFactors = FALSE)
    o$pw$end_label <- ifelse(o$pw$end == "Activated", "Activated end", "KO-like end")
  } else {
    o$pw <- data.frame(collection = character(), term = character(), padj = numeric(), stringsAsFactors = FALSE)
  }

  o$page_keys <- .rp_page_keys(o$pathway_mode, enrichment)
  o$hubs <- if (o$pathway_mode %in% c("ora", "gsea")) {
    src <- if (o$pathway_mode == "ora") o$e else o$g
    stats::setNames(lapply(c("Inhibited", "Activated"), function(k)
      tryCatch(pai_enrich_hubs(src, direction = k), error = function(e) conditionMessage(e))), c("Inhibited", "Activated"))
  } else list()

  say("Checks and tables ...")
  cov <- function(lib) tryCatch(mean(unique(get_l1000_library(lib)$gene) %in% o$prep$rank_input$gene), error = function(e) NA_real_)
  o$coverage <- c(crispr = cov("crispr"), chem = cov("chem"))
  o$prolif <- .deg_proliferation(o$prep$rank_input, .rp_prolif_sets(), seed)
  o$tabs <- list(crispr = .deg_target_table(o$x, o$prep), chem = .deg_compound_table(o$xc, .rp_reference()))
  o$classes <- .deg_class_table(o$tabs$chem)
  o$checks <- .deg_checks(o)
  o$text <- .deg_text(o)
  o$flags <- .deg_flags(o)
  o$families <- tryCatch(as.data.frame(attr(suppressMessages(pai_familyplot(o$xc)), "plot_data")), error = function(e) data.frame())

  vers <- vapply(c("L1000PAI", "L1000PAIdata", "fgsea", "org.Hs.eg.db"), function(p)
    tryCatch(as.character(utils::packageVersion(p)), error = function(e) "not installed"), "")
  src <- tryCatch(L1000PAIdata::l1000pai_data_sources(), error = function(e) NULL)
  sha <- if (is.null(src)) c(crispr = "not available", chem = "not available") else stats::setNames(src$content_sha256, src$library_id)
  tax <- get_chem_family_catalog("term_annotation")$taxonomy_version[1] %||% "3.0.1"
  prolif_version <- attr(.rp_prolif_sets(), "db_version")
  o$meta <- list(contrast = contrast, seed = seed, collections = collections, versions = vers, taxonomy = tax, sha = sha,
    reference_version = "v2 (12 contrasts, 9 groups)", msigdb = paste("MSigDB", prolif_version), n_pages = length(o$page_keys),
    byline = sprintf("Input: %s   \u00b7   generated %s   \u00b7   L1000PAI %s", input_label %||% "data.frame", format(Sys.Date()), vers[["L1000PAI"]]),
    footer = sprintf("L1000PAI %s  \u00b7  Chem taxonomy %s  \u00b7  %s  \u00b7  Scores measure similarity to L1000 cell-line signatures, not causation.",
                     vers[["L1000PAI"]], tax, format(Sys.Date())),
    files = files)
  o$params <- data.frame(setting = c("contrast", "input", "gene column", "identifier type", "ranking statistic", "genes used", "seed",
                                     "library mode", "gene-set size", "eps", "p value semantics", "pair adjustment", "activity threshold",
                                     "call P threshold", "pathway collections", "pathway view", "Chem taxonomy", "proliferation sets",
                                     "common-hit reference", paste("version", names(vers)), "CRISPR library sha256", "Chem library sha256"),
                         value = c(contrast, input_label %||% "data.frame", o$prep$gene_col, o$prep$id_type, o$prep$stat_method,
                                   n_used, seed, "paired", "10-500", "1e-10", "observed_tail", "BH", 1, 0.05,
                                   paste(collections, collapse = ", "), o$pathway_mode, tax, o$meta$msigdb,
                                   o$meta$reference_version, vers, sha[["crispr"]], sha[["chem"]]), stringsAsFactors = FALSE)

  say("Writing the PDF ...")
  cairo <- isTRUE(capabilities("cairo"))
  if (!cairo) warning("This R has no cairo support; the PDF is written with pdf() and some symbols may not show.", call. = FALSE)
  prev <- grDevices::dev.cur()
  o$figs <- list(); heights <- stats::setNames(numeric(length(o$page_keys)), o$page_keys)
  if (cairo) grDevices::cairo_pdf(files[["pdf"]], width = .rp_L$page_w / 25.4, height = .rp_L$page_h / 25.4, onefile = TRUE, family = .pai_font())
  else grDevices::pdf(files[["pdf"]], width = .rp_L$page_w / 25.4, height = .rp_L$page_h / 25.4, onefile = TRUE, useDingbats = FALSE)
  tryCatch({
    for (k in o$page_keys) {
      grid::grid.newpage()   # start the page before building it: text is measured on the open page
      spec <- .rp_specs[[k]]
      b <- tryCatch(.rp_fit_page(get(spec$fun, mode = "function"), o, spec$n, spec$n_min), error = function(e) {
        warning(sprintf("Report page %d could not be drawn: %s", .rp_pn(o, k), conditionMessage(e)), call. = FALSE)
        .rp_error_page(o, k, e)
      })
      heights[[k]] <- .rp_draw_page(b, o, .rp_pn(o, k))
    }
  }, finally = {
    grDevices::dev.off()
    if (prev > 1) grDevices::dev.set(prev)
  })
  if (isTRUE(excel)) { say("Writing the Excel workbook ..."); .deg_write_xlsx(o, files[["xlsx"]]) }

  files[] <- file.path(normalizePath(out_dir, winslash = "/"), basename(files))   # one separator on Windows too
  out <- structure(list(
    files = files, contrast = contrast, checks = o$checks, findings = o$text,
    tables = list(input = o$prep$qc, crispr_targets = o$tabs$crispr, chem_compounds = o$tabs$chem, pathways = o$pw,
                  chem_classes = o$classes, chem_families = o$families, proliferation = o$prolif, flags = o$flags),
    results = list(crispr = o$x, chem = o$xc, enrichment = o$e, pathway_gsea = o$g, hubs = o$hubs),
    params = o$params,
    pages = data.frame(page = seq_along(heights), key = o$page_keys, question = vapply(o$page_keys, function(k) .rp_specs[[k]]$q, ""),
                       height_mm = round(unname(heights), 1),
                       fits = unname(heights) <= .rp_L$H - .rp_L$foot)), class = "pai_deg_report")
  saveRDS(out, files[["rds"]])
  say("Done in %.0f s:\n  %s", as.numeric(difftime(Sys.time(), t0, units = "secs")), paste(files, collapse = "\n  "))
  invisible(out)
}

#' @export
print.pai_deg_report <- function(x, ...) {
  cat(sprintf("<pai_deg_report>  %s\n", x$contrast))
  for (f in names(x$files)) cat(sprintf("  %-5s %s\n", f, x$files[[f]]))
  cat("Checks:\n")
  for (i in seq_len(nrow(x$checks))) cat(sprintf("  %-8s %-21s %s\n", x$checks$status[i], x$checks$check[i], x$checks$result[i]))
  invisible(x)
}
