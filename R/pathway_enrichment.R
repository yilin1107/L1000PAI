# Pathway enrichment of L1000PAI calls: over-representation (pai_enrich) and ranked GSEA (pai_gsea_pathways).
#
# The Activated targets and the Inhibited targets of a gene-level library (CRISPR knockout,
# overexpression) are tested separately against GO BP/MF/CC, KEGG, Reactome and the Elsevier
# Pathway Collection. The background is the set of targets that L1000PAI actually scored,
# restricted to genes annotated in the collection, so that the composition of the L1000
# perturbation library does not masquerade as enrichment.

.pai_collections <- c(
  GO_BP = "GO:BP", GO_MF = "GO:MF", GO_CC = "GO:CC",
  KEGG = "KEGG_2026", Reactome = "Reactome_Pathways_2024", Elsevier = "Elsevier_Pathway_Collection"
)
.enrichr_url <- "https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=%s"

.geneset_dir <- function(cache_dir) file.path(cache_dir, "genesets")

# Where gene-set libraries are kept: option L1000PAI.cache_dir, else the L1000PAIdata cache, else
# the user cache directory of L1000PAI.
.pai_cache_dir <- function() {
  d <- getOption("L1000PAI.cache_dir")
  if (!is.null(d)) return(d)
  if (requireNamespace("L1000PAIdata", quietly = TRUE)) return(L1000PAIdata::l1000pai_data_cache_dir())
  tools::R_user_dir("L1000PAI", "cache")
}

# Collections as given by the user: names(.pai_collections), Enrichr library names or GMT files.
# A named vector sets the labels used in tables and plots.
.collection_labels <- function(collections) {
  lab <- names(collections) %||% rep("", length(collections))
  auto <- ifelse(file.exists(collections), tools::file_path_sans_ext(basename(collections)), collections)
  stats::setNames(unname(collections), ifelse(nzchar(lab), lab, auto))
}

.write_gmt <- function(t2g, t2n, path) {
  sp <- split(t2g$gene, t2g$term)
  desc <- t2n$name[match(names(sp), t2n$term)]
  lines <- vapply(seq_along(sp), function(i)
    paste(c(names(sp)[i], desc[i], sort(unique(sp[[i]]))), collapse = "\t"), character(1))
  writeLines(sort(lines), path, useBytes = TRUE)
}

.read_gmt2 <- function(path) {
  lines <- readLines(path, encoding = "UTF-8", warn = FALSE)
  lines <- lines[nzchar(lines)]
  parts <- strsplit(lines, "\t", fixed = TRUE)
  term <- vapply(parts, `[`, character(1), 1L)
  desc <- vapply(parts, function(p) if (length(p) >= 2L && nzchar(p[2L])) p[2L] else p[1L], character(1))
  genes <- lapply(parts, function(p) {
    g <- p[-(1:2)]
    g <- sub(",.*$", "", g)             # Enrichr weighted entries "GENE,1.0"
    unique(g[nzchar(g)])
  })
  list(
    TERM2GENE = data.frame(term = rep(term, lengths(genes)), gene = unlist(genes), stringsAsFactors = FALSE),
    TERM2NAME = data.frame(term = term, name = desc, stringsAsFactors = FALSE)
  )
}

# Gene sets for pathway enrichment
# 
# GO sets are built from org.Hs.eg.db (all evidence codes, propagated to ancestors, as in
# clusterProfiler::enrichGO). KEGG, Reactome and the Elsevier Pathway Collection are
# Enrichr gene-set libraries downloaded on first use; any other Enrichr library name also
# works. Every collection is kept as a GMT file in the cache, with its source, version,
# retrieval time and SHA-256 recorded in genesets/geneset_record.csv.
#' Gene sets for pathway enrichment
#'
#' GO sets are built from org.Hs.eg.db (all evidence codes, propagated to ancestors, as in
#' `clusterProfiler::enrichGO()`). KEGG, Reactome and the Elsevier Pathway Collection are Enrichr
#' gene-set libraries downloaded on first use; any other Enrichr library name also works, and so
#' does the path of a GMT file of your own. Downloaded and built collections are kept as GMT files
#' in the cache, with source, version, retrieval time and SHA-256 recorded in
#' `genesets/geneset_record.csv`.
#'
#' @param collection One of `"GO_BP"`, `"GO_MF"`, `"GO_CC"`, `"KEGG"`, `"Reactome"`,
#'   `"Elsevier"`, an Enrichr library name, or a GMT file.
#' @param cache_dir Directory of the cache: `getOption("L1000PAI.cache_dir")`, else the cache of
#'   L1000PAIdata when it is installed, else `tools::R_user_dir("L1000PAI", "cache")`.
#' @param refresh Download or build again even when the file is cached.
#' @return A list with `TERM2GENE` and `TERM2NAME` data frames; attributes `version`, `file` and
#'   `sha256` record where the sets came from.
#' @family pathway analysis
#' @export
get_pai_genesets <- function(collection, cache_dir = .pai_cache_dir(), refresh = FALSE) {
  if (length(collection) != 1L) stop("`collection` must be a single name or file.", call. = FALSE)
  if (file.exists(collection)) {                                  # a GMT file of the user's own
    out <- .read_gmt2(collection)
    if (!nrow(out$TERM2GENE)) stop("No gene sets in ", collection, call. = FALSE)
    attr(out, "collection") <- collection
    attr(out, "version") <- basename(collection)
    attr(out, "file") <- normalizePath(collection)
    attr(out, "sha256") <- digest::digest(collection, algo = "sha256", file = TRUE)
    return(out)
  }
  src <- if (collection %in% names(.pai_collections)) .pai_collections[[collection]] else collection
  dir <- .geneset_dir(cache_dir)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (startsWith(src, "GO:")) {
    ont <- sub("GO:", "", src)
    go_meta <- AnnotationDbi::metadata(GO.db::GO.db)
    go_date <- go_meta$value[go_meta$name == "GOSOURCEDATE"]
    version <- sprintf("org.Hs.eg.db %s; GO %s", utils::packageVersion("org.Hs.eg.db"), go_date)
    file <- file.path(dir, sprintf("GO_%s__orgHs%s_GO%s.gmt", ont, utils::packageVersion("org.Hs.eg.db"), go_date))
    if (refresh || !file.exists(file)) {
      # The distinct (GO term, symbol) pairs of the propagated annotation, read in SQL from the
      # org.Hs.eg.db tables that select(keytype = "SYMBOL", columns = "GOALL") reads: the same pairs
      # and GMT file, in seconds and under 1 GB instead of about 6 GB of memory.
      m <- DBI::dbGetQuery(AnnotationDbi::dbconn(org.Hs.eg.db::org.Hs.eg.db), sprintf(paste(
        "SELECT DISTINCT a.go_id AS GOALL, g.symbol AS SYMBOL FROM go_%s_all a",
        "JOIN gene_info g ON a._id = g._id"), tolower(ont)))
      t2g <- data.frame(term = m$GOALL, gene = m$SYMBOL, stringsAsFactors = FALSE)
      ids <- unique(t2g$term)
      nm <- AnnotationDbi::Term(ids)
      t2n <- data.frame(term = ids, name = ifelse(is.na(nm), ids, unname(nm)), stringsAsFactors = FALSE)
      .write_gmt(t2g, t2n, file)
      src_desc <- "org.Hs.eg.db GOALL"
    }
  } else {
    file <- file.path(dir, paste0(src, ".gmt"))
    version <- src
    if (refresh || !file.exists(file)) {
      raw <- tempfile(fileext = ".txt")
      st <- utils::download.file(sprintf(.enrichr_url, utils::URLencode(src, reserved = TRUE)), raw,
        mode = "wb", quiet = TRUE)
      if (!identical(as.integer(st), 0L)) stop("Download of Enrichr library ", src, " failed.")
      gs <- .read_gmt2(raw)
      if (!nrow(gs$TERM2GENE)) stop("Enrichr library ", src, " is empty or unknown.")
      .write_gmt(gs$TERM2GENE, gs$TERM2NAME, file)
      src_desc <- sprintf(.enrichr_url, src)
    }
  }
  rec_path <- file.path(dir, "geneset_record.csv")
  if (exists("src_desc", inherits = FALSE)) {
    rec <- data.frame(collection = collection, version = version, file = basename(file), source = src_desc,
      sha256 = digest::digest(file, algo = "sha256", file = TRUE),
      retrieved_utc = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))
    utils::write.table(rec, rec_path, sep = ",", row.names = FALSE, append = file.exists(rec_path),
      col.names = !file.exists(rec_path))
  }
  out <- .read_gmt2(file)
  attr(out, "collection") <- collection
  attr(out, "version") <- version
  attr(out, "file") <- file
  attr(out, "sha256") <- digest::digest(file, algo = "sha256", file = TRUE)
  out
}

# Over-representation of Activated and Inhibited targets
# 
# @param x An L1000PAIResult from a gene-level library (crispr, oe).
# @param collections Names from `names(.pai_collections)` or Enrichr library names.
# @param min_size,max_size Gene-set size limits, counted within the background.
# @return A list with `table` (one row per collection, direction and term), `results`
#   (clusterProfiler enrichResult objects, usable with enrichplot), `background` and `sets`.
#' Pathways over-represented among the Activated and the Inhibited targets
#'
#' The Activated and the Inhibited targets of a gene-level library (CRISPR knockout,
#' overexpression) are tested separately (hypergeometric test through
#' `clusterProfiler::enricher()`). The background is the set of targets L1000PAI actually scored,
#' so that the make-up of the perturbation library does not pass for enrichment. Testing the two
#' directions separately keeps them from cancelling out.
#'
#' The result is the input of every pathway-level plot: [pai_gobubble()], [pai_gocircle()],
#' [pai_armmap()], [pai_ladder()], [pai_enrichplot()], [pai_enrich_hubs()] and the others.
#'
#' @param x An L1000PAIResult from a gene-level library.
#' @param collections Gene-set collections: any of `"GO_BP"`, `"GO_MF"`, `"GO_CC"`, `"KEGG"`,
#'   `"Reactome"`, `"Elsevier"`, other Enrichr library names, or GMT files. Names of a named
#'   vector are used as the collection labels.
#' @param min_size,max_size Gene-set size limits, counted within the background.
#' @param p_adjust_method Multiple-testing correction within each collection and direction.
#' @param library The library of `x` when it cannot be read from the result (for example
#'   `"crispr"` for a knockout library of your own given as `term2gene`).
#' @inheritParams get_pai_genesets
#' @return A `pai_enrichment` object: `table` (one row per collection, direction and term),
#'   `terms` (the term table, see [pai_term_table()]), `results` (clusterProfiler `enrichResult`
#'   objects, usable with enrichplot), `background`, `lists`, `sets` (version, file and SHA-256
#'   of every collection) and the per-target values the plots need.
#' @family pathway analysis
#' @export
pai_enrich <- function(x, collections = names(.pai_collections), min_size = 10, max_size = 500,
                       p_adjust_method = "BH", library = NULL, cache_dir = .pai_cache_dir()) {
  lib <- .pai_library(x, library)
  if (!lib %in% c("crispr", "oe")) {
    stop("pai_enrich() needs a gene-level library (crispr or oe); `", lib, "` targets are not genes. ",
      "For a library of your own whose targets are genes, give library = \"crispr\" or \"oe\".", call. = FALSE)
  }
  sd <- .pai_score_direction(x, lib)
  sem <- .pai_semantics(lib, sd)
  df <- .pai_df(x)
  cols <- .collection_labels(collections)
  scored <- !is.na(df$activity_score)
  gene <- df$perturbagen %||% df$term
  background <- unique(gene[scored])
  lists <- list(
    Activated = unique(gene[scored & as.character(df$activity_class) == "Activated"]),
    Inhibited = unique(gene[scored & as.character(df$activity_class) == "Inhibited"]))
  tabs <- list(); res <- list(); sets <- list()
  for (cl in names(cols)) {
    gs <- get_pai_genesets(cols[[cl]], cache_dir = cache_dir)
    sets[[cl]] <- attributes(gs)[c("version", "file", "sha256")]
    for (d in names(lists)) {
      r <- suppressMessages(clusterProfiler::enricher(
        gene = lists[[d]], universe = background, TERM2GENE = gs$TERM2GENE, TERM2NAME = gs$TERM2NAME,
        minGSSize = min_size, maxGSSize = max_size, pvalueCutoff = 1, qvalueCutoff = 1,
        pAdjustMethod = p_adjust_method))
      res[[paste(cl, d, sep = "|")]] <- r
      t <- if (is.null(r)) data.frame() else as.data.frame(r)
      if (nrow(t)) {
        gr <- do.call(rbind, strsplit(t$GeneRatio, "/")); br <- do.call(rbind, strsplit(t$BgRatio, "/"))
        t <- data.frame(collection = cl, direction = d, direction_label = sem$short[[d]],
          ID = t$ID, Description = t$Description,
          k = as.integer(gr[, 1]), n = as.integer(gr[, 2]), K = as.integer(br[, 1]), N = as.integer(br[, 2]),
          fold_enrichment = (as.integer(gr[, 1]) / as.integer(gr[, 2])) / (as.integer(br[, 1]) / as.integer(br[, 2])),
          pvalue = t$pvalue, p.adjust = t$p.adjust, qvalue = t$qvalue, geneID = t$geneID,
          stringsAsFactors = FALSE)
        tabs[[length(tabs) + 1L]] <- t
      }
    }
  }
  table <- do.call(rbind, tabs)
  if (is.null(table)) stop("No gene set of ", min_size, "-", max_size, " scored targets was tested.", call. = FALSE)
  rownames(table) <- NULL
  e <- structure(list(table = table, results = res, background = background, lists = lists, sets = sets,
                 library = lib, semantics = sem, score_direction = sd, threshold = .pai_threshold(x),
                 targets = .pai_targets(df),
                 params = list(min_size = min_size, max_size = max_size, p_adjust_method = p_adjust_method,
                 cache_dir = cache_dir)), class = "pai_enrichment")
  e$terms <- .build_term_table(e)
  e
}

# The per-target values the pathway plots need, kept with the enrichment so that they take one object.
.pai_targets <- function(df) {
  keep <- intersect(c("term", "activity_score", "activity_class", "NES_Up", "NES_Down", "padj_coherent"), names(df))
  out <- df[, keep, drop = FALSE]
  out$activity_class <- as.character(out$activity_class)
  out$gene <- df$perturbagen %||% df$term
  rownames(out) <- NULL
  out
}

#' @export
print.pai_enrichment <- function(x, ...) {
  t <- x$table; sig <- t[t$p.adjust < 0.05, ]
  cat(sprintf("<pai_enrichment>  %s library, %s scored targets (%s Activated, %s Inhibited)\n", x$library,
    .fmt_n(length(x$background)), .fmt_n(length(x$lists$Activated)), .fmt_n(length(x$lists$Inhibited))))
  cat(sprintf("  collections: %s\n", paste(names(x$sets), collapse = ", ")))
  cat(sprintf("  terms with adj. P < 0.05: %s Activated, %s Inhibited\n",
    .fmt_n(sum(sig$direction == "Activated")), .fmt_n(sum(sig$direction == "Inhibited"))))
  cat("  x$table: one row per collection, direction and term; plot with pai_gobubble(x), pai_enrichplot(x), ...\n")
  invisible(x)
}

.wrap_text <- function(s, width = 42) vapply(s, function(z) paste(strwrap(z, width), collapse = "\n"), character(1))

# Butterfly plot: Activated targets to the left, Inhibited targets to the right
#' Butterfly plot of the direction-specific enrichment
#'
#' For each collection, the top terms of the Activated targets to the left and of the Inhibited
#' targets to the right, at their -log10 adjusted P; dot size is the fold enrichment.
#'
#' @inheritParams pai_plot_params
#' @param n Terms per collection and direction.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_enrichplot <- function(e, n = 4, collections = NULL, cutoff = 0.05, title = NULL) {
  t <- e$table
  if (!is.null(collections)) t <- t[t$collection %in% collections, ]
  cols <- unique(t$collection)
  pick <- unlist(lapply(cols, function(cl) unlist(lapply(c("Activated", "Inhibited"), function(d) {
    s <- t[t$collection == cl & t$direction == d & t$p.adjust < cutoff, ]
    s <- s[order(s$p.adjust, -s$fold_enrichment), ]
    paste(cl, utils::head(s$ID, n), sep = "|")
  }))))
  t$key <- paste(t$collection, t$ID, sep = "|")
  d <- t[t$key %in% pick, ]
  if (!nrow(d)) stop(sprintf(paste0("No term has adj. P < %s (smallest: %s). To draw the top terms descriptively, ",
    "use a larger cutoff (e.g. cutoff = 1); they are then not significant."), cutoff,
    if (nrow(t)) signif(min(t$p.adjust), 2) else "no term tested"), call. = FALSE)
  ord <- do.call(rbind, lapply(split(d, d$key), function(z) data.frame(key = z$key[1],
    best = min(z$p.adjust), dir = z$direction[which.min(z$p.adjust)])))
  ord <- ord[order(ord$dir == "Inhibited", -ord$best), ]
  lab_map <- stats::setNames(.wrap_text(ifelse(nchar(d$Description) > 58, paste0(substr(d$Description, 1, 55), "..."),
    d$Description), 40), d$key)
  pretty <- c(GO_BP = "GO BP", GO_MF = "GO MF", GO_CC = "GO CC", KEGG = "KEGG", Reactome = "Reactome", Elsevier = "Elsevier")
  bgN <- tapply(t$N, t$collection, max)
  strip <- stats::setNames(sprintf("%s\nN = %s", ifelse(cols %in% names(pretty), pretty[cols], cols), .fmt_n(as.integer(bgN[cols]))), cols)
  d$term <- factor(d$key, levels = ord$key)
  d$collection <- factor(d$collection, levels = cols)
  d$sig <- d$p.adjust < cutoff
  d$side <- ifelse(d$direction == "Activated", -1, 1)
  d$x <- d$side * .neglog10(d$p.adjust, 1e-30)
  thr <- -log10(cutoff)
  lim <- max(abs(d$x)) * 1.18
  sizes <- vapply(e$lists, length, integer(1))
  empty <- vapply(c("Activated", "Inhibited"), function(z) !any(d$sig & d$direction == z), logical(1))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = x, y = term)) +
    ggplot2::annotate("rect", xmin = -thr, xmax = thr, ymin = -Inf, ymax = Inf, fill = .tok$pill) +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = x, yend = term, colour = direction), linewidth = 0.55, alpha = 0.5) +
    ggplot2::geom_point(data = d[!d$sig, ], ggplot2::aes(size = fold_enrichment, colour = direction), shape = 21,
      fill = "white", stroke = 0.6) +
    ggplot2::geom_point(data = d[d$sig, ], ggplot2::aes(size = fold_enrichment, fill = direction), shape = 21,
      colour = "white", stroke = 0.55) +
    ggplot2::geom_text(ggplot2::aes(x = x + side * lim * 0.065, label = k, hjust = ifelse(side > 0, 0, 1)),
      size = 2, colour = .tok$ink2, family = .pai_font()) +
    ggplot2::facet_grid(collection ~ ., scales = "free_y", space = "free_y", labeller = ggplot2::as_labeller(strip)) +
    ggplot2::scale_y_discrete(labels = function(k) lab_map[k]) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.2, 3.8), name = "Fold enrichment", breaks = scales::breaks_pretty(3),
      guide = ggplot2::guide_legend(title.position = "left", nrow = 1, override.aes = list(fill = .tok$muted, colour = "white"))) +
    ggplot2::scale_x_continuous(limits = c(-lim, lim), labels = function(v) abs(v), breaks = scales::breaks_pretty(6)) +
    ggplot2::labs(x = "\u2212log10 adjusted P", y = NULL) +
    theme_pai(grid = "x") +
    ggplot2::theme(legend.position = "bottom", legend.justification = "left",
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold", hjust = 0, colour = .tok$ink, lineheight = 1.05),
      axis.text.y = ggplot2::element_text(colour = .tok$ink, size = ggplot2::rel(0.88), lineheight = 0.85),
      axis.ticks.y = ggplot2::element_blank(), panel.spacing.y = ggplot2::unit(2.2, "mm"))
  none <- function(z) if (empty[[z]]) "  (none significant)" else ""
  chips <- .pai_chips(c(stats::setNames(.pal$marks[["Activated"]], sprintf("\u2190 %s targets \u00b7 %s%s",
      e$semantics$short[["Activated"]], .fmt_n(sizes[["Activated"]]), none("Activated"))),
    stats::setNames(.pal$marks[["Inhibited"]], sprintf("%s targets \u00b7 %s%s \u2192",
      e$semantics$short[["Inhibited"]], .fmt_n(sizes[["Inhibited"]]), none("Inhibited"))),
    stats::setNames(.tok$axis, sprintf("open / shaded: adj. P \u2265 %s", cutoff))))
  out <- patchwork::wrap_plots(list(
    .pai_header(title %||% "Pathways over-represented among the called targets",
      sprintf("Background: the %s scored %s targets, kept if annotated in the collection (N).\nNumbers beside the dots: targets in the term.",
        .fmt_n(length(e$background)), .pai_libname(e$library))),
    chips, p,
    .pai_badge("ORA", sprintf("hypergeometric test  \u00b7  %s within each collection and direction  \u00b7  sets of %d\u2013%d genes",
      e$params$p_adjust_method, e$params$min_size, e$params$max_size))), ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(2) + 1, 5.4, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- d
  .with_size(out, 130, round(min(260, 52 + length(unique(d$key)) * 3.3 + length(cols) * 4), 1))
}

# Ranked pathway GSEA on activity scores ----------------------------------------------------
# Pathway GSEA on the activity scores of all scored targets
# 
# No call list is cut: every scored target of a gene-level library is ranked by its activity
# score, and each pathway is tested for concentration at either end (fgsea through
# clusterProfiler::GSEA). Pathways are restricted to the scored targets. A positive NES means
# the pathway's genes sit among the Activated end (activity up); negative, the Inhibited end.
#' Pathway GSEA on the activity scores of all scored targets
#'
#' No call list is cut: every scored target is ranked by its activity score and each pathway is
#' tested for concentration at either end (fgsea through `clusterProfiler::GSEA()`). A positive
#' NES means the pathway's genes sit at the Activated end, a negative NES at the Inhibited end.
#'
#' @inheritParams pai_enrich
#' @param eps,seed Passed to fgsea; the seed makes the P values reproducible.
#' @return A `pai_pathway_gsea` object with `table`, `results`, `stats` (the activity ranking) and
#'   `sets`.
#' @family pathway analysis
#' @export
pai_gsea_pathways <- function(x, collections = names(.pai_collections), min_size = 10, max_size = 500,
                              eps = 1e-10, seed = 20261023L, p_adjust_method = "BH", library = NULL,
                              cache_dir = .pai_cache_dir()) {
  lib <- .pai_library(x, library)
  if (!lib %in% c("crispr", "oe")) {
    stop("pai_gsea_pathways() needs a gene-level library (crispr or oe); `", lib, "` targets are not genes. ",
      "For a library of your own whose targets are genes, give library = \"crispr\" or \"oe\".", call. = FALSE)
  }
  cols <- .collection_labels(collections)
  sem <- .pai_semantics(lib, .pai_score_direction(x, lib))
  df <- .pai_df(x)
  df <- df[!is.na(df$activity_score), ]
  gene <- df$perturbagen %||% df$term
  if (anyDuplicated(gene)) stop("Duplicated target genes in the result table.")
  stats <- sort(stats::setNames(df$activity_score, gene), decreasing = TRUE)
  tabs <- list(); res <- list(); sets <- list()
  for (cl in names(cols)) {
    gs <- get_pai_genesets(cols[[cl]], cache_dir = cache_dir)
    sets[[cl]] <- attributes(gs)[c("version", "file", "sha256")]
    set.seed(seed)
    r <- suppressMessages(suppressWarnings(clusterProfiler::GSEA(stats, TERM2GENE = gs$TERM2GENE,
      TERM2NAME = gs$TERM2NAME, exponent = 1, minGSSize = min_size, maxGSSize = max_size, eps = eps,
      pvalueCutoff = 1, pAdjustMethod = p_adjust_method, seed = TRUE, by = "fgsea", verbose = FALSE)))
    res[[cl]] <- r
    t <- as.data.frame(r)
    if (nrow(t)) {
      tabs[[cl]] <- data.frame(collection = cl, ID = t$ID, Description = t$Description, size = t$setSize,
        ES = t$enrichmentScore, NES = t$NES, pvalue = t$pvalue, p.adjust = t$p.adjust,
        end = ifelse(t$NES > 0, sem$short[["Activated"]], sem$short[["Inhibited"]]),
        core_enrichment = t$core_enrichment, stringsAsFactors = FALSE)
    }
  }
  table <- do.call(rbind, tabs); rownames(table) <- NULL
  structure(list(table = table, results = res, stats = stats, sets = sets, library = lib, semantics = sem,
                 params = list(min_size = min_size, max_size = max_size, eps = eps, seed = seed,
                 p_adjust_method = p_adjust_method, cache_dir = cache_dir)), class = "pai_pathway_gsea")
}

# Pathways at either end of the activity ranking, with the ORA result alongside
#' Pathways at either end of the activity ranking
#'
#' The top pathways of [pai_gsea_pathways()] at each end of the ranking, with a tag on those also
#' found by the over-representation analysis.
#'
#' @inheritParams pai_plot_params
#' @param g Result of [pai_gsea_pathways()].
#' @param ora Optional result of [pai_enrich()], to tag pathways it also finds.
#' @param n Pathways per end and collection.
#' @return A patchwork (or ggplot) object. The values drawn are attached as
#'   `attr(p, "plot_data")` and a recommended size in millimetres as `attr(p, "size_mm")`,
#'   which [pai_save()] uses.
#' @family pathway plots
#' @export
pai_pathwayplot <- function(g, ora = NULL, n = 4, collections = NULL, cutoff = 0.05, title = NULL) {
  t <- g$table
  if (!is.null(collections)) t <- t[t$collection %in% collections, ]
  cols <- unique(t$collection)
  sel <- do.call(rbind, lapply(cols, function(cl) do.call(rbind, lapply(c(1, -1), function(s) {
    z <- t[t$collection == cl & sign(t$NES) == s & t$p.adjust < cutoff, ]
    utils::head(z[order(z$p.adjust, -abs(z$NES)), ], n)
  }))))
  if (is.null(sel) || !nrow(sel)) stop("No pathway passes the cutoff.")
  sel$key <- paste(sel$collection, sel$ID, sep = "|")
  sel$term <- factor(sel$key, levels = sel$key[order(sel$NES)])
  lab_map <- stats::setNames(.wrap_text(ifelse(nchar(sel$Description) > 58, paste0(substr(sel$Description, 1, 55), "..."),
    sel$Description), 40), sel$key)
  pretty <- c(GO_BP = "GO BP", GO_MF = "GO MF", GO_CC = "GO CC", KEGG = "KEGG", Reactome = "Reactome", Elsevier = "Elsevier")
  strip <- stats::setNames(ifelse(cols %in% names(pretty), pretty[cols], cols), cols)
  sel$collection <- factor(sel$collection, levels = cols)
  sel$cls <- ifelse(sel$NES > 0, "Activated", "Inhibited")
  sel$ora <- NA
  if (!is.null(ora)) {
    o <- ora$table
    hit <- paste(o$collection, o$ID, o$direction, sep = "|")[o$p.adjust < cutoff]
    sel$ora <- paste(sel$collection, sel$ID, sel$cls, sep = "|") %in% hit
  }
  lim <- max(abs(sel$NES)) * 1.12
  ox <- lim * 1.2
  p <- ggplot2::ggplot(sel, ggplot2::aes(x = NES, y = term)) +
    ggplot2::geom_vline(xintercept = 0, colour = .tok$axis, linewidth = 0.35) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = term, colour = cls), linewidth = 0.6, alpha = 0.5) +
    ggplot2::geom_point(ggplot2::aes(size = size, fill = cls), shape = 21, colour = "white", stroke = 0.55) +
    ggplot2::facet_grid(collection ~ ., scales = "free_y", space = "free_y", labeller = ggplot2::as_labeller(strip)) +
    ggplot2::scale_y_discrete(labels = function(k) lab_map[k]) +
    ggplot2::scale_colour_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_fill_manual(values = .pal$marks, guide = "none") +
    ggplot2::scale_size_continuous(range = c(1.3, 3.8), name = "Pathway size", breaks = scales::breaks_pretty(3),
      guide = ggplot2::guide_legend(title.position = "left", nrow = 1, override.aes = list(fill = .tok$muted))) +
    ggplot2::scale_x_continuous(limits = c(-lim, if (is.null(ora)) lim else ox + lim * 0.12),
      breaks = function(l) { b <- pretty(c(-lim, lim), 5); b[abs(b) <= lim] }, labels = .lab_minus) +
    ggplot2::labs(x = "NES on the activity ranking", y = NULL) +
    theme_pai(grid = "x") +
    ggplot2::theme(legend.position = "bottom", legend.justification = "left",
      strip.text.y = ggplot2::element_text(angle = 0, face = "bold", hjust = 0, colour = .tok$ink),
      axis.text.y = ggplot2::element_text(colour = .tok$ink, size = ggplot2::rel(0.88), lineheight = 0.85),
      axis.ticks.y = ggplot2::element_blank(), panel.spacing.y = ggplot2::unit(2.2, "mm"))
  if (!is.null(ora)) {
    p <- p + ggplot2::geom_label(data = sel[sel$ora, ], ggplot2::aes(x = ox, label = "ORA"), size = 1.85,
      family = .pai_font(), fontface = "bold", colour = .tok$ink2, fill = .tok$pill, linewidth = 0,
      label.r = grid::unit(0.9, "mm"), label.padding = grid::unit(0.5, "mm")) +
      ggplot2::geom_point(data = sel[!sel$ora, ], ggplot2::aes(x = ox), shape = 16, size = 0.6, colour = .tok$axis)
  }
  chips <- .pai_chips(c(stats::setNames(.pal$marks[["Inhibited"]], sprintf("\u2190 %s end \u00b7 activity down", g$semantics$short[["Inhibited"]])),
    stats::setNames(.pal$marks[["Activated"]], sprintf("%s end \u00b7 activity up \u2192", g$semantics$short[["Activated"]])),
    if (!is.null(ora)) c(`ORA tag: also significant by ORA` = .tok$axis)))
  parts <- list(.pai_header(title %||% "Pathways at either end of the activity ranking",
      sprintf("All %s scored %s targets ranked by activity score; no list is cut. Top %d per end and collection.",
        .fmt_n(length(g$stats)), .pai_libname(g$library), n)), chips, p,
    .pai_badge("GSEA", sprintf("fgsea on the activity ranking  \u00b7  %s within each collection  \u00b7  adj. P < %s",
      g$params$p_adjust_method, cutoff)))
  out <- patchwork::wrap_plots(parts, ncol = 1) +
    patchwork::plot_layout(heights = grid::unit(c(.header_mm(1) + 1, 5.4, 1, 6.5), c("mm", "mm", "null", "mm")))
  attr(out, "plot_data") <- sel
  .with_size(out, 130, round(min(260, 52 + nrow(sel) * 3.3 + length(cols) * 4), 1))
}
