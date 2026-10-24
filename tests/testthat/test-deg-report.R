test_that("the input reader recognises limma, DESeq2 and edgeR tables", {
  limma <- data.frame(ID = c("TP53", "MDM2", "CDKN1A", "MDM2"), logFC = c(1, -2, 3, 0.5), AveExpr = 5,
                      t = c(2.5, -4, 6, 1), P.Value = c(0.01, 0.001, 1e-5, 0.2), adj.P.Val = 0.05)
  p <- L1000PAI:::.deg_prepare(limma)
  expect_identical(p$gene_col, "ID"); expect_identical(p$stat_col, "t"); expect_identical(p$id_type, "symbol")
  expect_identical(nrow(p$rank_input), 3L)
  expect_equal(p$rank_input$stat[p$rank_input$gene == "MDM2"], -4)     # duplicate: largest |t| kept
  expect_identical(unname(p$counts[["duplicates"]]), 1L)
  expect_true("duplicate gene, smaller |statistic| dropped" %in% p$qc$status)

  deseq <- data.frame(baseMean = 100, log2FoldChange = c(1, -1), lfcSE = 0.2, stat = c(5, -3), pvalue = 1e-4, padj = 1e-3,
                      row.names = c("ENSG00000141510.18", "ENSG00000135679"))
  p <- L1000PAI:::.deg_prepare(deseq)
  expect_identical(p$gene_col, ".row_names"); expect_identical(p$id_type, "ensembl"); expect_identical(p$stat_col, "stat")
  expect_setequal(p$rank_input$gene, c("TP53", "MDM2"))

  edger <- data.frame(logFC = c(2, -1), logCPM = 4, PValue = c(1e-3, 1e-2), FDR = 0.05, row.names = c("7157", "4193"))
  p <- L1000PAI:::.deg_prepare(edger)
  expect_identical(p$id_type, "entrez")
  expect_equal(p$rank_input$stat, c(3, -2))                            # sign(logFC) * -log10(PValue)
  expect_identical(p$rank_input$gene, c("TP53", "MDM2"))

  expect_error(L1000PAI:::.deg_prepare(data.frame(gene = "TP53", score = 1)), "No signed statistic")
})

test_that("quality checks follow the documented thresholds", {
  o <- new.env()
  o$coverage <- c(crispr = 0.95, chem = 0.95); o$pathway_mode <- "gsea"
  o$prolif <- data.frame(gene_set = c("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT"), nes = c(-2.5, -1), padj = c(0.001, 0.2),
                         strong = c(TRUE, FALSE))
  calls <- c(rep("Activated", 6), rep("Inhibited (KO-like)", 6), rep("Discordant", 20), rep("Not called", 68))
  o$tabs <- list(
    crispr = data.frame(call = calls, activity_score = seq(5, -5, length.out = 100), target_measured = c(rep(FALSE, 4), rep(TRUE, 96))),
    chem = data.frame(call = calls, activity_score = seq(5, -5, length.out = 100), common_hit = c(rep(TRUE, 7), rep(FALSE, 93))))
  ch <- L1000PAI:::.deg_checks(o)
  st <- stats::setNames(ch$status, ch$check)
  expect_identical(st[["Gene coverage"]], "OK")
  expect_identical(st[["Formal calls"]], "OK")            # 12 calls
  expect_identical(st[["Discordant share"]], "Check")     # 20%
  expect_identical(st[["Proliferation signal"]], "Check")
  expect_identical(st[["Common hits"]], "Caution")        # 7 of the 20 strongest compound calls
  expect_identical(st[["Unmeasured targets"]], "Check")   # 4 of the 20 strongest gene calls
  expect_identical(st[["Pathways"]], "Check")
  o$coverage[] <- NA
  expect_identical(L1000PAI:::.deg_checks(o)$status[1], "Check")
})

test_that("the shipped reference files match their documentation", {
  d <- system.file("extdata", "report", package = "L1000PAI")
  ref <- utils::read.csv(file.path(d, "chem_common_response_reference_v2.csv"), stringsAsFactors = FALSE)
  expect_identical(nrow(ref), 5425L)
  expect_identical(anyDuplicated(ref$compound), 0L)
  expect_identical(sum(ref$common_response == "usually_top"), 193L)
  expect_identical(sum(ref$common_response == "usually_bottom"), 138L)
  expect_true(all(ref$n_groups == 9L))
  expect_true(all(ref$reference_percentile[ref$common_response == "usually_top"] >= 0.8))
  sets <- L1000PAI:::.rp_prolif_sets()
  expect_identical(lengths(sets)[["HALLMARK_E2F_TARGETS"]], 200L)
  expect_identical(lengths(sets)[["HALLMARK_G2M_CHECKPOINT"]], 200L)
  readme <- readLines(file.path(d, "README.md"))
  for (f in c("chem_common_response_reference_v2.csv", "hallmark_proliferation_sets.csv"))
    expect_true(any(grepl(digest::digest(file = file.path(d, f), algo = "sha256"), readme, fixed = TRUE)), info = f)
})

test_that("the workbook has no relationship to a missing part", {
  f <- tempfile(fileext = ".xlsx")
  wb <- openxlsx::createWorkbook(); openxlsx::addWorksheet(wb, "a"); openxlsx::writeData(wb, "a", data.frame(x = 1:3))
  openxlsx::saveWorkbook(wb, f)
  L1000PAI:::.rp_xlsx_drop_dangling_rels(f)
  td <- tempfile(); utils::unzip(f, exdir = td)
  for (r in list.files(file.path(td, "xl", "worksheets", "_rels"), full.names = TRUE)) {
    x <- paste(readLines(r, warn = FALSE), collapse = "")
    tg <- regmatches(x, gregexpr('Target="[^"]+"', x))[[1]]
    tg <- sub('Target="([^"]+)"', "\\1", tg)
    expect_true(all(file.exists(file.path(td, "xl", "worksheets", tg))), info = basename(r))
  }
  expect_equal(openxlsx::read.xlsx(f, "a")$x, 1:3)
})

test_that("pai_deg_report() writes the PDF, the workbook and the objects", {
  fx <- toy_plot_fixture()
  crispr <- fx$x; crispr@metadata$library <- "crispr"
  chem <- fx$x; chem@metadata$library <- "chem"     # the toy library stands in for both libraries
  wd <- tempfile("report_"); dir.create(wd); old <- setwd(wd)
  tryCatch(expect_warning(                            # a relative out_dir, as most users give it
    r <- pai_deg_report(fx$ranks, contrast = "Toy contrast", out_dir = "out", results = list(crispr = crispr, chem = chem),
                        collections = c(toy = fx$gmt), verbose = FALSE),
    "Only 1500 genes"), finally = setwd(old))
  expect_identical(unname(r$files[["pdf"]]), normalizePath(file.path(wd, "out", "Toy_contrast_L1000PAI_report.pdf"), winslash = "/"))
  expect_s3_class(r, "pai_deg_report")
  expect_true(all(file.exists(r$files)))
  expect_identical(basename(r$files[["pdf"]]), "Toy_contrast_L1000PAI_report.pdf")
  view <- r$params$value[r$params$setting == "pathway view"]
  expect_identical(view, "ora")                             # the toy pathways are over-represented among the calls
  expect_identical(r$pages$key, c("overview", "input", "crispr", "evidence", "pathways", "enrich_overview", "enrich_terms",
    "enrich_circle", "enrich_chord", "enrich_bars", "enrich_heat", "enrich_armmap", "enrich_hubs", "enrich_hubnet",
    "chem", "chem_evidence", "classes", "cautions", "methods"))
  expect_true(all(r$pages$fits), info = paste(r$pages$height_mm, collapse = ", "))
  expect_identical(readBin(r$files[["pdf"]], "raw", 5), charToRaw("%PDF-"))
  expect_setequal(openxlsx::getSheetNames(r$files[["xlsx"]]),
    c("README", "Summary", "Checks", "Input_QC", "Proliferation", "CRISPR_targets", "Pathways", "Pathway_terms", "Hub_genes", "Hub_modules",
      "Chem_compounds", "Chem_classes", "Chem_families", "Flags", "Parameters"))
  expect_true(all(c("Inhibited", "Activated") %in% names(r$results$hubs)))
  expect_identical(nrow(r$checks), 7L)
  expect_identical(nrow(r$tables$crispr_targets), nrow(as.data.frame(result_table(fx$x))))
  expect_true(all(is.na(r$tables$proliferation$nes)))      # toy genes are not in the Hallmark sets
  expect_identical(readRDS(r$files[["rds"]])$checks, r$checks)
  expect_output(print(r), "Checks")
})

test_that("GO gene sets hold the same gene-term pairs as select(GOALL)", {
  cache <- tempfile("go_")
  g <- get_pai_genesets("GO_MF", cache_dir = cache)
  genes <- c("TP53", "MDM2", "CDKN1A", "EGFR")
  ref <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = genes, columns = c("GOALL", "ONTOLOGYALL"),
                                                keytype = "SYMBOL"))
  ref <- unique(ref[!is.na(ref$GOALL) & ref$ONTOLOGYALL == "MF", c("GOALL", "SYMBOL")])
  got <- g$TERM2GENE[g$TERM2GENE$gene %in% genes, ]
  expect_setequal(paste(got$term, got$gene), paste(ref$GOALL, ref$SYMBOL))
})

test_that("the enrichment pages follow the pathway result", {
  expect_identical(L1000PAI:::.rp_page_keys("gsea")[6:7], c("enrich_hubs", "enrich_hubnet"))
  expect_false(any(grepl("^enrich", L1000PAI:::.rp_page_keys("none"))))
  expect_false(any(grepl("^enrich", L1000PAI:::.rp_page_keys("ora", "none"))))
  expect_length(L1000PAI:::.rp_page_keys("ora"), 19L)
  expect_length(L1000PAI:::.rp_page_keys("ora", "core"), 13L)
  expect_length(L1000PAI:::.rp_page_keys("none"), 10L)
  expect_true(all(vapply(L1000PAI:::.rp_specs, function(s) exists(s$fun, envir = asNamespace("L1000PAI")), TRUE)))
})

test_that("patchwork is loaded with the package, so plots that use | and / work in a fresh session", {
  expect_true("patchwork" %in% names(getNamespaceImports("L1000PAI")))
})
