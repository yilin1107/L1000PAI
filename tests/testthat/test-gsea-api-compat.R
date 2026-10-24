.gsea_api_args <- function() {
  ranks <- stats::setNames(seq(5, -5, length.out = 40), paste0("G", seq_len(40)))
  list(geneList = ranks, exponent = 1, minGSSize = 2, maxGSSize = 20,
       eps = 1e-10, pvalueCutoff = 1, pAdjustMethod = "BH",
       TERM2GENE = data.frame(term = rep(c("POS", "NEG"), each = 6),
                              gene = names(ranks)[c(1:6, 35:40)]),
       verbose = FALSE, seed = FALSE)
}

test_that("legacy GSEA receives the original explicit backend and permutation contract", {
  calls <- list()
  old <- function(by = "fgsea", ...) {
    calls[[length(calls) + 1L]] <<- c(list(by = by), list(...))
    data.frame(ID = "legacy")
  }
  args <- .gsea_api_args()
  expect_equal(.call_l1000_gsea(args, "fgsea", 77, gsea_function = old)$ID, "legacy")
  expect_identical(calls[[1]]$by, "fgsea")
  expect_false("nPerm" %in% names(calls[[1]]))
  expect_equal(calls[[1]][names(args)], args)
  .call_l1000_gsea(args, "DOSE", gsea_function = old)
  expect_equal(calls[[2]]$nPerm, 1000)
  .call_l1000_gsea(args, "DOSE", 77, gsea_function = old)
  expect_equal(calls[[3]]$nPerm, 77)
  expect_identical(calls[[3]]$by, "DOSE")
})

test_that("modern method API never receives by or substitutes the enrichit engine", {
  calls <- NULL
  modern <- function(method = "multilevel", ...) stop("Modern engine must not be called")
  spy <- function(...) {
    calls <<- list(...)
    data.frame(pathway = c("NEG", "POS"), pval = c(0.02, 0.01),
               ES = c(-1, 1), NES = c(-1.2, 1.4), size = c(6L, 6L))
  }
  args <- .gsea_api_args()
  out <- .call_l1000_gsea(args, "fgsea", gsea_function = modern, fgsea_function = spy)
  expect_false(any(c("by", "method", "nperm", "nPerm") %in% names(calls)))
  expect_identical(calls$stats, args$geneList)
  expect_equal(calls[c("minSize", "maxSize", "eps", "gseaParam", "nproc")],
               list(minSize = 2, maxSize = 20, eps = 1e-10, gseaParam = 1, nproc = 0))
  expect_equal(out$ID, c("POS", "NEG"))
  expect_equal(out$p.adjust, c(0.02, 0.02))
  expect_equal(out$NES, c(1.4, -1.2))
  expect_equal(out$rank, c(6L, 7L))
  expect_equal(out$leading_edge, c("tags=100%, list=15%, signal=100%", "tags=100%, list=18%, signal=97%"))
  expect_equal(out$core_enrichment, c(paste0("G", 1:6, collapse = "/"), paste0("G", 36:40, collapse = "/")))
  args$TERM2NAME <- data.frame(term = c("POS", "NEG"), label = c("Positive pathway", NA))
  named <- .call_l1000_gsea(args, "fgsea", gsea_function = modern, fgsea_function = spy)
  expect_equal(named$Description, c("Positive pathway", "NEG"))
})

test_that("modern API explicitly rejects DOSE rather than remapping it", {
  modern <- function(method = "multilevel", ...) stop("Must not be called")
  expect_error(.call_l1000_gsea(.gsea_api_args(), "DOSE", gsea_function = modern),
               'no alternative engine has been substituted')
  expect_error(.call_l1000_gsea(.gsea_api_args(), "fgsea", gsea_function = function(...) NULL),
               'Unrecognized clusterProfiler')
})

test_that("the direct adapter adjusts its full universe before both cutoff filters", {
  args <- .gsea_api_args()
  args$pvalueCutoff <- 0.05
  args$pAdjustMethod <- "bonferroni"
  spy <- function(...) data.frame(pathway = c("POS", "NEG"), pval = c(0.02, 0.04),
                                  ES = c(1, -1), NES = c(2, -3), size = c(6L, 6L))
  out <- .direct_l1000_fgsea(args, spy)
  expect_equal(out$ID, "POS")
  expect_equal(out$p.adjust, 0.04)
  args$pvalueCutoff <- 0.01
  empty <- .direct_l1000_fgsea(args, spy)
  expect_equal(nrow(empty), 0L)
  expect_true(all(c("ID", "NES", "pvalue", "p.adjust", "qvalue", "rank", "leading_edge", "core_enrichment") %in% names(empty)))
})

test_that("qvalue failure is non-estimable and NA P-values are not significant", {
  args <- .gsea_api_args()
  spy <- function(...) data.frame(pathway = c("POS", "NEG"), pval = c(0.01, NA_real_),
                                  ES = c(1, -1), NES = c(2, NA_real_), size = c(6L, 6L))
  out <- .direct_l1000_fgsea(args, spy)
  expect_equal(out$ID, "POS")
  expect_true(is.na(out$qvalue))
  expect_equal(out$pvalue, 0.01)
})

test_that("estimable qvalues preserve the documented legacy estimation settings", {
  args <- .gsea_api_args()
  labels <- paste0("TERM", seq_len(10))
  args$TERM2GENE <- data.frame(term = rep(labels, each = 6),
                              gene = rep(names(args$geneList)[1:6], 10))
  values <- c(0.01, 0.08, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.8, 0.9)
  spy <- function(...) data.frame(pathway = labels, pval = values, ES = 1, NES = 2, size = 6L)
  set.seed(221L)
  expected <- qvalue::qvalue(values, lambda = 0.05, pi0.method = "bootstrap")$qvalues
  set.seed(221L)
  actual <- .direct_l1000_fgsea(args, spy)
  expect_true(all(is.finite(actual$qvalue)))
  expect_equal(actual$qvalue[match(labels, actual$ID)], expected)
  expect_equal(actual$p.adjust[match(labels, actual$ID)], stats::p.adjust(values, "BH"))
})

test_that("the direct fgsea adapter matches the installed legacy wrapper contract", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("fgsea")
  skip_if_not_installed("qvalue")
  if (!"by" %in% names(formals(clusterProfiler::GSEA))) skip("Legacy wrapper unavailable on this runtime")
  args <- .gsea_api_args()
  # Two small pathways make an inexpensive real-engine comparison, including
  # seeded state, adjusted-P universe, sorting and diagnostic membership.
  for (seeded in c(FALSE, TRUE)) {
    args$seed <- seeded
    args$TERM2NAME <- data.frame(term = c("POS", "NEG"), label = c("Positive pathway", "Negative pathway"))
    set.seed(164L)
    legacy <- suppressWarnings(as.data.frame(do.call(clusterProfiler::GSEA, c(args, list(by = "fgsea")))))
    set.seed(164L)
    direct <- suppressWarnings(.direct_l1000_fgsea(args))
    expect_equal(direct, legacy, tolerance = 1e-12)
  }
})

test_that("ranked pathway analysis uses the same modern adapter and retains descriptions", {
  fx <- toy_plot_fixture()
  original <- .call_l1000_gsea
  modern <- function(method = "multilevel", ...) stop("Modern engine must not be called")
  calls <- list()
  spy <- function(pathways, ...) {
    labels <- names(pathways)
    data.frame(pathway = labels, pval = seq_along(labels) / 1000,
      ES = rep(c(0.4, -0.4), length.out = length(labels)),
      NES = rep(c(2, -2), length.out = length(labels)), size = lengths(pathways))
  }
  local_mocked_bindings(.call_l1000_gsea = function(args, by, nPerm = NULL) {
    calls[[length(calls) + 1L]] <<- list(args = args, by = by)
    original(args, by, nPerm, gsea_function = modern, fgsea_function = spy)
  })
  out <- pai_gsea_pathways(fx$x, collections = c(Toy = fx$gmt), min_size = 5,
                           library = "crispr", cache_dir = tempdir())
  expect_s3_class(out, "pai_pathway_gsea")
  expect_identical(calls[[1]]$by, "fgsea")
  expect_true(is.data.frame(calls[[1]]$args$TERM2NAME))
  expect_s3_class(out$results$Toy, "data.frame")
  expected <- get_pai_genesets(fx$gmt, cache_dir = tempdir())$TERM2NAME
  expect_equal(out$table$Description, expected[[2L]][match(out$table$ID, expected[[1L]])])
  expect_true(all(c("NES", "pvalue", "p.adjust", "core_enrichment") %in% names(out$table)))
})
