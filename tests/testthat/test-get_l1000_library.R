test_that("available_l1000_libraries returns expected names", {
  libs <- available_l1000_libraries()

  expect_type(libs, "character")
  expect_true(all(c("crispr", "chem", "kinase_gpcr") %in% libs))
})

test_that("get_l1000_library works when the crispr library is available", {
  skip_if_no_l1000_library("crispr")

  res <- get_l1000_library("crispr")

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("term", "gene"))
  expect_true(nrow(res) > 0)
})

test_that("get_l1000_library reads the kinase_gpcr consensus file", {
  skip_if_no_l1000_library("kinase_gpcr")

  res <- get_l1000_library("kinase_gpcr")

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("term", "gene"))
  expect_true(nrow(res) > 0)
})

test_that("get_l1000_library errors for unsupported package", {
  expect_error(
    get_l1000_library("crispr", package = "this_package_should_not_exist_123"),
    "required but not installed"
  )
})

test_that("core refuses an unverified cached library instead of reading it", {
  skip_if_not_installed("L1000PAIdata", minimum_version = "1.0.0")
  exports <- getNamespaceExports("L1000PAIdata")
  skip_if_not(all(c("get_l1000pai_data_path", "l1000pai_data_sources") %in% exports),
              "The tiny CI companion fixture does not implement verified downloads")
  cache <- tempfile("l1000pai-corrupt-cache-")
  dir.create(file.path(cache, "libraries"), recursive = TRUE)
  on.exit(unlink(cache, recursive = TRUE), add = TRUE)
  old_options <- options(L1000PAIdata.cache_dir = cache)
  on.exit(options(old_options), add = TRUE)
  source_row <- L1000PAIdata::l1000pai_data_sources()
  name <- source_row$cache_filename[source_row$library_id == "crispr"]
  writeLines("TP53 up\t\tNOT_A_REFERENCE_GENE\t", file.path(cache, "libraries", name))
  expect_error(get_l1000_library("crispr"), "failed content verification")
})
