test_that("available_l1000_libraries returns expected names", {
  libs <- available_l1000_libraries()

  expect_type(libs, "character")
  expect_true(all(c("crispr", "chem", "kinase_gpcr") %in% libs))
})

test_that("get_l1000_library works when L1000PAIdata is installed", {
  testthat::skip_if_not_installed("L1000PAIdata")

  res <- get_l1000_library("crispr")

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("term", "gene"))
  expect_true(nrow(res) > 0)
})

test_that("get_l1000_library can merge kinase_gpcr up/down files", {
  testthat::skip_if_not_installed("L1000PAIdata")

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
