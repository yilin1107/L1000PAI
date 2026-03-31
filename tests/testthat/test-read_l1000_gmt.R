test_that("read_l1000_gmt reads a valid GMT file", {
  gmt_file <- tempfile(fileext = ".gmt")

  writeLines(
    c(
      "TERM_A\tna\tGene1\tGene2",
      "TERM_B\tna\tGene3\tGene4"
    ),
    con = gmt_file
  )

  res <- read_l1000_gmt(gmt_file)

  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("term", "gene"))
  expect_equal(nrow(res), 4)
  expect_equal(unique(res$term), c("TERM_A", "TERM_B"))
})

test_that("read_l1000_gmt removes duplicated term-gene pairs", {
  gmt_file <- tempfile(fileext = ".gmt")

  writeLines(
    c(
      "TERM_A\tna\tGene1\tGene1\tGene2"
    ),
    con = gmt_file
  )

  res <- read_l1000_gmt(gmt_file)

  expect_equal(nrow(res), 2)
  expect_true(all(c("Gene1", "Gene2") %in% res$gene))
})

test_that("read_l1000_gmt errors for missing file", {
  expect_error(
    read_l1000_gmt("this_file_should_not_exist_123.gmt"),
    "does not exist"
  )
})

test_that("read_l1000_gmt errors for non-character input", {
  expect_error(
    read_l1000_gmt(123),
    "single character string"
  )
})
