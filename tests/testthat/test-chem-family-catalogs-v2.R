test_that("Chem family v2 catalogs satisfy frozen schema", {
  term <- get_chem_family_catalog("term_registry")
  family <- get_chem_family_catalog("family_registry")
  membership <- get_chem_family_catalog("membership")
  source <- get_chem_family_catalog("source_registry")

  expect_identical(available_chem_family_versions(), "2.0.0")
  expect_equal(nrow(term), 5425L)
  expect_equal(sum(term$curation_status == "mechanism_curated_v2"), 80L)
  expect_identical(anyDuplicated(term$library_term_id), 0L)
  expect_identical(anyDuplicated(term$library_term), 0L)
  expect_identical(anyDuplicated(family$family_id), 0L)
  expect_true(all(term$pair_complete))
  expect_true(all(membership$analysis_unit_id %in% term$analysis_unit_id))
  expect_true(all(membership$family_id %in% family$family_id))
  expect_true(all(membership$source_id %in% source$source_id))

  membership_key <- paste(membership$analysis_unit_id, membership$family_id, sep = "\r")
  expect_identical(anyDuplicated(membership_key), 0L)
})

test_that("Chem family v2 keeps outcome fields outside frozen taxonomy", {
  tables <- lapply(
    c("term_registry", "family_registry", "membership", "source_registry"),
    get_chem_family_catalog
  )
  forbidden <- paste0(
    "(^|_)(activity_score|activity_gls_z|activity_class|nes(_up|_down)?|",
    "pvalue|p_value|p_combined|padj.*|fdr|qvalue|neg_log10.*|evaluated|",
    "recovered|recovery.*|is_significant|n_recovered.*|best_recovered.*)(_|$)"
  )
  for (tbl in tables) {
    expect_length(
      grep(forbidden, names(tbl), ignore.case = TRUE, perl = TRUE, value = TRUE),
      0L
    )
  }
})

test_that("Chem family v2 freezes canonical analysis units", {
  term <- get_chem_family_catalog("term_registry")
  unit <- setNames(term$analysis_unit_id, term$library_term)
  preferred <- setNames(term$preferred_for_analysis, term$library_term)

  expect_identical(unit[["17-AAG"]], unit[["Tanespimycin"]])
  expect_true(preferred[["17-AAG"]])
  expect_false(preferred[["Tanespimycin"]])
  expect_identical(unit[["Ixazomib"]], unit[["MLN-9708"]])
  expect_true(preferred[["Ixazomib"]])
  expect_false(preferred[["MLN-9708"]])
  expect_true("PLX-4720" %in% term$library_term)

  preferred_n <- tapply(as.integer(term$preferred_for_analysis), term$analysis_unit_id, sum)
  expect_true(all(preferred_n == 1L))
})

test_that("direct mechanism and pathway context remain distinct", {
  input <- data.frame(
    perturbagen = c("Erlotinib", "Cobimetinib", "Alisertib", "Docetaxel", "Vandetanib"),
    activity_score = c(2, 1, 3, -1, 0.5),
    stringsAsFactors = FALSE
  )

  direct <- annotate_chem_families(
    input,
    family_level = "direct_moa",
    evidence_tier = c("A", "B")
  )
  context <- annotate_chem_families(
    input,
    family_level = "pathway_context",
    evidence_tier = c("A", "B")
  )

  expect_true("ERBB_TKI_CORE" %in% direct$family_id[direct$compound_name == "Erlotinib"])
  expect_false("MEK1_2_INHIBITOR" %in% direct$family_id[direct$compound_name == "Erlotinib"])
  expect_true("MEK1_2_INHIBITOR" %in% direct$family_id[direct$compound_name == "Cobimetinib"])
  expect_true("AURORA_KINASE_INHIBITOR" %in% direct$family_id[direct$compound_name == "Alisertib"])
  expect_false("MICROTUBULE_STABILIZER" %in% direct$family_id[direct$compound_name == "Alisertib"])
  expect_true("MICROTUBULE_STABILIZER" %in% direct$family_id[direct$compound_name == "Docetaxel"])
  expect_equal(
    direct$mapping_status[direct$compound_name == "Vandetanib"],
    "no_membership_at_requested_level"
  )
  expect_true("RAS_RAF_MEK_ERK_AXIS" %in% context$family_id[context$compound_name == "Vandetanib"])
})

test_that("many-to-many family annotation does not select the first match", {
  input <- data.frame(
    perturbagen = "Erlotinib",
    activity_score = 2,
    stringsAsFactors = FALSE
  )
  annotated <- annotate_chem_families(
    input,
    family_level = "all",
    evidence_tier = c("A", "B")
  )

  expect_gt(nrow(annotated), 1L)
  expect_identical(anyDuplicated(names(annotated)), 0L)
  expect_true(all(c("ERBB_TKI_CORE", "ERBB_RTK_AXIS", "RAS_RAF_MEK_ERK_AXIS") %in% annotated$family_id))
  expect_true(all(annotated$result_row_id == 1L))
})

test_that("packaged Chem family v2 verification passes", {
  qa <- verify_chem_family_catalogs(
    verify_hash = requireNamespace("digest", quietly = TRUE),
    stop_on_error = FALSE
  )
  expect_true(all(qa$pass), info = paste(qa$check_id[!qa$pass], collapse = ", "))
})

test_that("legacy Chem annotation behavior remains one-to-one", {
  legacy <- get_l1000_annotation_catalog("chem")
  expect_equal(nrow(legacy), 5425L)
  expect_identical(anyDuplicated(legacy$compound_name), 0L)

  input <- data.frame(
    perturbagen = "Erlotinib",
    activity_score = 2,
    activity_class = "Activated",
    stringsAsFactors = FALSE
  )
  annotated <- annotate_l1000_result(input, library = "chem")
  expect_equal(nrow(annotated), 1L)
})
