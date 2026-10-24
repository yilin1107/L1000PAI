test_that("available_l1000_annotation_catalogs returns expected names", {
  expect_equal(
    available_l1000_annotation_catalogs(),
    c("chem", "ligand", "kinase_gpcr")
  )
})

test_that("get_l1000_annotation_catalog loads packaged catalogs", {
  chem <- get_l1000_annotation_catalog("chem")
  ligand <- get_l1000_annotation_catalog("ligand")
  kg <- get_l1000_annotation_catalog("kinase_gpcr")

  expect_true(all(c("compound_name", "BroadClass", "FamilyClass", "TargetClass", "EvidenceTier") %in% names(chem)))
  legacy <- get_l1000_annotation_catalog("chem", version = "legacy")
  expect_true(all(c("compound_name", "BigClass", "MechanismClass", "TargetFamily") %in% names(legacy)))
  expect_equal(nrow(chem), nrow(legacy))
  expect_true(all(c("base_ligand", "LigandGroup", "FunctionalTheme") %in% names(ligand)))
  expect_true(all(c("target", "SignalFamily", "TargetClass") %in% names(kg)))
  expect_true(nrow(chem) > 0)
  expect_true(nrow(ligand) > 0)
  expect_true(nrow(kg) > 0)
})

test_that("annotate_l1000_result merges chem annotations", {
  chem_df <- data.frame(
    term = "Erlotinib",
    perturbagen = "Erlotinib",
    activity_score = 2,
    is_significant = TRUE,
    activity_class = "Activated",
    stringsAsFactors = FALSE
  )

  res <- annotate_l1000_result(chem_df, library = "chem")

  expect_true(all(c("compound_name", "BroadClass", "FamilyClass", "TargetClass") %in% names(res)))
  expect_equal(res$compound_name[[1]], "Erlotinib")
  expect_match(res$TargetGenes[[1]], "EGFR")
  expect_identical(res$EvidenceTier[[1]], "A")

  legacy <- annotate_l1000_result(chem_df, library = "chem", version = "legacy")
  expect_true(all(c("BigClass", "MechanismClass", "TargetFamily") %in% names(legacy)))
  expect_false(is.na(legacy$BigClass[[1]]))
})

test_that("build_l1000_annotation_tables returns ligand weighted summaries", {
  ligand_df <- data.frame(
    term = c("AREG(A375)", "AREG(HA1E)", "VIP(A375)"),
    perturbagen = c("AREG(A375)", "AREG(HA1E)", "VIP(A375)"),
    activity_score = c(2.5, 1.9, -1.7),
    p_combined = c(0.001, 0.01, 0.02),
    padj_combined = c(0.01, 0.03, 0.04),
    is_significant = c(TRUE, TRUE, TRUE),
    activity_class = c("Activated", "Activated", "Inhibited"),
    stringsAsFactors = FALSE
  )

  tables <- build_l1000_annotation_tables(ligand_df, library = "ligand")

  expect_true(all(c("annotated_result", "base_ligand_summary", "ligand_group_summary", "functional_theme_summary") %in% names(tables)))
  expect_true("base_ligand" %in% names(tables$base_ligand_summary))
  expect_true("LigandGroup" %in% names(tables$base_ligand_summary))
  expect_true(nrow(tables$ligand_group_summary) > 0)
})

test_that("annotation summaries prefer GLS scores and coherent adjusted p-values", {
  ligand_df <- data.frame(
    term = c("AREG(A375)", "AREG(HA1E)"),
    perturbagen = c("AREG(A375)", "AREG(HA1E)"),
    activity_score = c(-100, 100),
    activity_gls_z = c(3, -1),
    p_combined = c(0.9, 0.001),
    padj_combined = c(0.9, 0.001),
    p_gls_global = c(0.001, 0.2),
    padj_gls_global = c(0.01, 0.2),
    p_coherent = c(0.001, 0.2),
    padj_coherent = c(0.01, 0.2),
    activity_class = c("Activated", "Inhibited"),
    stringsAsFactors = FALSE
  )

  tables <- build_l1000_annotation_tables(ligand_df, library = "ligand")
  areg <- tables$base_ligand_summary[
    tables$base_ligand_summary$base_ligand == "AREG",
    ,
    drop = FALSE
  ]

  expect_equal(nrow(areg), 1L)
  expect_equal(areg$n_sig, 1L)
  expect_equal(areg$n_activated_sig, 1L)
  expect_equal(areg$n_inhibited_sig, 0L)
  expect_equal(areg$supporting_direction, "Activated")
  expect_equal(areg$activity_score_source, "activity_gls_z")
  expect_equal(areg$weighted_mean_activity, 1)
  expect_equal(areg$meta_z, sqrt(2), tolerance = 1e-12)
})

test_that("build_l1000_annotation_tables returns kinase_gpcr summaries", {
  kg_df <- data.frame(
    term = c("ABL1(A375)", "ABL1(HA1E)", "ADORA2A(A375)"),
    perturbagen = c("ABL1(A375)", "ABL1(HA1E)", "ADORA2A(A375)"),
    activity_score = c(2.2, 1.8, -1.4),
    p_combined = c(0.002, 0.02, 0.03),
    padj_combined = c(0.01, 0.03, 0.04),
    is_significant = c(TRUE, TRUE, TRUE),
    activity_class = c("Activated", "Activated", "Inhibited"),
    stringsAsFactors = FALSE
  )

  tables <- build_l1000_annotation_tables(kg_df, library = "kinase_gpcr")

  expect_true(all(c("annotated_result", "target_summary", "signal_family_summary", "target_class_summary", "fine_category_summary") %in% names(tables)))
  expect_true("FineCategory" %in% names(tables$target_summary))
  expect_true(nrow(tables$signal_family_summary) > 0)
})
