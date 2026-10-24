test_that("Chem taxonomy 3.0.1 is the default and passes verification", {
  expect_identical(formals(get_chem_family_catalog)$version, "3.0.1")
  expect_identical(formals(annotate_chem_families)$version, "3.0.1")
  qa <- verify_chem_family_catalogs(
    verify_hash = requireNamespace("digest", quietly = TRUE),
    stop_on_error = FALSE
  )
  expect_true(all(qa$pass), info = paste(qa$check_id[!qa$pass], collapse = ", "))
})

test_that("per-term v3 calls carry targets, tiers and record identifiers", {
  ann <- get_chem_family_catalog("term_annotation")
  expect_equal(nrow(ann), 5425L)
  expect_identical(anyDuplicated(ann$library_term), 0L)
  expect_true(all(ann$tier %in% c("A", "B", "C", "D", "U", "X")))
  call <- function(term) ann[ann$library_term == term, , drop = FALSE]

  expect_identical(call("Gefitinib")$tier, "A")
  expect_true("EGFR" %in% strsplit(call("Gefitinib")$primary_genes, "|", fixed = TRUE)[[1]])
  expect_true("MDM2" %in% strsplit(call("Nutlin-3")$primary_genes, "|", fixed = TRUE)[[1]])
  expect_true("HDAC1" %in% strsplit(call("Vorinostat")$primary_genes, "|", fixed = TRUE)[[1]])
  expect_match(call("Gefitinib")$evidence_record_ids, "ChEMBL mec_id")
  expect_identical(call("GALR1 GALANIN")$tier, "X")
  # legacy labels are kept for comparison, not used for the call
  expect_identical(call("Dacomitinib")$legacy_Class1, "PI3K / AKT / mTOR pathway inhibitor")
  expect_identical(call("Dacomitinib")$legacy_vs_v3, "conflict")
  # class columns only for tier A/B calls
  expect_true(all(is.na(ann$class_family[!ann$tier %in% c("A", "B")])))
})

test_that("evidence table keeps identifiers but no LINCS source text", {
  ev <- get_chem_family_catalog("evidence")
  expect_true(all(ev$source_family %in% c("AU", "CM", "GP", "LM", "BA")))
  lm <- ev$source_family == "LM"
  expect_true(any(lm))
  expect_true(all(is.na(ev$label[lm])))
  expect_true(all(is.na(ev$detail[lm])))
})

test_that("default family tiers keep curated v2 members and corroborated additions", {
  input <- data.frame(
    perturbagen = c("Nutlin-3", "RG7112", "CGM-097", "RITA", "AZD-8931", "Gefitinib", "Docetaxel", "Vinblastine"),
    activity_score = c(3, 2, 2, 1, 1, 1, -1, -1),
    stringsAsFactors = FALSE
  )
  fam <- function(x, term) x$family_id[x$compound_name == term]
  default <- annotate_chem_families(input, family_level = "direct_moa", include_unmapped = FALSE)
  with_c <- annotate_chem_families(input, family_level = "direct_moa", evidence_tier = c("A", "B", "C"), include_unmapped = FALSE)

  expect_true("MDM2_ANTAGONIST" %in% fam(default, "Nutlin-3"))
  expect_true("MDM2_ANTAGONIST" %in% fam(default, "RG7112"))      # curated v2 member, retained (B)
  expect_true("MDM2_ANTAGONIST" %in% fam(default, "CGM-097"))
  expect_false("MDM2_ANTAGONIST" %in% fam(default, "RITA"))       # single-source candidate (C)
  expect_true("MDM2_ANTAGONIST" %in% fam(with_c, "RITA"))
  expect_true("ERBB_TKI_CORE" %in% fam(default, "AZD-8931"))      # new, corroborated (A)
  expect_true("ERBB_TKI_CORE" %in% fam(default, "Gefitinib"))
  expect_true("MICROTUBULE_STABILIZER" %in% fam(default, "Docetaxel"))
  expect_true("MICROTUBULE_TARGETING" %in% fam(default, "Docetaxel"))
  expect_true("MICROTUBULE_TARGETING" %in% fam(default, "Vinblastine"))
})

test_that("every v2 direct family keeps its curated members at the default tiers", {
  v2 <- get_chem_family_catalog("membership", version = "2.0.0")
  v3 <- get_chem_family_catalog("membership")
  v2d <- v2[v2$membership_role == "primary", c("analysis_unit_id", "family_id")]
  v3d <- v3[v3$membership_role == "primary" & v3$evidence_tier %in% c("A", "B"), c("analysis_unit_id", "family_id")]
  key <- function(d) paste(d$analysis_unit_id, d$family_id)
  # a v2 member can drop out only when a corroborated v3 call contradicts it (tier U)
  dropped <- setdiff(key(v2d), key(v3d))
  u <- v3[v3$evidence_tier == "U", , drop = FALSE]
  expect_true(all(dropped %in% key(u)))
})

test_that("Chem annotation tables use the v3 class levels by default", {
  df <- data.frame(
    term = c("Gefitinib", "Erlotinib", "Vorinostat"),
    perturbagen = c("Gefitinib", "Erlotinib", "Vorinostat"),
    activity_score = c(2, 1.5, -2),
    is_significant = TRUE,
    activity_class = c("Activated", "Activated", "Inhibited"),
    stringsAsFactors = FALSE
  )
  tables <- build_l1000_annotation_tables(df, library = "chem")
  expect_true(all(c("annotated_result", "broad_class_summary", "family_class_summary", "target_class_summary") %in% names(tables)))
  expect_true(any(grepl("Kinase", tables$broad_class_summary$category)))
  legacy <- build_l1000_annotation_tables(df, library = "chem", version = "legacy")
  expect_true("big_class_summary" %in% names(legacy))
})

test_that("author decisions override the database evidence", {
  log <- get_chem_family_catalog("curation_log")
  expect_true("Mubritinib" %in% log$library_term)
  ann <- get_chem_family_catalog("term_annotation")
  mub <- ann[ann$library_term == "Mubritinib", , drop = FALSE]
  expect_identical(mub$primary_source, "AU")
  expect_identical(mub$confidence, "author_curated")
  expect_match(mub$class_target, "complex I")
  expect_false("ERBB2" %in% strsplit(mub$primary_genes, "|", fixed = TRUE)[[1]])
  expect_true("MT-ND1" %in% strsplit(mub$primary_genes, "|", fixed = TRUE)[[1]])
  fam <- annotate_chem_families(data.frame(perturbagen = "Mubritinib", activity_score = 1),
    family_level = "all", evidence_tier = c("A", "B", "C"), include_unmapped = FALSE)
  expect_false("ERBB_TKI_CORE" %in% fam$family_id)
  mem <- get_chem_family_catalog("membership")
  unit <- ann$analysis_unit_id[ann$library_term == "Mubritinib"]
  removed <- mem[mem$analysis_unit_id == unit & mem$family_id == "ERBB_TKI_CORE", , drop = FALSE]
  expect_identical(removed$evidence_tier, "U")
  expect_identical(removed$curation_status, "author_override_excluded")
  # rotenone: LINCS writes MT-ND1 as ND1, an alias HGNC also lists for IVNS1ABP
  expect_identical(ann$primary_genes[ann$library_term == "Rotenone"], "MT-ND1")
})
