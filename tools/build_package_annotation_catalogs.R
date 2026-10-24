`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(normalizePath(".", winslash = "/", mustWork = TRUE))
  }
  normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = TRUE)
}

repo_dir <- normalizePath(
  file.path(dirname(get_script_path()), ".."),
  winslash = "/",
  mustWork = TRUE
)

# Maintainer tool: rebuilds inst/extdata/annotation_catalogs from the curated
# annotation tables, which are not part of the repository. Point
# L1000PAI_ANNOTATION_SOURCE at the folder that holds
# Chem_default_library_annotation/ and the ligand and kinase/GPCR tables.
annotation_dir <- Sys.getenv("L1000PAI_ANNOTATION_SOURCE", unset = "")
if (!nzchar(annotation_dir) || !dir.exists(annotation_dir)) {
  stop("Set L1000PAI_ANNOTATION_SOURCE to the curated annotation source folder.")
}
out_dir <- file.path(repo_dir, "inst", "extdata", "annotation_catalogs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_gz_csv <- function(df, path) {
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  utils::write.csv(df, con, row.names = FALSE, na = "")
}

chem <- utils::read.csv(
  file.path(
    annotation_dir,
    "Chem_default_library_annotation",
    "chem_default_library_mechanism_classes_hierarchical_v2.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

chem_keep <- c(
  "compound_name",
  "Class1",
  "Class2",
  "Class3",
  "BigClass",
  "MechanismClass",
  "TargetFamily",
  "class_source",
  "literature_citation",
  "literature_url",
  "cfde_suggested_subtype",
  "cfde_confidence"
)

chem_catalog <- unique(chem[, chem_keep, drop = FALSE])
write_gz_csv(chem_catalog, file.path(out_dir, "chem_annotation_catalog.csv.gz"))

ligand <- utils::read.csv(
  file.path(
    annotation_dir,
    "Ligand_library_annotation",
    "l1000_ligand_library_base_ligand_annotated.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

ligand_keep <- c(
  "base_ligand",
  "n_exact_ligand_entries",
  "n_terms",
  "n_cell_lines",
  "LigandGroup",
  "LigandSubclass",
  "ReceptorAxis",
  "SignalFamily",
  "AnnotationConfidence",
  "FunctionalTheme",
  "FunctionalPolarity"
)

ligand_catalog <- unique(ligand[, ligand_keep, drop = FALSE])
write_gz_csv(ligand_catalog, file.path(out_dir, "ligand_annotation_catalog.csv.gz"))

kinase_gpcr <- utils::read.csv(
  file.path(
    annotation_dir,
    "Kinase_GPCR_library_annotation",
    "l1000_kinase_gpcr_library_targets_annotated.csv"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

kg_keep <- c(
  "target",
  "n_terms",
  "n_cell_lines",
  "TargetClass",
  "TargetSubclass",
  "SignalFamily",
  "FunctionalTheme",
  "AnnotationConfidence",
  "KinaseFamily",
  "GPCRFamily"
)

kg_catalog <- unique(kinase_gpcr[, kg_keep, drop = FALSE])
write_gz_csv(kg_catalog, file.path(out_dir, "kinase_gpcr_annotation_catalog.csv.gz"))

message("Wrote annotation catalogs to: ", out_dir)
