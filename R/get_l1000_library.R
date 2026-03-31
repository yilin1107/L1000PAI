#' List supported bundled L1000 libraries
#'
#' @return A character vector of supported library names.
#' @export
available_l1000_libraries <- function() {
  c(
    "crispr",
    "chem",
    "ligand",
    "oe",
    "single_drug",
    "kinase_gpcr",
    "kinase_gpcr_up",
    "kinase_gpcr_down"
  )
}

#' Get a bundled L1000 library from L1000PAIdata
#'
#' Load a bundled L1000 perturbation signature library distributed in the
#' companion data package `L1000PAIdata`.
#'
#' @param library One of `available_l1000_libraries()`.
#' @param package Name of the companion data package. Default is
#'   `"L1000PAIdata"`.
#'
#' @return A two-column data.frame with columns `term` and `gene`.
#' @export
get_l1000_library <- function(
    library = c(
      "crispr",
      "chem",
      "ligand",
      "oe",
      "single_drug",
      "kinase_gpcr",
      "kinase_gpcr_up",
      "kinase_gpcr_down"
    ),
    package = "L1000PAIdata"
) {
  library <- match.arg(library)

  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Package `", package, "` is required but not installed. ",
      "Please install it first."
    )
  }

  supported <- available_l1000_libraries()
  if (!library %in% supported) {
    stop(
      "`library` must be one of: ",
      paste(supported, collapse = ", ")
    )
  }

  if (library == "kinase_gpcr") {
    up_path <- system.file(
      "extdata",
      "L1000_Kinase_and_GPCR_Perturbations_up_v2.txt",
      package = package
    )
    down_path <- system.file(
      "extdata",
      "L1000_Kinase_and_GPCR_Perturbations_down_v2.txt",
      package = package
    )

    if (identical(up_path, "") || identical(down_path, "")) {
      stop("Bundled kinase/gpcr library files were not found in package `", package, "`.")
    }

    up_tbl <- read_l1000_gmt(up_path)
    down_tbl <- read_l1000_gmt(down_path)
    out <- rbind(up_tbl, down_tbl)
    rownames(out) <- NULL
    return(out)
  }

  fname <- switch(
    library,
    crispr = "LINCS_L1000_CRISPR_KO_Consensus_Sigs.txt",
    chem = "LINCS_L1000_Chem_Pert_Consensus_Sigs.txt",
    ligand = "LINCS_L1000_Ligand_Perturbations_All.txt",
    oe = "l1000_oe.gmt",
    single_drug = "single_drug_perturbations_v1_0.gmt",
    kinase_gpcr_up = "L1000_Kinase_and_GPCR_Perturbations_up_v2.txt",
    kinase_gpcr_down = "L1000_Kinase_and_GPCR_Perturbations_down_v2.txt"
  )

  path <- system.file("extdata", fname, package = package)

  if (identical(path, "")) {
    stop(
      "Bundled file for `", library, "` was not found in package `", package, "`."
    )
  }

  read_l1000_gmt(path)
}
