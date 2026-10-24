#' List supported L1000 libraries
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

#' Get an L1000 library from L1000PAIdata
#'
#' Load an L1000 perturbation signature library provided by the companion
#' package `L1000PAIdata`. Since `L1000PAIdata` 1.0.0 the libraries are not
#' bundled: download and verify them once with
#' `L1000PAIdata::download_l1000pai_data()`.
#'
#' @param library One of `available_l1000_libraries()`.
#' @param package Name of the companion data package. Default is
#'   `"L1000PAIdata"`. Version 1.0.0 or newer is required for the default
#'   provider contract.
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

  if (
    identical(package, "L1000PAIdata") &&
      utils::packageVersion(package) < package_version("1.0.0")
  ) {
    stop(
      "Installed package `L1000PAIdata` is too old. ",
      "Version 1.0.0 or newer is required; found ",
      as.character(utils::packageVersion(package)), "."
    )
  }

  supported <- available_l1000_libraries()
  if (!library %in% supported) {
    stop(
      "`library` must be one of: ",
      paste(supported, collapse = ", ")
    )
  }

  provider <- tryCatch(
    getExportedValue(package, "list_l1000pai_data_files"),
    error = function(e) NULL
  )
  path_provider <- tryCatch(
    getExportedValue(package, "get_l1000pai_data_path"),
    error = function(e) NULL
  )
  path <- ""
  if (is.function(path_provider)) {
    # The current companion validates the requested cached content before
    # returning its path. Do not swallow integrity errors or fall back to an
    # unchecked file. The list-only contract remains for synthetic fixtures.
    path <- path_provider(library)
  } else if (is.function(provider)) {
    files <- provider()
    if (is.character(files) && library %in% names(files)) {
      path <- unname(files[[library]])
    }
  } else if (!identical(package, "L1000PAIdata")) {
    legacy_filename <- switch(
      library,
      crispr = "LINCS_L1000_CRISPR_KO_Consensus_Sigs.txt",
      chem = "LINCS_L1000_Chem_Pert_Consensus_Sigs.txt",
      ligand = "L1000_Ligand_Consensus_Sigs.gmt",
      oe = "L1000_OverExpression_Consensus_Sigs.gmt",
      single_drug = "single_drug_perturbations_v1_0.gmt",
      kinase_gpcr = "L1000_Kinase_and_GPCR_Consensus_Sigs.txt",
      kinase_gpcr_up = "L1000_Kinase_and_GPCR_Perturbations_up.txt",
      kinase_gpcr_down = "L1000_Kinase_and_GPCR_Perturbations_down.txt"
    )
    path <- system.file("extdata", legacy_filename, package = package)
  }

  if (length(path) != 1L || identical(path, "") || !file.exists(path)) {
    hint <- if (identical(package, "L1000PAIdata")) {
      paste0(
        " Download and verify it with L1000PAIdata::download_l1000pai_data(\"",
        library, "\"), or add a local copy with ",
        "L1000PAIdata::register_l1000pai_data_file()."
      )
    } else {
      ""
    }
    stop(
      "The companion provider did not supply a readable file for library `",
      library, "` from package `", package, "`.", hint
    )
  }

  read_l1000_gmt(path)
}
