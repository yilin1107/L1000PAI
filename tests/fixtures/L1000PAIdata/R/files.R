#' List synthetic companion fixture files
#'
#' @return A named character vector.
#' @export
list_l1000pai_data_files <- function() {
  filenames <- c(
    crispr = "LINCS_L1000_CRISPR_KO_Consensus_Sigs.txt",
    chem = "LINCS_L1000_Chem_Pert_Consensus_Sigs.txt",
    ligand = "L1000_Ligand_Consensus_Sigs.gmt",
    oe = "L1000_OverExpression_Consensus_Sigs.gmt",
    single_drug = "single_drug_perturbations_v1_0.gmt",
    kinase_gpcr = "L1000_Kinase_and_GPCR_Consensus_Sigs.txt",
    kinase_gpcr_up = "L1000_Kinase_and_GPCR_Perturbations_up.txt",
    kinase_gpcr_down = "L1000_Kinase_and_GPCR_Perturbations_down.txt"
  )
  paths <- file.path(system.file("extdata", package = "L1000PAIdata"), filenames)
  structure(paths, names = names(filenames))
}
