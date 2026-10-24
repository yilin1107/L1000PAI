# Tests that read the L1000 libraries run only when L1000PAIdata has them in
# its local cache (see L1000PAIdata::download_l1000pai_data()).
skip_if_no_l1000_library <- function(library) {
  testthat::skip_if_not_installed("L1000PAIdata")
  # Absence may skip an integration test; corrupted content must fail it.
  files <- L1000PAIdata::list_l1000pai_data_files()
  testthat::skip_if_not(
    library %in% names(files),
    paste0("L1000 library '", library, "' is not available locally")
  )
}
