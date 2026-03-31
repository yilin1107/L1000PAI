#' Extract the ranked gene vector
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A named numeric vector.
#' @export
rank_table <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }
  x@input_rank
}

#' Extract the raw GSEA table
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A data.frame.
#' @export
gsea_table <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }
  x@gsea_table
}

#' Extract the paired activity table
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A data.frame.
#' @export
pair_table <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }
  x@pair_table
}

#' Extract metadata
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A list.
#' @export
pai_metadata <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }
  x@metadata
}


#' Extract the final result table
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A data.frame.
#' @export
result_table <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }
  x@pair_table
}

#' Extract the library mode used
#'
#' @param x An `L1000PAIResult` object.
#'
#' @return A character scalar.
#' @export
result_mode <- function(x) {
  if (!methods::is(x, "L1000PAIResult")) {
    stop("`x` must be an L1000PAIResult object.")
  }
  x@metadata$library_mode_used
}
