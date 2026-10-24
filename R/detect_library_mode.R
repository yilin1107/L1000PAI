#' Detect L1000 library mode from term names
#'
#' @param term2gene A two-column TERM2GENE data.frame.
#'
#' @return One of `"paired"` or `"single"`.
#' @export
detect_library_mode <- function(term2gene) {
  if (!is.data.frame(term2gene) || ncol(term2gene) < 2) {
    stop("`term2gene` must be a data.frame with at least two columns.")
  }

  term <- as.character(term2gene[[1]])

  has_up <- any(grepl("([[:space:]_.-]+)(Up|up)$", term))
  has_down <- any(grepl("([[:space:]_.-]+)(Down|down|Dn|dn)$", term))

  if (has_up && has_down) {
    return("paired")
  }

  "single"
}
