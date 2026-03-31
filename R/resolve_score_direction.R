#' Resolve paired-score direction from library type
#'
#' @param library Optional bundled library name.
#' @param score_direction One of `"auto"`, `"down_minus_up"`, or
#'   `"up_minus_down"`.
#'
#' @return A resolved score-direction string.
#' @export
resolve_score_direction <- function(
    library = NULL,
    score_direction = c("auto", "down_minus_up", "up_minus_down")
) {
  score_direction <- match.arg(score_direction)

  if (score_direction != "auto") {
    return(score_direction)
  }

  if (is.null(library)) {
    # 對 user-supplied TERM2GENE，保留舊版預設
    return("down_minus_up")
  }

  # perturbation concordance-oriented libraries
  if (library %in% c(
    "single_drug",
    "chem",
    "ligand",
    "oe",
    "kinase_gpcr",
    "kinase_gpcr_up",
    "kinase_gpcr_down"
  )) {
    return("up_minus_down")
  }

  # target activity inference-oriented library
  if (library %in% c("crispr")) {
    return("down_minus_up")
  }

  # fallback
  "down_minus_up"
}
