#' Strip direction suffix from perturbation terms
#'
#' Remove common directional suffixes such as `-up`, `-dn`, `_Up`, `_Down`,
#' ` Up`, or ` Down` from perturbation terms and return the base perturbagen
#' or drug name.
#'
#' @param x A character vector of perturbation terms.
#'
#' @return A character vector with directional suffixes removed.
#' @export
strip_direction_suffix <- function(x) {
  x <- as.character(x)
  x <- trimws(x)

  # single-drug / chem style: -up / -dn
  x <- sub("(?i)-up$", "", x, perl = TRUE)
  x <- sub("(?i)-dn$", "", x, perl = TRUE)

  # paired library style: _Up / _Down / _Dn / -Up / -Down / -Dn / " Up" / " Down" / " Dn"
  x <- sub("([[:space:]_.-]+)(?i:up)$", "", x, perl = TRUE)
  x <- sub("([[:space:]_.-]+)(?i:down)$", "", x, perl = TRUE)
  x <- sub("([[:space:]_.-]+)(?i:dn)$", "", x, perl = TRUE)

  trimws(x)
}


#' Add a Cluster column to L1000 result tables
#'
#' Create a `Cluster` column for downstream plotting or row-binding workflows.
#' For single-signature results, `term` is typically used. For paired results,
#' `perturbagen` is usually preferred when available.
#'
#' @param df A data.frame containing L1000 result rows.
#' @param term_col Name of the term column. Default is `"term"`.
#' @param perturbagen_col Name of the perturbagen column. Default is `"perturbagen"`.
#' @param cluster_col Name of the output cluster column. Default is `"Cluster"`.
#' @param prefer_perturbagen Logical; if `TRUE` and `perturbagen_col` exists,
#'   use it preferentially. Default is `TRUE`.
#'
#' @return A data.frame with an added `Cluster` column.
#' @export
add_cluster_column <- function(
    df,
    term_col = "term",
    perturbagen_col = "perturbagen",
    cluster_col = "Cluster",
    prefer_perturbagen = TRUE
) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }

  source_col <- NULL

  if (prefer_perturbagen && perturbagen_col %in% colnames(df)) {
    source_col <- perturbagen_col
  } else if (term_col %in% colnames(df)) {
    source_col <- term_col
  } else if (perturbagen_col %in% colnames(df)) {
    source_col <- perturbagen_col
  } else {
    stop("Neither `term_col` nor `perturbagen_col` was found in `df`.")
  }

  df[[cluster_col]] <- strip_direction_suffix(df[[source_col]])
  df
}


#' Merge multiple L1000 result tables
#'
#' Safely combine multiple L1000 result tables using `dplyr::bind_rows()`.
#' Missing columns are automatically filled with `NA`.
#'
#' @param ... Data.frames to merge.
#' @param .list Optional list of data.frames to merge.
#' @param .id Optional name of an identifier column passed to `bind_rows()`.
#'
#' @return A merged data.frame.
#' @export
merge_l1000_results <- function(..., .list = NULL, .id = NULL) {
  inputs <- list(...)

  if (!is.null(.list)) {
    if (!is.list(.list)) {
      stop("`.list` must be a list of data.frames.")
    }
    inputs <- c(inputs, .list)
  }

  if (length(inputs) == 0) {
    stop("No input data.frames were provided.")
  }

  ok <- vapply(inputs, is.data.frame, logical(1))
  if (!all(ok)) {
    stop("All inputs must be data.frames.")
  }

  dplyr::bind_rows(inputs, .id = .id)
}


#' Prepare result tables for compareCluster-like downstream plots
#'
#' Add a `Cluster` column and optionally merge with additional result tables.
#'
#' @param df A data.frame containing L1000 result rows.
#' @param other Optional another data.frame, or a list of data.frames, to merge.
#' @param term_col Name of the term column. Default is `"term"`.
#' @param perturbagen_col Name of the perturbagen column. Default is `"perturbagen"`.
#' @param cluster_col Name of the output cluster column. Default is `"Cluster"`.
#' @param prefer_perturbagen Logical; if `TRUE` and `perturbagen_col` exists,
#'   use it preferentially. Default is `TRUE`.
#' @param .id Optional identifier column name for merged outputs.
#'
#' @return A data.frame ready for downstream merging/plotting.
#' @export
prepare_cluster_result <- function(
    df,
    other = NULL,
    term_col = "term",
    perturbagen_col = "perturbagen",
    cluster_col = "Cluster",
    prefer_perturbagen = TRUE,
    .id = NULL
) {
  df2 <- add_cluster_column(
    df = df,
    term_col = term_col,
    perturbagen_col = perturbagen_col,
    cluster_col = cluster_col,
    prefer_perturbagen = prefer_perturbagen
  )

  if (is.null(other)) {
    return(df2)
  }

  if (is.data.frame(other)) {
    return(merge_l1000_results(df2, other, .id = .id))
  }

  if (is.list(other)) {
    return(merge_l1000_results(.list = c(list(df2), other), .id = .id))
  }

  stop("`other` must be NULL, a data.frame, or a list of data.frames.")
}
