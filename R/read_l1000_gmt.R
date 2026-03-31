#' Read an L1000 GMT file
#'
#' Read a GMT file and standardize it into a two-column TERM2GENE table
#' with columns `term` and `gene`.
#'
#' @param gmt_file Path to a GMT file.
#'
#' @return A data.frame with two columns: `term` and `gene`.
#' @export
read_l1000_gmt <- function(gmt_file) {
  if (!is.character(gmt_file) || length(gmt_file) != 1) {
    stop("`gmt_file` must be a single character string.")
  }

  if (!file.exists(gmt_file)) {
    stop("`gmt_file` does not exist: ", gmt_file)
  }

  lines <- readLines(gmt_file, warn = FALSE, encoding = "UTF-8")

  if (length(lines) == 0) {
    stop("`gmt_file` is empty.")
  }

  parsed_list <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]

    # GMT 至少應有: term, description, gene1
    if (length(fields) < 3) {
      return(NULL)
    }

    term <- fields[1]
    genes <- fields[-c(1, 2)]

    genes <- genes[!is.na(genes) & genes != ""]

    if (length(genes) == 0) {
      return(NULL)
    }

    data.frame(
      term = rep(term, length(genes)),
      gene = genes,
      stringsAsFactors = FALSE
    )
  })

  parsed_list <- parsed_list[!vapply(parsed_list, is.null, logical(1))]

  if (length(parsed_list) == 0) {
    stop("No valid term-gene pairs were found in `gmt_file`.")
  }

  gmt <- do.call(rbind, parsed_list)
  gmt <- unique(gmt)
  rownames(gmt) <- NULL

  gmt
}
