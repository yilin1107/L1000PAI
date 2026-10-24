# The plotting functions combine ggplot2 objects with patchwork's `|` and `/` operators. Those
# methods exist only once the patchwork namespace is loaded, so it is imported here: otherwise
# pai_gocircle(), pai_gochord(), pai_goheat() and others fail when called before any other
# patchwork function in a session.
#' @importFrom patchwork wrap_plots
NULL
