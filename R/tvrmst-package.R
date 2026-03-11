#' tvrmst: Matrix-first dynamic restricted mean survival time utilities
#'
#' Tools for dynamic restricted mean survival time estimation, restricted mean
#' survival time contrasts, bootstrap uncertainty, and plotting from survival
#' probability matrices evaluated on a common time grid.
#'
#' @name tvrmst
NULL

#' @keywords internal
"_PACKAGE"

if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "delta", "estimate", "group", "hi", "id", "lo", "mean", "t", "tau"
  ))
}

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.", call. = FALSE)
  }
}
