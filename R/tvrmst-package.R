#' tvrmst: Matrix-first dynamic RMST utilities
#'
#' Tools for dynamic RMST estimation, RMST deltas, bootstrap uncertainty,
#' and plotting from survival probability matrices evaluated on a common time grid.
#'
#' @name tvrmst
NULL

#' @keywords internal
"_PACKAGE"

.require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function.", call. = FALSE)
  }
}
