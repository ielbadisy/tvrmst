#' tvrmst: Time-Varying RMST from Survival Matrices
#'
#' Model-agnostic utilities for RMST and time-varying RMST quantities computed
#' from survival curves provided on a time grid. Functions accept a time vector
#' and survival matrices, returning RMST-based quantities and bootstrap summaries.
#'
#' @keywords internal
"_PACKAGE"

utils::globalVariables(c(
  ".data", "arm", "estimate", "lower", "s", "survival", "time", "upper", "value"
))
