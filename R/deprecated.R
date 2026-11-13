#' Deprecated functions in tvrmst
#'
#' These functions were renamed to glued lowercase names and will be removed
#' in a future release. Use the new names instead.
#'
#' @param ... Arguments passed to the new function.
#' @return The result of the new function.
#' @name tvrmst-deprecated
#' @keywords internal
NULL

#' @rdname tvrmst-deprecated
#' @export
as_survmat <- function(...) {
  .Deprecated("assurvmat")
  assurvmat(...)
}

#' @rdname tvrmst-deprecated
#' @export
as_survprob_matrix <- function(...) {
  .Deprecated("assurvprobmatrix")
  assurvprobmatrix(...)
}

#' @rdname tvrmst-deprecated
#' @export
bind_survmat <- function(...) {
  .Deprecated("bindsurvmat")
  bindsurvmat(...)
}

#' @rdname tvrmst-deprecated
#' @export
boot_rmst_delta <- function(...) {
  .Deprecated("bootrmstdelta")
  bootrmstdelta(...)
}

#' @rdname tvrmst-deprecated
#' @export
bootstrap_curve <- function(...) {
  .Deprecated("bootstrapcurve")
  bootstrapcurve(...)
}

#' @rdname tvrmst-deprecated
#' @export
nobs_survmat <- function(...) {
  .Deprecated("nobssurvmat")
  nobssurvmat(...)
}

#' @rdname tvrmst-deprecated
#' @export
plot_boot_curve <- function(...) {
  .Deprecated("plotbootcurve")
  plotbootcurve(...)
}

#' @rdname tvrmst-deprecated
#' @export
plot_delta_curve <- function(...) {
  .Deprecated("plotdeltacurve")
  plotdeltacurve(...)
}

#' @rdname tvrmst-deprecated
#' @export
plot_rmst_individual_by_group <- function(...) {
  .Deprecated("plotrmstindividualbygroup")
  plotrmstindividualbygroup(...)
}

#' @rdname tvrmst-deprecated
#' @export
plot_rmst_two_arms <- function(...) {
  .Deprecated("plotrmsttwoarms")
  plotrmsttwoarms(...)
}

#' @rdname tvrmst-deprecated
#' @export
rmst_delta <- function(...) {
  .Deprecated("rmstdelta")
  rmstdelta(...)
}

#' @rdname tvrmst-deprecated
#' @export
rmst_dynamic <- function(...) {
  .Deprecated("rmstdynamic")
  rmstdynamic(...)
}
