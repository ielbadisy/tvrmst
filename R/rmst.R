#' RMST at a single horizon tau
#'
#' RMST(tau) = \eqn{\int_0^\tau S(u)\,du} computed by trapezoid rule on the grid.
#'
#' @param t Numeric time grid (length n_time).
#' @param S Survival matrix (n_time x n_series).
#' @param tau Nonnegative scalar horizon.
#' @return Numeric vector of length ncol(S).
#' @export
rmst_tau <- function(t, S, tau) {
  .check_survmat(t, S, "S")
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau < 0) {
    .stop("`tau` must be a single nonnegative number.")
  }

  ss <- .sort_survmat(t, S); t <- ss$t; S <- ss$S
  ez <- .extend_to_zero(t, S); t <- ez$t; S <- ez$S

  keep <- t <= tau
  t0 <- t[keep]
  S0 <- S[keep, , drop = FALSE]

  if (length(t0) < 2) {
    # degenerate: tau is before the first provided time > 0
    out <- rep(tau, ncol(S))
    names(out) <- colnames(S)
    return(out)
  }

  et <- .extend_to_tau(t0, S0, tau)
  out <- .trapz_cols(et$t, et$S)
  names(out) <- colnames(S)
  out
}

#' RMST curve evaluated at tau = `t[-1]`
#'
#' Returns a data.frame with tau and RMST(tau) for each series column.
#'
#' @param t Numeric time grid.
#' @param S Survival matrix (n_time x n_series).
#' @return data.frame with column tau and one column per series.
#' @export
rmst_curve <- function(t, S) {
  .check_survmat(t, S, "S")
  ss <- .sort_survmat(t, S); t <- ss$t; S <- ss$S
  ez <- .extend_to_zero(t, S); t <- ez$t; S <- ez$S

  rmst_mat <- .cumtrapz_cols(t, S)          # length(t)-1 x ncol(S)
  out <- data.frame(tau = t[-1], rmst_mat, check.names = FALSE)
  if (!is.null(colnames(S))) names(out)[-1] <- colnames(S)
  out
}

#' Individual (dynamic) RMST curves aligned to the time grid
#'
#' Returns RMST values for each series at every grid time, including tau = 0.
#'
#' @param t Numeric time grid.
#' @param S Survival matrix (n_time x n_series).
#' @return data.frame with column tau (including 0) and one column per series.
#' @export
rmst_dynamic <- function(t, S) {
  .check_survmat(t, S, "S")
  ss <- .sort_survmat(t, S); t <- ss$t; S <- ss$S
  ez <- .extend_to_zero(t, S); t <- ez$t; S <- ez$S

  rmst_mat <- rbind(rep(0, ncol(S)), .cumtrapz_cols(t, S))
  out <- data.frame(tau = t, rmst_mat, check.names = FALSE)
  if (!is.null(colnames(S))) names(out)[-1] <- colnames(S)
  out
}

#' Delta RMST curve between two arms in an rmst_curve() output
#'
#' @param rmst_curve_df Output of rmst_curve().
#' @param arm1 Column name for arm 1.
#' @param arm0 Column name for arm 0.
#' @return data.frame(tau, delta)
#' @export
rmst_delta_curve <- function(rmst_curve_df, arm1, arm0) {
  if (!is.data.frame(rmst_curve_df) || !"tau" %in% names(rmst_curve_df)) {
    .stop("`rmst_curve_df` must come from rmst_curve().")
  }
  if (!arm1 %in% names(rmst_curve_df) || !arm0 %in% names(rmst_curve_df)) {
    .stop("`arm1` and `arm0` must be columns in rmst_curve_df.")
  }
  data.frame(tau = rmst_curve_df$tau, delta = rmst_curve_df[[arm1]] - rmst_curve_df[[arm0]])
}
