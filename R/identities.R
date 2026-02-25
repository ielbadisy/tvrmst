#' Check RMST identity
#'
#' Checks `RMST(s+tau)-RMST(s) = S(s) * mu_c(s,tau)`.
#'
#' @param S Survival matrix/data.frame (rows = units, cols = time).
#' @param time Numeric strictly increasing time vector.
#' @param s_grid Landmark times.
#' @param tau Positive window length.
#' @param eps Positivity threshold.
#' @param statistic One of `"mean"` or `"unit"`.
#' @return `statistic = "mean"`: data.frame with `s`, `max_abs_error`.
#'   `statistic = "unit"`: data.frame with `s`, `unit`, `max_abs_error`.
#' @export
check_identities <- function(S, time, s_grid, tau, eps = 1e-8,
                             statistic = c("mean", "unit")) {
  S <- validate_surv_input(S, time, name = "S")
  statistic <- match.arg(statistic)
  time <- as.numeric(time)

  if (!is.numeric(s_grid) || length(s_grid) < 1 || any(!is.finite(s_grid)) || any(s_grid < 0)) {
    .stop("`s_grid` must be a finite numeric vector with nonnegative values.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    .stop("`tau` must be a single positive number.")
  }

  rmst_mat <- rmst_dynamic(S, time)
  mu_cond <- tvrmst_cond(S, time, s_grid, tau, eps = eps, statistic = "unit")

  rmst_s <- .linear_at(time, rmst_mat, s_grid)
  rmst_s_tau <- .linear_at(time, rmst_mat, s_grid + tau)
  surv_s <- .surv_at(time, S, s_grid)

  err <- abs((rmst_s_tau - rmst_s) - (surv_s * mu_cond))

  if (statistic == "mean") {
    return(data.frame(s = s_grid, max_abs_error = apply(err, 2, max, na.rm = TRUE)))
  }

  unit_names <- rownames(S)
  if (is.null(unit_names)) unit_names <- paste0("unit", seq_len(nrow(S)))

  data.frame(
    s = rep(s_grid, each = nrow(S)),
    unit = rep(unit_names, times = length(s_grid)),
    max_abs_error = as.vector(err)
  )
}
