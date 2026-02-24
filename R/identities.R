#' Identity checks for correctness
#'
#' Checks: mu(s+tau) - mu(s) = S(s) * mu_c(s,tau)
#'
#' @param S Survival matrix (n_units x n_time)
#' @param time Time grid
#' @param s_grid Landmark grid
#' @param tau Window length
#' @param eps Stability threshold for mu_c
#' @return data.frame(series, max_abs_error)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' s_grid <- c(0, 1, 2)
#' check_identities(S, t, s_grid, tau = 2)
#' @export
check_identities <- function(S, time, s_grid, tau, eps = 1e-8) {
  S <- .coerce_unit_time(S, time, "S")
  time <- as.numeric(time)

  rd <- rmst_dynamic(S, time)
  unit_names <- setdiff(names(rd), "tau")
  mu_mat <- t(as.matrix(rd[, unit_names, drop = FALSE]))

  mu_s <- .surv_at(time, mu_mat, s_grid)
  mu_s_tau <- .surv_at(time, mu_mat, s_grid + tau)
  S_s <- .surv_at(time, S, s_grid)

  muc <- tvrmst_cond(S, time, s_grid, tau, eps = eps)
  muc_mat <- t(as.matrix(muc[, unit_names, drop = FALSE]))

  lhs <- mu_s_tau - mu_s
  rhs <- S_s * muc_mat

  out <- data.frame(
    series = unit_names,
    max_abs_error = apply(abs(lhs - rhs), 1, max, na.rm = TRUE)
  )

  out
}
