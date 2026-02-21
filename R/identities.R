#' Identity checks for correctness
#'
#' Checks: mu(s+tau) - mu(s) = S(s) * mu_c(s,tau)
#'
#' @param t Time grid
#' @param S Survival matrix
#' @param s_grid Landmark grid
#' @param tau Window length
#' @param eps Stability threshold for mu_c
#' @return data.frame(series, max_abs_error)
#' @export
check_identities <- function(t, S, s_grid, tau, eps = 0.05) {
  .check_survmat(t, S, "S")
  ss <- .sort_survmat(t, S); t <- ss$t; S <- ss$S
  ez <- .extend_to_zero(t, S); t <- ez$t; S <- ez$S

  rc <- rmst_curve(t, S)
  series <- setdiff(names(rc), "tau")

  mu_mat <- rbind(rep(0, length(series)), as.matrix(rc[, series, drop = FALSE]))

  out <- lapply(seq_along(series), function(j) {
    Sj <- S[, j, drop = FALSE]
    muj <- mu_mat[, j, drop = FALSE]

    mu_at <- function(x) as.numeric(.surv_at(t, muj, x))
    S_at  <- function(x) as.numeric(.surv_at(t, Sj,  x))

    muc <- tvrmst_cond(t, Sj, s_grid, tau, eps = eps)[[2]]

    lhs <- mu_at(s_grid + tau) - mu_at(s_grid)
    rhs <- S_at(s_grid) * muc

    data.frame(
      series = series[j],
      max_abs_error = max(abs(lhs - rhs), na.rm = TRUE)
    )
  })

  do.call(rbind, out)
}
