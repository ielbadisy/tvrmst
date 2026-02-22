#' Window RMST difference W(s, tau) = \eqn{\int_s^{s+\tau} (S_1 - S_0)\,du}
#'
#' @param t Time grid.
#' @param S1 Survival matrix (n_time x n_series) for arm 1.
#' @param S0 Survival matrix (n_time x n_series) for arm 0.
#' @param s_grid Nonnegative landmark grid.
#' @param tau Positive window length.
#' @return data.frame(s, <series...>)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' S1 <- cbind(A = exp(-0.15 * t), B = exp(-0.25 * t))
#' s_grid <- c(0, 1, 2)
#' rmst_window(t, S1, S0, s_grid, tau = 2)
#' @export
rmst_window <- function(t, S1, S0, s_grid, tau) {
  .check_survmat(t, S1, "S1"); .check_survmat(t, S0, "S0")
  if (ncol(S1) != ncol(S0)) .stop("S1 and S0 must have same number of columns (paired series).")
  if (!is.numeric(s_grid) || any(!is.finite(s_grid)) || any(s_grid < 0)) {
    .stop("`s_grid` must be nonnegative finite numeric.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    .stop("`tau` must be a single positive number.")
  }

  ss1 <- .sort_survmat(t, S1); t <- ss1$t; S1 <- ss1$S
  ss0 <- .sort_survmat(t, S0);           S0 <- ss0$S

  out <- lapply(s_grid, function(s) {
    # build explicit endpoints (s, s+tau) even if not in t
    t_win <- sort(unique(c(t[t >= s & t <= (s + tau)], s, s + tau)))
    if (length(t_win) < 2) return(rep(NA_real_, ncol(S1)))

    S1w <- .surv_at(t, S1, t_win)
    S0w <- .surv_at(t, S0, t_win)

    .trapz_cols(t_win, S1w - S0w)
  })

  mat <- do.call(rbind, out)
  colnames(mat) <- if (!is.null(colnames(S1))) colnames(S1) else paste0("series", seq_len(ncol(S1)))
  data.frame(s = s_grid, mat, check.names = FALSE)
}

#' Conditional tvRMST mu_c(s, tau) = \eqn{(1 / S(s)) \int_s^{s+\tau} S(u)\,du}
#'
#' @param t Time grid.
#' @param S Survival matrix (n_time x n_series).
#' @param s_grid Nonnegative landmark grid.
#' @param tau Positive window length.
#' @param eps Stability threshold: require S(s) >= eps.
#' @return data.frame(s, <series...>)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' s_grid <- c(0, 1, 2)
#' tvrmst_cond(t, S, s_grid, tau = 2)
#' @export
tvrmst_cond <- function(t, S, s_grid, tau, eps = 0.05) {
  .check_survmat(t, S, "S")
  if (!is.numeric(s_grid) || any(!is.finite(s_grid)) || any(s_grid < 0)) {
    .stop("`s_grid` must be nonnegative finite numeric.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    .stop("`tau` must be a single positive number.")
  }
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0 || eps >= 1) {
    .stop("`eps` must be in (0,1).")
  }

  ss <- .sort_survmat(t, S); t <- ss$t; S <- ss$S

  out <- lapply(s_grid, function(s) {
    Ss <- .surv_at(t, S, s)              # 1 x ncol(S)
    Ss <- as.numeric(Ss)
    ok <- Ss >= eps

    t_win <- sort(unique(c(t[t >= s & t <= (s + tau)], s, s + tau)))
    if (length(t_win) < 2) return(rep(NA_real_, ncol(S)))

    Sw <- .surv_at(t, S, t_win)
    integ <- .trapz_cols(t_win, Sw)

    mu <- rep(NA_real_, ncol(S))
    mu[ok] <- integ[ok] / Ss[ok]
    mu
  })

  mat <- do.call(rbind, out)
  colnames(mat) <- if (!is.null(colnames(S))) colnames(S) else paste0("series", seq_len(ncol(S)))
  data.frame(s = s_grid, mat, check.names = FALSE)
}

#' Conditional tvRMST difference Delta_c(s, tau) = mu_c1(s, tau) - mu_c0(s, tau)
#'
#' @param t Time grid.
#' @param S1 Survival matrix for arm 1.
#' @param S0 Survival matrix for arm 0.
#' @param s_grid Landmark grid.
#' @param tau Window length.
#' @param eps Stability threshold.
#' @return data.frame(s, delta_c) if single series, else data.frame(s, <series...>)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' S1 <- cbind(A = exp(-0.15 * t), B = exp(-0.25 * t))
#' s_grid <- c(0, 1, 2)
#' tvrmst_diff(t, S1, S0, s_grid, tau = 2)
#' @export
tvrmst_diff <- function(t, S1, S0, s_grid, tau, eps = 0.05) {
  .check_survmat(t, S1, "S1"); .check_survmat(t, S0, "S0")
  if (ncol(S1) != ncol(S0)) .stop("S1 and S0 must have same number of columns (paired series).")

  mu1 <- tvrmst_cond(t, S1, s_grid, tau, eps = eps)
  mu0 <- tvrmst_cond(t, S0, s_grid, tau, eps = eps)

  delta <- mu1
  delta[-1] <- mu1[-1] - mu0[-1]

  if (ncol(delta) == 2) names(delta)[2] <- "delta_c"
  delta
}
