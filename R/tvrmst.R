#' Window RMST difference W(s, tau) = \eqn{\int_s^{s+\tau} (S_1 - S_0)\,du}
#'
#' @param S1 Survival matrix (n_units x n_time) for arm 1.
#' @param S0 Survival matrix (n_units x n_time) for arm 0.
#' @param time Time grid.
#' @param s_grid Nonnegative landmark grid.
#' @param tau Positive window length.
#' @return data.frame(s, <unit...>)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' S1 <- rbind(A = exp(-0.15 * t), B = exp(-0.25 * t))
#' s_grid <- c(0, 1, 2)
#' rmst_window(S1, S0, t, s_grid, tau = 2)
#' @export
rmst_window <- function(S1, S0, time, s_grid, tau) {
  S1 <- .coerce_unit_time(S1, time, "S1")
  S0 <- .coerce_unit_time(S0, time, "S0")
  if (nrow(S1) != nrow(S0)) .stop("S1 and S0 must have same number of rows (paired units).")
  if (!is.numeric(s_grid) || any(!is.finite(s_grid)) || any(s_grid < 0)) {
    .stop("`s_grid` must be nonnegative finite numeric.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    .stop("`tau` must be a single positive number.")
  }

  time <- as.numeric(time)

  out <- lapply(s_grid, function(s) {
    # build explicit endpoints (s, s+tau) even if not in t
    t_win <- sort(unique(c(time[time >= s & time <= (s + tau)], s, s + tau)))
    if (length(t_win) < 2) return(rep(NA_real_, nrow(S1)))

    S1w <- .surv_at(time, S1, t_win)
    S0w <- .surv_at(time, S0, t_win)

    .trapz_cols(t_win, S1w - S0w)
  })

  mat <- do.call(rbind, out)
  colnames(mat) <- if (!is.null(rownames(S1))) rownames(S1) else paste0("unit", seq_len(nrow(S1)))
  data.frame(s = s_grid, mat, check.names = FALSE)
}

#' Conditional tvRMST mu_c(s, tau) = \eqn{(1 / S(s)) \int_s^{s+\tau} S(u)\,du}
#'
#' @param S Survival matrix (n_units x n_time).
#' @param time Time grid.
#' @param s_grid Nonnegative landmark grid.
#' @param tau Positive window length.
#' @param eps Stability threshold: require S(s) >= eps.
#' @return data.frame(s, <unit...>)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' s_grid <- c(0, 1, 2)
#' tvrmst_cond(S, t, s_grid, tau = 2)
#' @export
tvrmst_cond <- function(S, time, s_grid, tau, eps = 1e-8) {
  S <- .coerce_unit_time(S, time, "S")
  if (!is.numeric(s_grid) || any(!is.finite(s_grid)) || any(s_grid < 0)) {
    .stop("`s_grid` must be nonnegative finite numeric.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    .stop("`tau` must be a single positive number.")
  }
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0 || eps >= 1) {
    .stop("`eps` must be in (0,1).")
  }

  time <- as.numeric(time)

  out <- lapply(s_grid, function(s) {
    Ss <- as.numeric(.surv_at(time, S, s))
    ok <- Ss >= eps

    t_win <- sort(unique(c(time[time >= s & time <= (s + tau)], s, s + tau)))
    if (length(t_win) < 2) return(rep(NA_real_, nrow(S)))

    Sw <- .surv_at(time, S, t_win)
    integ <- .trapz_cols(t_win, Sw)

    mu <- rep(NA_real_, nrow(S))
    mu[ok] <- integ[ok] / Ss[ok]
    mu
  })

  mat <- do.call(rbind, out)
  colnames(mat) <- if (!is.null(rownames(S))) rownames(S) else paste0("unit", seq_len(nrow(S)))
  data.frame(s = s_grid, mat, check.names = FALSE)
}

#' Conditional tvRMST difference Delta_c(s, tau) = mu_c1(s, tau) - mu_c0(s, tau)
#'
#' @param time Time grid.
#' @param S1 Survival matrix for arm 1.
#' @param S0 Survival matrix for arm 0.
#' @param s_grid Landmark grid.
#' @param tau Window length.
#' @param eps Stability threshold.
#' @param statistic Aggregation: "mean", "median", or "unit".
#' @return data.frame(s, estimate) for statistic != "unit";
#'   data.frame(s, <unit...>) for statistic = "unit".
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' S1 <- rbind(A = exp(-0.15 * t), B = exp(-0.25 * t))
#' s_grid <- c(0, 1, 2)
#' tvrmst_diff(S1, S0, t, s_grid, tau = 2)
#' @export
tvrmst_diff <- function(S1, S0, time, s_grid, tau, eps = 1e-8,
                        statistic = c("mean", "median", "unit")) {
  statistic <- match.arg(statistic)
  S1 <- .coerce_unit_time(S1, time, "S1")
  S0 <- .coerce_unit_time(S0, time, "S0")

  mu1 <- tvrmst_cond(S1, time, s_grid, tau, eps = eps)
  mu0 <- tvrmst_cond(S0, time, s_grid, tau, eps = eps)

  mu1_mat <- as.matrix(mu1[-1])
  mu0_mat <- as.matrix(mu0[-1])

  if (statistic == "unit") {
    if (nrow(S1) != nrow(S0)) .stop("S1 and S0 must have same number of rows for statistic = \"unit\".")
    delta <- mu1
    delta[-1] <- mu1_mat - mu0_mat
    if (ncol(delta) == 2) names(delta)[2] <- "delta_c"
    return(delta)
  }

  if (statistic == "mean") {
    est <- rowMeans(mu1_mat) - rowMeans(mu0_mat)
  } else {
    est <- apply(mu1_mat, 1, stats::median) - apply(mu0_mat, 1, stats::median)
  }

  data.frame(s = s_grid, estimate = est)
}
