.validate_landmark <- function(s_grid, tau) {
  if (!is.numeric(s_grid) || length(s_grid) < 1 || any(!is.finite(s_grid)) || any(s_grid < 0)) {
    .stop("`s_grid` must be a finite numeric vector with nonnegative values.")
  }
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau <= 0) {
    .stop("`tau` must be a single positive number.")
  }
}

.window_integral_matrix <- function(S, time, s_grid, tau) {
  n_units <- nrow(S)
  out <- matrix(NA_real_, nrow = n_units, ncol = length(s_grid))

  for (j in seq_along(s_grid)) {
    s <- s_grid[j]
    t_win <- sort(unique(c(s, s + tau, time[time >= s & time <= s + tau])))
    if (length(t_win) < 2) next
    Sw <- .surv_at(time, S, t_win)
    out[, j] <- .trapz_cols(t_win, Sw)
  }

  if (!is.null(rownames(S))) rownames(out) <- rownames(S)
  colnames(out) <- as.character(s_grid)
  out
}

#' Window RMST difference
#'
#' @param S1 Arm-1 survival matrix/data.frame (rows = units, cols = time).
#' @param S0 Arm-0 survival matrix/data.frame (rows = units, cols = time).
#' @param time Numeric strictly increasing time vector.
#' @param s_grid Landmark times.
#' @param tau Positive window length.
#' @param statistic One of `"mean"`, `"median"`, `"unit"`.
#' @return For summary statistics: data.frame `s`, `estimate`.
#'   For `statistic = "unit"`: numeric matrix n_units x length(s_grid).
#' @export
rmst_window <- function(S1, S0, time, s_grid, tau,
                        statistic = c("mean", "median", "unit")) {
  S1 <- validate_surv_input(S1, time, name = "S1")
  S0 <- validate_surv_input(S0, time, name = "S0")
  statistic <- match.arg(statistic)
  time <- as.numeric(time)
  .validate_landmark(s_grid, tau)

  int1 <- .window_integral_matrix(S1, time, s_grid, tau)
  int0 <- .window_integral_matrix(S0, time, s_grid, tau)

  if (statistic == "unit") {
    if (nrow(S1) != nrow(S0)) {
      .stop("S1 and S0 must have same number of rows for statistic = \"unit\".")
    }
    out <- int1 - int0
    if (!is.null(rownames(S1))) rownames(out) <- rownames(S1)
    return(out)
  }

  if (statistic == "mean") {
    estimate <- colMeans(int1, na.rm = TRUE) - colMeans(int0, na.rm = TRUE)
  } else {
    estimate <- apply(int1, 2, stats::median, na.rm = TRUE) -
      apply(int0, 2, stats::median, na.rm = TRUE)
  }

  data.frame(s = s_grid, estimate = as.numeric(estimate))
}

#' Conditional tvRMST
#'
#' @param S Survival matrix/data.frame (rows = units, cols = time).
#' @param time Numeric strictly increasing time vector.
#' @param s_grid Landmark times.
#' @param tau Positive window length.
#' @param eps Positivity threshold for `S(s)`.
#' @param statistic One of `"mean"`, `"median"`, `"unit"`.
#' @return For summary statistics: data.frame `s`, `estimate`.
#'   For `statistic = "unit"`: numeric matrix n_units x length(s_grid).
#' @export
tvrmst_cond <- function(S, time, s_grid, tau, eps = 1e-8,
                        statistic = c("mean", "median", "unit")) {
  S <- validate_surv_input(S, time, name = "S")
  statistic <- match.arg(statistic)
  time <- as.numeric(time)
  .validate_landmark(s_grid, tau)

  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0 || eps >= 1) {
    .stop("`eps` must be a single number in (0,1).")
  }

  num <- .window_integral_matrix(S, time, s_grid, tau)
  denom <- .surv_at(time, S, s_grid)

  mu <- matrix(NA_real_, nrow = nrow(S), ncol = length(s_grid))
  ok <- denom >= eps
  mu[ok] <- num[ok] / denom[ok]

  if (!is.null(rownames(S))) rownames(mu) <- rownames(S)
  colnames(mu) <- as.character(s_grid)

  if (statistic == "unit") {
    return(mu)
  }

  if (statistic == "mean") {
    estimate <- colMeans(mu, na.rm = TRUE)
  } else {
    estimate <- apply(mu, 2, stats::median, na.rm = TRUE)
  }

  data.frame(s = s_grid, estimate = as.numeric(estimate))
}

#' Difference in conditional tvRMST between arms
#'
#' @param S1 Arm-1 survival matrix/data.frame (rows = units, cols = time).
#' @param S0 Arm-0 survival matrix/data.frame (rows = units, cols = time).
#' @param time Numeric strictly increasing time vector.
#' @param s_grid Landmark times.
#' @param tau Positive window length.
#' @param eps Positivity threshold.
#' @param statistic One of `"mean"`, `"median"`, `"unit"`.
#' @return For summary statistics: data.frame `s`, `estimate`.
#'   For `statistic = "unit"`: numeric matrix n_units x length(s_grid).
#' @export
tvrmst_diff <- function(S1, S0, time, s_grid, tau, eps = 1e-8,
                        statistic = c("mean", "median", "unit")) {
  S1 <- validate_surv_input(S1, time, name = "S1")
  S0 <- validate_surv_input(S0, time, name = "S0")
  statistic <- match.arg(statistic)
  time <- as.numeric(time)

  mu1 <- tvrmst_cond(S1, time, s_grid, tau, eps = eps, statistic = "unit")
  mu0 <- tvrmst_cond(S0, time, s_grid, tau, eps = eps, statistic = "unit")

  if (statistic == "unit") {
    if (nrow(mu1) != nrow(mu0)) {
      .stop("S1 and S0 must have same number of rows for statistic = \"unit\".")
    }
    out <- mu1 - mu0
    if (!is.null(rownames(S1))) rownames(out) <- rownames(S1)
    return(out)
  }

  if (statistic == "mean") {
    estimate <- colMeans(mu1, na.rm = TRUE) - colMeans(mu0, na.rm = TRUE)
  } else {
    estimate <- apply(mu1, 2, stats::median, na.rm = TRUE) -
      apply(mu0, 2, stats::median, na.rm = TRUE)
  }

  data.frame(s = s_grid, estimate = as.numeric(estimate))
}
