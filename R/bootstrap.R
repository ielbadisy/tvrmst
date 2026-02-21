#' Bootstrap Delta_c(s, tau) by resampling columns (units) within each arm
#'
#' @param t Time grid
#' @param S1_units n_time x n_units survival curves for treatment units
#' @param S0_units n_time x n_units survival curves for control units
#' @param s_grid Landmark grid
#' @param tau Window length
#' @param R Bootstrap replicates
#' @param eps Stability threshold
#' @param conf Confidence level
#' @param statistic Aggregator across units: "mean" or "median"
#' @return data.frame(s, estimate, lower, upper)
#' @export
bootstrap_tvrmst_diff_cols <- function(t, S1_units, S0_units, s_grid, tau,
                                       R = 500, eps = 0.05, conf = 0.95,
                                       statistic = c("mean","median")) {
  statistic <- match.arg(statistic)
  .check_survmat(t, S1_units, "S1_units")
  .check_survmat(t, S0_units, "S0_units")
  n1 <- ncol(S1_units); n0 <- ncol(S0_units)
  if (n1 < 5 || n0 < 5) .stop("Need at least 5 unit curves per arm for bootstrap.")
  if (!is.numeric(R) || length(R) != 1 || R < 10) .stop("`R` must be >= 10.")
  if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) .stop("`conf` must be in (0,1).")

  agg <- function(M) {
    if (statistic == "mean") rowMeans(M) else apply(M, 1, stats::median)
  }

  S1_hat <- matrix(agg(S1_units), ncol = 1)
  S0_hat <- matrix(agg(S0_units), ncol = 1)
  est_df <- tvrmst_diff(t, S1_hat, S0_hat, s_grid, tau, eps = eps)
  est_vec <- est_df[[2]]

  boot_mat <- matrix(NA_real_, nrow = R, ncol = length(s_grid))
  for (b in seq_len(R)) {
    id1 <- sample.int(n1, size = n1, replace = TRUE)
    id0 <- sample.int(n0, size = n0, replace = TRUE)
    S1b <- matrix(agg(S1_units[, id1, drop = FALSE]), ncol = 1)
    S0b <- matrix(agg(S0_units[, id0, drop = FALSE]), ncol = 1)
    boot_mat[b, ] <- tvrmst_diff(t, S1b, S0b, s_grid, tau, eps = eps)[[2]]
  }

  alpha <- 1 - conf
  lower <- apply(boot_mat, 2, stats::quantile, probs = alpha/2, na.rm = TRUE, names = FALSE)
  upper <- apply(boot_mat, 2, stats::quantile, probs = 1-alpha/2, na.rm = TRUE, names = FALSE)

  data.frame(s = s_grid, estimate = est_vec, lower = lower, upper = upper)
}

#' Summarize Delta_c(s, tau) from a list of replicate survival matrices
#'
#' @param t Time grid
#' @param reps List of length R, each element list(S1=..., S0=...)
#' @param s_grid Landmark grid
#' @param tau Window length
#' @param eps Stability threshold
#' @param conf Confidence level
#' @return data.frame(s, estimate, lower, upper)
#' @export
bootstrap_tvrmst_diff_reps <- function(t, reps, s_grid, tau, eps = 0.05, conf = 0.95) {
  if (!is.list(reps) || length(reps) < 10) .stop("`reps` must be a list of >= 10 replicates.")
  if (is.null(reps[[1]]$S1) || is.null(reps[[1]]$S0)) .stop("Each replicate must be list(S1=..., S0=...).")

  .check_survmat(t, reps[[1]]$S1, "reps[[1]]$S1")
  .check_survmat(t, reps[[1]]$S0, "reps[[1]]$S0")

  R <- length(reps)
  boot_mat <- matrix(NA_real_, nrow = R, ncol = length(s_grid))
  for (r in seq_len(R)) {
    S1 <- reps[[r]]$S1
    S0 <- reps[[r]]$S0
    if (ncol(S1) != ncol(S0)) .stop("Replicate %d: S1 and S0 must have same number of columns.", r)
    boot_mat[r, ] <- tvrmst_diff(t, S1, S0, s_grid, tau, eps = eps)[[2]]
  }

  S1_mean <- Reduce(`+`, lapply(reps, `[[`, "S1")) / R
  S0_mean <- Reduce(`+`, lapply(reps, `[[`, "S0")) / R
  est <- tvrmst_diff(t, S1_mean, S0_mean, s_grid, tau, eps = eps)[[2]]

  alpha <- 1 - conf
  lower <- apply(boot_mat, 2, stats::quantile, probs = alpha/2, na.rm = TRUE, names = FALSE)
  upper <- apply(boot_mat, 2, stats::quantile, probs = 1-alpha/2, na.rm = TRUE, names = FALSE)

  data.frame(s = s_grid, estimate = est, lower = lower, upper = upper)
}
