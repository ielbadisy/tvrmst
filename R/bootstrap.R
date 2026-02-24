#' Bootstrap Delta_c(s, tau) by resampling rows (units) within each arm
#'
#' @param time Time grid
#' @param S1_units n_units x n_time survival curves for treatment units
#' @param S0_units n_units x n_time survival curves for control units
#' @param s_grid Landmark grid
#' @param tau Window length
#' @param R Bootstrap replicates
#' @param eps Stability threshold
#' @param conf Confidence level
#' @param statistic Aggregator across units: "mean" or "median"
#' @return data.frame(s, estimate, lower, upper)
#' @examples
#' set.seed(1)
#' t <- seq(0, 5, by = 1)
#' s_grid <- c(0, 1, 2)
#' S0_units <- t(sapply(1:6, function(i) exp(-(0.20 + runif(1, -0.02, 0.02)) * t)))
#' S1_units <- t(sapply(1:6, function(i) exp(-(0.15 + runif(1, -0.02, 0.02)) * t)))
#' bootstrap_tvrmst_diff_cols(S1_units, S0_units, t, s_grid, tau = 2, R = 20, conf = 0.9)
#' @export
bootstrap_tvrmst_diff_cols <- function(S1_units, S0_units, time, s_grid, tau,
                                       R = 500, eps = 1e-8, conf = 0.95,
                                       statistic = c("mean","median")) {
  statistic <- match.arg(statistic)
  S1_units <- .coerce_unit_time(S1_units, time, "S1_units")
  S0_units <- .coerce_unit_time(S0_units, time, "S0_units")
  n1 <- nrow(S1_units); n0 <- nrow(S0_units)
  if (n1 < 5 || n0 < 5) .stop("Need at least 5 unit curves per arm for bootstrap.")
  if (!is.numeric(R) || length(R) != 1 || R < 10) .stop("`R` must be >= 10.")
  if (!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) .stop("`conf` must be in (0,1).")

  est_df <- tvrmst_diff(S1_units, S0_units, time, s_grid, tau, eps = eps, statistic = statistic)
  est_vec <- est_df$estimate

  boot_mat <- matrix(NA_real_, nrow = R, ncol = length(s_grid))
  for (b in seq_len(R)) {
    id1 <- sample.int(n1, size = n1, replace = TRUE)
    id0 <- sample.int(n0, size = n0, replace = TRUE)
    S1b <- S1_units[id1, , drop = FALSE]
    S0b <- S0_units[id0, , drop = FALSE]
    boot_mat[b, ] <- as.numeric(
      tvrmst_diff(S1b, S0b, time, s_grid, tau, eps = eps, statistic = statistic)$estimate
    )
  }

  alpha <- 1 - conf
  lower <- apply(boot_mat, 2, stats::quantile, probs = alpha/2, na.rm = TRUE, names = FALSE)
  upper <- apply(boot_mat, 2, stats::quantile, probs = 1-alpha/2, na.rm = TRUE, names = FALSE)

  data.frame(s = s_grid, estimate = est_vec, lower = lower, upper = upper)
}

#' Summarize Delta_c(s, tau) from a list of replicate survival matrices
#'
#' @param time Time grid
#' @param reps List of length R, each element list(S1=..., S0=...)
#' @param s_grid Landmark grid
#' @param tau Window length
#' @param eps Stability threshold
#' @param conf Confidence level
#' @return data.frame(s, estimate, lower, upper)
#' @examples
#' set.seed(1)
#' t <- seq(0, 5, by = 1)
#' s_grid <- c(0, 1, 2)
#' reps <- lapply(1:10, function(i) {
#'   list(
#'     S1 = matrix(exp(-(0.15 + runif(1, -0.02, 0.02)) * t), nrow = 1),
#'     S0 = matrix(exp(-(0.20 + runif(1, -0.02, 0.02)) * t), nrow = 1)
#'   )
#' })
#' bootstrap_tvrmst_diff_reps(reps, t, s_grid, tau = 2, conf = 0.9)
#' @export
bootstrap_tvrmst_diff_reps <- function(reps, time, s_grid, tau, eps = 1e-8, conf = 0.95) {
  if (!is.list(reps) || length(reps) < 10) .stop("`reps` must be a list of >= 10 replicates.")
  if (is.null(reps[[1]]$S1) || is.null(reps[[1]]$S0)) .stop("Each replicate must be list(S1=..., S0=...).")

  S1_first <- .coerce_unit_time(reps[[1]]$S1, time, "reps[[1]]$S1")
  S0_first <- .coerce_unit_time(reps[[1]]$S0, time, "reps[[1]]$S0")
  if (nrow(S1_first) != nrow(S0_first)) .stop("Replicate 1: S1 and S0 must have same number of rows.")
  if (nrow(S1_first) != 1) .stop("Replicates must contain a single curve per arm.")

  R <- length(reps)
  boot_mat <- matrix(NA_real_, nrow = R, ncol = length(s_grid))
  for (r in seq_len(R)) {
    S1 <- .coerce_unit_time(reps[[r]]$S1, time, sprintf("reps[[%d]]$S1", r))
    S0 <- .coerce_unit_time(reps[[r]]$S0, time, sprintf("reps[[%d]]$S0", r))
    if (nrow(S1) != nrow(S0)) .stop("Replicate %d: S1 and S0 must have same number of rows.", r)
    if (nrow(S1) != 1) .stop("Replicate %d: each arm must contain a single curve.", r)
    boot_mat[r, ] <- as.numeric(
      tvrmst_diff(S1, S0, time, s_grid, tau, eps = eps, statistic = "unit")[[2]]
    )
  }

  S1_mean <- Reduce(`+`, lapply(reps, function(rp) .coerce_unit_time(rp$S1, time, "rep$S1"))) / R
  S0_mean <- Reduce(`+`, lapply(reps, function(rp) .coerce_unit_time(rp$S0, time, "rep$S0"))) / R
  est <- as.numeric(
    tvrmst_diff(S1_mean, S0_mean, time, s_grid, tau, eps = eps, statistic = "unit")[[2]]
  )

  alpha <- 1 - conf
  lower <- apply(boot_mat, 2, stats::quantile, probs = alpha/2, na.rm = TRUE, names = FALSE)
  upper <- apply(boot_mat, 2, stats::quantile, probs = 1-alpha/2, na.rm = TRUE, names = FALSE)

  data.frame(s = s_grid, estimate = est, lower = lower, upper = upper)
}
