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
