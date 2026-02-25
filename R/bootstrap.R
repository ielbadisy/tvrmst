#' Bootstrap summary for tvRMST differences
#'
#' @param reps List of bootstrap replicates. Each element must be
#'   `list(S1 = <matrix/data.frame>, S0 = <matrix/data.frame>)`.
#' @param time Numeric strictly increasing time vector.
#' @param s_grid Landmark times.
#' @param tau Positive window length.
#' @param eps Positivity threshold passed to `tvrmst_diff()`.
#' @param conf Confidence level in (0,1).
#' @param statistic Summary statistic passed to `tvrmst_diff()`.
#' @return data.frame with columns `s`, `estimate`, `lower`, `upper`.
#' @export
bootstrap_tvrmst_diff_reps <- function(reps, time, s_grid, tau,
                                       eps = 1e-8, conf = 0.95,
                                       statistic = "mean") {
  if (!is.list(reps) || length(reps) < 2) {
    .stop("`reps` must be a list with at least 2 replicates.")
  }
  if (!is.numeric(conf) || length(conf) != 1 || !is.finite(conf) || conf <= 0 || conf >= 1) {
    .stop("`conf` must be a single number in (0,1).")
  }

  statistic <- match.arg(statistic, c("mean", "median"))

  R <- length(reps)
  boot_mat <- matrix(NA_real_, nrow = R, ncol = length(s_grid))

  for (r in seq_len(R)) {
    rep_r <- reps[[r]]
    if (!is.list(rep_r) || is.null(rep_r$S1) || is.null(rep_r$S0)) {
      .stop("Each replicate must be list(S1 = <matrix/data.frame>, S0 = <matrix/data.frame>).")
    }

    validate_surv_input(rep_r$S1, time, name = sprintf("reps[[%d]]$S1", r))
    validate_surv_input(rep_r$S0, time, name = sprintf("reps[[%d]]$S0", r))

    dr <- tvrmst_diff(
      S1 = rep_r$S1,
      S0 = rep_r$S0,
      time = time,
      s_grid = s_grid,
      tau = tau,
      eps = eps,
      statistic = statistic
    )

    if (!is.data.frame(dr) || !all(c("s", "estimate") %in% names(dr))) {
      .stop("Each bootstrap replicate must produce a data.frame with columns `s` and `estimate`.")
    }

    boot_mat[r, ] <- as.numeric(dr$estimate)
  }

  alpha <- 1 - conf
  data.frame(
    s = s_grid,
    estimate = colMeans(boot_mat, na.rm = TRUE),
    lower = apply(boot_mat, 2, stats::quantile, probs = alpha / 2, na.rm = TRUE, names = FALSE),
    upper = apply(boot_mat, 2, stats::quantile, probs = 1 - alpha / 2, na.rm = TRUE, names = FALSE)
  )
}
