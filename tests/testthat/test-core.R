library(testthat)

test_that("rmst_dynamic returns n x m matrix and accepts data.frame", {
  time <- seq(0, 6, by = 1)
  S <- rbind(A = exp(-0.2 * time), B = exp(-0.1 * time))

  out <- rmst_dynamic(S, time)
  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(S))
  expect_equal(colnames(out), as.character(time))

  pred_df <- as.data.frame(S)
  out_df <- rmst_dynamic(pred_df, time)
  expect_true(is.matrix(out_df))
  expect_equal(dim(out_df), dim(S))
})

test_that("rmst_tau snaps tau to nearest grid and stores tau_used", {
  time <- seq(0, 10, by = 2)
  S <- rbind(A = exp(-0.2 * time), B = exp(-0.1 * time))

  out <- rmst_tau(S, time, tau = 5.1)
  expect_equal(attr(out, "tau_used"), 6)
  expect_length(out, nrow(S))
})

test_that("rmst_curve and rmst_delta_curve new schema", {
  time <- seq(0, 4, by = 1)
  S0 <- matrix(exp(-0.2 * time), nrow = 1)
  S1 <- matrix(exp(-0.1 * time), nrow = 1)

  rc0 <- rmst_curve(S0, time, statistic = "mean", probs = NULL)
  rc1 <- rmst_curve(S1, time, statistic = "mean", probs = NULL)

  wide <- data.frame(tau = time, Control = rc0$estimate, Treatment = rc1$estimate)
  dr <- rmst_delta_curve(wide, arm1 = "Treatment", arm0 = "Control")

  expect_equal(names(rc0), c("tau", "estimate"))
  expect_equal(names(dr), c("tau", "estimate"))
  expect_true(all(dr$estimate >= 0))
})

test_that("rmst_window/tvrmst_cond/tvrmst_diff support summary and unit modes", {
  time <- seq(0, 8, by = 1)
  s_grid <- seq(0, 6, by = 2)
  tau <- 2

  S0 <- matrix(exp(-0.2 * time), nrow = 2, ncol = length(time), byrow = TRUE)
  S1 <- matrix(exp(-0.15 * time), nrow = 3, ncol = length(time), byrow = TRUE)

  w_mean <- rmst_window(S1, S0, time, s_grid, tau, statistic = "mean")
  expect_equal(names(w_mean), c("s", "estimate"))

  mu_med <- tvrmst_cond(S1, time, s_grid, tau, statistic = "median")
  expect_equal(names(mu_med), c("s", "estimate"))

  d_mean <- tvrmst_diff(S1, S0, time, s_grid, tau, statistic = "mean")
  expect_equal(names(d_mean), c("s", "estimate"))

  S0_pair <- S0
  S1_pair <- matrix(exp(-0.15 * time), nrow = 2, ncol = length(time), byrow = TRUE)
  w_unit <- rmst_window(S1_pair, S0_pair, time, s_grid, tau, statistic = "unit")
  d_unit <- tvrmst_diff(S1_pair, S0_pair, time, s_grid, tau, statistic = "unit")
  mu_unit <- tvrmst_cond(S1_pair, time, s_grid, tau, statistic = "unit")

  expect_true(is.matrix(w_unit))
  expect_true(is.matrix(d_unit))
  expect_true(is.matrix(mu_unit))
  expect_equal(ncol(w_unit), length(s_grid))
})

test_that("vector inputs are rejected for compute functions", {
  time <- seq(0, 5, by = 1)
  S0_vec <- exp(-0.2 * time)
  S1_vec <- exp(-0.1 * time)

  expect_error(
    tvrmst_diff(S1_vec, S0_vec, time, s_grid = c(0, 1), tau = 1),
    "must be a numeric matrix or data.frame"
  )
})

test_that("plot helpers consume precomputed objects", {
  skip_if_not_installed("ggplot2")

  time <- seq(0, 5, by = 1)
  S0 <- matrix(exp(-0.2 * time), nrow = 1)
  S1 <- matrix(exp(-0.1 * time), nrow = 1)

  p_surv <- plot_survival_curves(S0, S1, time, show = "mean")
  expect_s3_class(p_surv, "ggplot")

  rc <- rmst_curve(rbind(S0, S1), time, statistic = "mean", probs = NULL)
  p_rc <- plot_rmst_curve(rc)
  expect_s3_class(p_rc, "ggplot")

  dr <- data.frame(tau = rc$tau, estimate = rc$estimate - rc$estimate)
  p_dr <- plot_rmst_delta(dr)
  expect_s3_class(p_dr, "ggplot")

  s_grid <- seq(0, 4, by = 1)
  mu0 <- tvrmst_cond(S0, time, s_grid, tau = 1, statistic = "mean")
  mu1 <- tvrmst_cond(S1, time, s_grid, tau = 1, statistic = "mean")
  delta_c <- tvrmst_diff(S1, S0, time, s_grid, tau = 1, statistic = "mean")

  p_mu <- plot_tvrmst(mu0, mu1)
  p_delta <- plot_tvrmst_diff(delta_c)
  expect_s3_class(p_mu, "ggplot")
  expect_s3_class(p_delta, "ggplot")

  RMST_mat <- rmst_dynamic(rbind(S0, S1), time)
  p_ind <- plot_rmst_individual(RMST_mat, group = c("Control", "Treatment"))
  p_mean <- plot_rmst_mean(RMST_mat, group = c("Control", "Treatment"), statistic = "mean")
  expect_s3_class(p_ind, "ggplot")
  expect_s3_class(p_mean, "ggplot")
})
