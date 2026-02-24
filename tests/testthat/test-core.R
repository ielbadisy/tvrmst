library(testthat)


test_that("rmst_tau handles basic and edge cases", {
  t <- seq(0, 10, by = 1)
  S <- matrix(1, nrow = 2, ncol = length(t))

  out <- rmst_tau(S, t, tau = 5)
  expect_equal(unname(out), c(5, 5))

  # tau beyond max(t) should extend with constant survival
  out2 <- rmst_tau(S, t, tau = 12)
  expect_equal(unname(out2), c(12, 12))

  # tau before first positive time when t does not include 0
  t2 <- seq(1, 5, by = 1)
  S2 <- matrix(1, nrow = 1, ncol = length(t2))
  out3 <- rmst_tau(S2, t2, tau = 0.5)
  expect_equal(unname(out3), 0.5)

  expect_error(rmst_tau(S, t, tau = -1))
})


test_that("rmst_curve and rmst_delta_curve behave as expected", {
  t <- seq(0, 4, by = 1)
  S <- rbind(rep(1, length(t)), rep(1, length(t)))

  rc <- rmst_curve(S, t)
  expect_equal(rc$tau, t[-1])
  expect_equal(rc[[1 + 1]], t[-1])
  expect_equal(rc[[2 + 1]], t[-1])

  delta <- rmst_delta_curve(rc, arm1 = colnames(rc)[2], arm0 = colnames(rc)[3])
  expect_true(all(delta$delta == 0))
})


test_that("rmst_window integrates windowed differences", {
  t <- seq(0, 10, by = 1)
  S1 <- matrix(1, nrow = 1, ncol = length(t))
  S0 <- matrix(0.5, nrow = 1, ncol = length(t))
  s_grid <- seq(0, 5, by = 1)
  tau <- 2

  w <- rmst_window(S1, S0, t, s_grid, tau)
  expect_equal(nrow(w), length(s_grid))
  expect_equal(w[[2]], rep(1, length(s_grid)))
})


test_that("tvrmst_cond returns tau for constant survival", {
  t <- seq(0, 10, by = 1)
  S <- matrix(1, nrow = 1, ncol = length(t))
  s_grid <- seq(0, 8, by = 1)
  tau <- 3

  mu <- tvrmst_cond(S, t, s_grid, tau, eps = 1e-6)
  expect_equal(mu[[2]], rep(tau, length(s_grid)))
})


test_that("tvrmst_diff is zero when arms are identical", {
  t <- seq(0, 8, by = 1)
  S <- matrix(exp(-0.2 * t), nrow = 1)
  s_grid <- seq(0, 6, by = 1)
  tau <- 2

  delta <- tvrmst_diff(S, S, t, s_grid, tau, eps = 1e-6)
  expect_true(all(abs(delta[[2]]) < 1e-12))
})


test_that("plot helpers return ggplot objects when ggplot2 is available", {
  skip_if_not_installed("ggplot2")

  t <- seq(0, 4, by = 1)
  S0 <- matrix(exp(-0.2 * t), nrow = 1)
  S1 <- matrix(exp(-0.1 * t), nrow = 1)

  p1 <- plot_survival_curves(S0, S1, t)
  expect_s3_class(p1, "ggplot")

  s_grid <- seq(0, 3, by = 1)
  mu0 <- tvrmst_cond(S0, t, s_grid, tau = 1)
  mu1 <- tvrmst_cond(S1, t, s_grid, tau = 1)

  p2 <- plot_tvrmst(mu0, mu1)
  expect_s3_class(p2, "ggplot")

  delta <- tvrmst_diff(S1, S0, t, s_grid, tau = 1)
  p3 <- plot_tvrmst_diff(delta)
  expect_s3_class(p3, "ggplot")
})
