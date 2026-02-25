library(testthat)

test_that("1 x m inputs work across rmst/tvrmst/window", {
  time <- seq(0, 8, by = 1)
  s_grid <- seq(0, 6, by = 2)
  tau <- 2

  S0 <- exp(-0.2 * time)
  S1 <- exp(-0.1 * time)

  S0_mat <- rbind(Control = S0)
  S1_mat <- rbind(Treatment = S1)
  S_both <- rbind(Control = S0, Treatment = S1)

  rd <- rmst_dynamic(S_both, time)
  rt <- rmst_tau(S_both, time, tau = 3.2)
  rw <- rmst_window(S1_mat, S0_mat, time, s_grid, tau, statistic = "mean")
  mu0 <- tvrmst_cond(S0_mat, time, s_grid, tau, statistic = "mean")
  mu1 <- tvrmst_cond(S1_mat, time, s_grid, tau, statistic = "mean")
  dc <- tvrmst_diff(S1_mat, S0_mat, time, s_grid, tau, statistic = "mean")

  expect_true(is.matrix(rd))
  expect_length(rt, 2)
  expect_equal(names(rw), c("s", "estimate"))
  expect_equal(names(mu0), c("s", "estimate"))
  expect_equal(names(mu1), c("s", "estimate"))
  expect_equal(names(dc), c("s", "estimate"))
})

test_that("validator catches ncol mismatch", {
  time <- seq(0, 5, by = 1)
  S <- matrix(runif(10), nrow = 2, ncol = 5)
  expect_error(rmst_dynamic(S, time), "ncol\\(S\\) == length\\(time\\)")
})
