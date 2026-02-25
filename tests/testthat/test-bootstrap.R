library(testthat)
test_that("bootstrap_tvrmst_diff_reps collapses for identical reps", {
  t <- seq(0, 6, by = 1)
  s_grid <- seq(0, 4, by = 1)
  tau <- 2

  S0 <- matrix(exp(-0.20 * t), nrow = 1)
  S1 <- matrix(exp(-0.15 * t), nrow = 1)

  reps <- replicate(10, list(S1 = S1, S0 = S0), simplify = FALSE)
  out <- bootstrap_tvrmst_diff_reps(reps, t, s_grid, tau, eps = 1e-6, conf = 0.95)

  expect_equal(names(out), c("s", "estimate", "lower", "upper"))
  expect_true(all(abs(out$lower - out$estimate) < 1e-12))
  expect_true(all(abs(out$upper - out$estimate) < 1e-12))
})
