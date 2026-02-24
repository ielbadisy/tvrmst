library(testthat)
test_that("bootstrap_tvrmst_diff_cols returns finite CI", {
  t <- seq(0, 6, by = 1)
  s_grid <- seq(0, 4, by = 1)
  tau <- 2

  S0_units <- rbind(
    exp(-0.20 * t), exp(-0.22 * t), exp(-0.18 * t), exp(-0.21 * t), exp(-0.19 * t)
  )
  S1_units <- rbind(
    exp(-0.15 * t), exp(-0.16 * t), exp(-0.14 * t), exp(-0.155 * t), exp(-0.145 * t)
  )

  set.seed(1)
  out <- bootstrap_tvrmst_diff_cols(S1_units, S0_units, t, s_grid, tau, R = 30, eps = 1e-4)

  expect_equal(names(out), c("s", "estimate", "lower", "upper"))
  expect_equal(nrow(out), length(s_grid))
  expect_true(all(is.finite(out$estimate)))
  expect_true(all(is.finite(out$lower)))
  expect_true(all(is.finite(out$upper)))
  expect_true(all(out$lower <= out$upper))
})


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
