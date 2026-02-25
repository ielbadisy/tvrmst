library(testthat)

test_that("bootstrap_tvrmst_diff_reps returns required columns and ribbon-ready output", {
  time <- seq(0, 6, by = 1)
  s_grid <- seq(0, 4, by = 1)
  tau <- 2

  S0 <- matrix(exp(-0.20 * time), nrow = 1)
  S1 <- matrix(exp(-0.15 * time), nrow = 1)

  reps <- replicate(12, list(S1 = S1, S0 = S0), simplify = FALSE)
  out <- bootstrap_tvrmst_diff_reps(reps, time, s_grid, tau, eps = 1e-6, conf = 0.95)

  expect_equal(names(out), c("s", "estimate", "lower", "upper"))
  expect_equal(nrow(out), length(s_grid))

  skip_if_not_installed("ggplot2")
  p <- plot_tvrmst_diff(out)
  expect_s3_class(p, "ggplot")
})
