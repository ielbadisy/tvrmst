library(testthat)
test_that("check_identities holds numerically", {
  t <- seq(0, 10, by = 1)
  S <- rbind(exp(-0.10 * t), exp(-0.20 * t))
  s_grid <- seq(0, 8, by = 1)
  tau <- 2

  res <- check_identities(S, t, s_grid, tau, eps = 1e-6)
  expect_equal(nrow(res), 2)
  expect_true(all(res$max_abs_error < 1e-8))
})
