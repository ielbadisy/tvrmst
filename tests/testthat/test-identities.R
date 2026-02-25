library(testthat)

test_that("check_identities returns expected schema for mean and unit", {
  time <- seq(0, 10, by = 1)
  S <- rbind(A = exp(-0.10 * time), B = exp(-0.20 * time))
  s_grid <- seq(0, 8, by = 2)
  tau <- 2

  res_mean <- check_identities(S, time, s_grid, tau, eps = 1e-6, statistic = "mean")
  expect_equal(names(res_mean), c("s", "max_abs_error"))
  expect_equal(nrow(res_mean), length(s_grid))

  res_unit <- check_identities(S, time, s_grid, tau, eps = 1e-6, statistic = "unit")
  expect_true(all(c("s", "unit", "max_abs_error") %in% names(res_unit)))
  expect_equal(nrow(res_unit), nrow(S) * length(s_grid))
})
