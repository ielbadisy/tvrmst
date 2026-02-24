library(testthat)

test_that("matrix input works with unit-first orientation", {
  time <- seq(0, 10, length.out = 20)
  S <- matrix(runif(50 * length(time)), nrow = 50, ncol = length(time))
  expect_silent(rmst_dynamic(S, time))
})

test_that("data.frame input works with unit-first orientation", {
  time <- seq(0, 10, length.out = 20)
  S <- matrix(runif(10 * length(time)), nrow = 10, ncol = length(time))
  S_df <- as.data.frame(S)
  out <- rmst_dynamic(S_df, time)
  expect_true(is.data.frame(out))
  expect_true("tau" %in% names(out))
})

test_that("error on mismatched ncol and time length", {
  time <- seq(0, 10, length.out = 20)
  S <- matrix(runif(50 * (length(time) - 1)), nrow = 50, ncol = length(time) - 1)
  expect_error(rmst_dynamic(S, time), "ncol\\(S\\) == length\\(time\\)")
})

test_that("two-arm mean diff works with unequal rows", {
  time <- seq(0, 6, by = 1)
  s_grid <- seq(0, 4, by = 1)
  tau <- 2
  S1 <- matrix(exp(-0.15 * time), nrow = 3, ncol = length(time), byrow = TRUE)
  S0 <- matrix(exp(-0.20 * time), nrow = 5, ncol = length(time), byrow = TRUE)

  out <- tvrmst_diff(S1, S0, time, s_grid, tau, statistic = "mean")
  expect_equal(names(out), c("s", "estimate"))
  expect_equal(nrow(out), length(s_grid))
})

test_that("two-arm unit diff requires equal rows", {
  time <- seq(0, 6, by = 1)
  s_grid <- seq(0, 4, by = 1)
  tau <- 2
  S1 <- matrix(exp(-0.15 * time), nrow = 3, ncol = length(time), byrow = TRUE)
  S0 <- matrix(exp(-0.20 * time), nrow = 5, ncol = length(time), byrow = TRUE)

  expect_error(
    tvrmst_diff(S1, S0, time, s_grid, tau, statistic = "unit"),
    "same number of rows"
  )
})
