library(testthat)


test_that("as_survmat handles canonical matrix", {
  t <- seq(0, 9, by = 1)
  S <- matrix(runif(50), nrow = length(t), ncol = 5)

  sm <- as_survmat(S, t = t)
  expect_s3_class(sm, "survmat")
  expect_identical(sm$t, t)
  expect_equal(nrow(sm$S), length(t))
})


test_that("as_survmat accepts matrix with times attribute", {
  t <- seq(0, 4, by = 1)
  S <- matrix(runif(15), nrow = 3, ncol = length(t))
  attr(S, "times") <- t

  sm <- as_survmat(S)
  expect_equal(sm$t, t)
  expect_equal(dim(sm$S), c(length(t), 3))
})


test_that("as_survmat accepts matrix with numeric dimnames", {
  t <- seq(0, 4, by = 1)
  S <- matrix(runif(length(t) * 3), nrow = length(t), ncol = 3)
  rownames(S) <- as.character(t)

  sm <- as_survmat(S)
  expect_equal(sm$t, t)
  expect_equal(dim(sm$S), c(length(t), 3))
})


test_that("as_survmat accepts list(time, surv) and list(times, survival)", {
  t <- seq(0, 3, by = 1)
  surv1 <- matrix(runif(length(t) * 2), nrow = length(t), ncol = 2)
  surv2 <- t(surv1)

  sm1 <- as_survmat(list(time = t, surv = surv1))
  sm2 <- as_survmat(list(times = t, survival = surv2))

  expect_equal(sm1$t, t)
  expect_equal(sm2$t, t)
  expect_equal(dim(sm1$S), c(length(t), 2))
  expect_equal(dim(sm2$S), c(length(t), 2))
})


test_that("as_survmat parses survdnn-style data.frame", {
  df <- data.frame(
    `t=0` = runif(5),
    `t=1` = runif(5),
    `t=2` = runif(5),
    check.names = FALSE
  )

  sm <- as_survmat(df)
  expect_equal(sm$t, c(0, 1, 2))
  expect_equal(dim(sm$S), c(3, 5))
})


test_that("as_survmat accepts data.frame with explicit t when no t= names", {
  t <- seq(0, 3, by = 1)
  df <- as.data.frame(matrix(runif(12), nrow = 3, ncol = length(t)))

  sm <- as_survmat(df, t = t)
  expect_equal(sm$t, t)
  expect_equal(dim(sm$S), c(length(t), 3))
})


test_that("as_survmat errors when time cannot be determined", {
  M <- matrix(runif(12), nrow = 3, ncol = 4)
  df <- as.data.frame(M)

  expect_error(as_survmat(M), "Could not determine a valid time grid")
  expect_error(as_survmat(df), "For data.frame inputs")
})


test_that("as_survmat auto uses deterministic defaults on square objects", {
  t <- seq(0, 3, by = 1)

  M <- matrix(runif(16), nrow = 4, ncol = 4)
  sm_m <- as_survmat(M, t = t)
  expect_equal(dim(sm_m$S), c(4, 4))

  df <- as.data.frame(M)
  sm_df <- as_survmat(df, t = t)
  expect_equal(dim(sm_df$S), c(4, 4))
})


test_that("as_survmat strict mismatch errors for explicit vs encoded time", {
  df <- data.frame(
    `t=0` = runif(3),
    `t=1` = runif(3),
    `t=2` = runif(3),
    check.names = FALSE
  )

  expect_error(
    as_survmat(df, t = c(0, 1, 3), strict = TRUE),
    "does not match times encoded"
  )
})

