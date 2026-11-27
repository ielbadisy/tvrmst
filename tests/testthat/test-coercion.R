library(testthat)

test_that("assurvprobmatrix handles supported input shapes", {
  m <- matrix(runif(12), nrow = 3)
  expect_equal(assurvprobmatrix(m), m)

  df <- as.data.frame(m)
  out_df <- assurvprobmatrix(df)
  expect_true(is.matrix(out_df))
  expect_equal(dim(out_df), dim(m))

  out_list <- assurvprobmatrix(list(survival = df))
  expect_true(is.matrix(out_list))
  expect_equal(dim(out_list), dim(m))

  expect_error(assurvprobmatrix(list(foo = 1)), "Unsupported predict\\(\\) output")
})

test_that("rmst functions accept a survival matrix with a times attribute", {
  fx <- make_final_scope_fixture()
  tm <- fx$time[-1]
  sm <- function(x) structure(x$S[, -1], times = tm, class = c("survmat", "matrix", "array"))
  mA <- sm(fx$xA); mB <- sm(fx$xB)
  expect_message(r <- rmstdynamic(mA, tau = 2), "prepending")
  expect_equal(r$mean, rmstdynamic(fx$xA, by = NULL)$mean)
  expect_null(r$by)
  d <- suppressMessages(rmstdelta(mA, mB, tau = 2))
  expect_equal(d$delta, rmstdelta(fx$xA, fx$xB)$delta)
  expect_equal(suppressMessages(nobssurvmat(mA)), fx$nA)
  expect_equal(suppressMessages(nobssurvmat(bindsurvmat(mA, mB))), fx$nA + fx$nB)
  b <- suppressMessages(bootrmstdelta(mA, mB, R = 5, seed = 1))
  expect_equal(b, bootrmstdelta(fx$xA, fx$xB, R = 5, seed = 1))
  expect_false(inherits(suppressMessages(tosm(mA))$S, "survmat"))
})
