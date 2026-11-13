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
