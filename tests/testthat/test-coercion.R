library(testthat)

test_that("as_survprob_matrix handles supported input shapes", {
  m <- matrix(runif(12), nrow = 3)
  expect_equal(as_survprob_matrix(m), m)

  df <- as.data.frame(m)
  out_df <- as_survprob_matrix(df)
  expect_true(is.matrix(out_df))
  expect_equal(dim(out_df), dim(m))

  out_list <- as_survprob_matrix(list(survival = df))
  expect_true(is.matrix(out_list))
  expect_equal(dim(out_list), dim(m))

  expect_error(as_survprob_matrix(list(foo = 1)), "Unsupported predict\\(\\) output")
})
