library(testthat)

test_that("survmat constructors and binders behave correctly", {
  f <- make_final_scope_fixture()

  expect_s3_class(f$xA, "survmat")
  expect_s3_class(f$xB, "survmat")
  expect_equal(nrow(f$xA$S), f$nA)
  expect_equal(nrow(f$xB$S), f$nB)
  expect_equal(length(f$xA$time), length(f$time))
  expect_true(all(f$xA$S >= 0 & f$xA$S <= 1))
  expect_true(all(f$xB$S >= 0 & f$xB$S <= 1))

  expect_equal(nobs_survmat(f$xA), f$nA)
  expect_equal(nobs_survmat(f$xB), f$nB)

  expect_s3_class(f$x_all, "survmat")
  expect_equal(nobs_survmat(f$x_all), f$nA + f$nB)
  expect_equal(length(f$x_all$group), f$nA + f$nB)
  expect_true(all(levels(f$x_all$group) %in% c("A", "B")))
})
