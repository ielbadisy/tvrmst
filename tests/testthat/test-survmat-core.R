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

  expect_equal(nobssurvmat(f$xA), f$nA)
  expect_equal(nobssurvmat(f$xB), f$nB)

  expect_s3_class(f$x_all, "survmat")
  expect_equal(nobssurvmat(f$x_all), f$nA + f$nB)
  expect_equal(length(f$x_all$group), f$nA + f$nB)
  expect_true(all(levels(f$x_all$group) %in% c("A", "B")))
})

test_that("assurvmat prepends t = 0 with S = 1 and messages", {
  expect_message(x <- assurvmat(rbind(c(0.9, 0.8)), c(1, 2)), "prepending")
  expect_equal(x$time, c(0, 1, 2))
  expect_equal(x$S[1, ], c(1, 0.9, 0.8))
  expect_silent(assurvmat(rbind(c(1, 0.8)), c(0, 2)))
})

test_that("assurvmat rejects invalid survival curves", {
  expect_error(assurvmat(rbind(c(0.95, 0.8)), c(0, 2)), "S\\(0\\) must equal 1")
  expect_error(assurvmat(rbind(c(1, 0.6, 0.7)), c(0, 1, 2)), "non-increasing")
  expect_error(assurvmat(rbind(c(1, 1.2)), c(0, 1)), "\\[0, 1\\]")
  expect_error(assurvmat(rbind(c(1, -0.1)), c(0, 1)), "\\[0, 1\\]")
  # numerical noise within tolerance is accepted and clamped
  x <- assurvmat(rbind(c(1 + 1e-12, 0.8, 0.8 + 1e-12)), c(0, 1, 2))
  expect_equal(x$S[1, 1], 1)
  expect_true(all(x$S <= 1))
})

test_that("rmstdelta is B minus A", {
  time <- c(0, 1, 2)
  xA <- assurvmat(rbind(c(1, 0.5, 0.2)), time)
  xB <- assurvmat(rbind(c(1, 0.9, 0.8)), time)
  expect_true(all(rmstdelta(xA, xB)$delta >= 0))
})
