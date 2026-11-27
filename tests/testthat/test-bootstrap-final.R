library(testthat)

test_that("bootstrapcurve scalar smoke test passes", {
  f <- make_final_scope_fixture()
  resample_idx <- getFromNamespace("resample_idx", "tvrmst")

  boot_simple <- bootstrapcurve(function(r) {
    if (is.null(r)) {
      mean(rmstdynamic(f$xA, by = NULL)$mean)
    } else {
      idx <- resample_idx(f$nA, f$xA$group)
      xb <- assurvmat(f$xA$S[idx, , drop = FALSE], f$time, group = f$xA$group[idx])
      mean(rmstdynamic(xb, by = NULL)$mean)
    }
  }, R = 100, seed = 1)

  expect_equal(length(boot_simple$estimate), 1)
  expect_gte(boot_simple$hi, boot_simple$lo)
})

test_that("bootrmstdelta returns aligned confidence bands", {
  f <- make_final_scope_fixture()

  boot_d <- bootrmstdelta(f$xA, f$xB, R = 200, seed = 1)

  expect_equal(length(boot_d$estimate), length(f$time))
  expect_equal(length(boot_d$lo), length(f$time))
  expect_equal(length(boot_d$hi), length(f$time))
  expect_true(all(boot_d$hi >= boot_d$lo))
})

test_that("bootrmstdelta(paired = TRUE) resamples subjects jointly", {
  fx <- make_final_scope_fixture()
  xA <- fx$xA
  xB <- assurvmat(xA$S^0.7, fx$time, group = rep("B", fx$nA))
  bp <- bootrmstdelta(xA, xB, R = 50, seed = 2, paired = TRUE)
  set.seed(2)
  k <- ncol(xA$S); d <- (trapz_cum(fx$time, xB$S) - trapz_cum(fx$time, xA$S))[, k]
  ref <- replicate(50, { i <- sample.int(fx$nA, fx$nA, replace = TRUE); mean(d[i]) })
  expect_equal(bp$lo[k], unname(quantile(ref, 0.025)))
  bi <- bootrmstdelta(xA, xB, R = 50, seed = 2)
  expect_lt(bp$hi[k] - bp$lo[k], (bi$hi[k] - bi$lo[k]) / 2)
  expect_equal(bi, bootrmstdelta(xA, xB, R = 50, seed = 2, paired = FALSE))
  expect_error(bootrmstdelta(fx$xA, fx$xB, R = 5, paired = TRUE), "same number of rows")
})
