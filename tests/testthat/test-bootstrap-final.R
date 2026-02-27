library(testthat)

test_that("bootstrap_curve scalar smoke test passes", {
  f <- make_final_scope_fixture()
  resample_idx <- getFromNamespace("resample_idx", "tvrmst")

  boot_simple <- bootstrap_curve(function(r) {
    if (is.null(r)) {
      mean(rmst_dynamic(f$xA, by = NULL)$mean)
    } else {
      idx <- resample_idx(f$nA, f$xA$group)
      xb <- as_survmat(f$xA$S[idx, , drop = FALSE], f$time, group = f$xA$group[idx])
      mean(rmst_dynamic(xb, by = NULL)$mean)
    }
  }, R = 100, seed = 1)

  expect_equal(length(boot_simple$estimate), 1)
  expect_gte(boot_simple$hi, boot_simple$lo)
})

test_that("boot_rmst_delta returns aligned confidence bands", {
  f <- make_final_scope_fixture()

  boot_d <- boot_rmst_delta(f$xA, f$xB, R = 200, seed = 1)

  expect_equal(length(boot_d$estimate), length(f$time))
  expect_equal(length(boot_d$lo), length(f$time))
  expect_equal(length(boot_d$hi), length(f$time))
  expect_true(all(boot_d$hi >= boot_d$lo))
})
