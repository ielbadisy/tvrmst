library(testthat)

test_that("rmst_dynamic supports individual, grouped, and at_tau checks", {
  f <- make_final_scope_fixture()
  tau_test <- c(1.2, 2.7, 4.6)

  resA <- rmst_dynamic(f$xA, tau = tau_test)
  resB <- rmst_dynamic(f$xB, tau = tau_test)
  resAll <- rmst_dynamic(f$x_all)

  expect_false(is.null(resA$individual))
  expect_false(is.null(resB$individual))
  expect_false(is.null(resAll$individual))
  expect_equal(length(resA$mean), length(f$time))
  expect_equal(length(resB$mean), length(f$time))
  expect_equal(nrow(resAll$individual), f$nA + f$nB)

  expect_true(all(apply(resAll$individual, 1, function(v) all(diff(v) >= -1e-10))))

  gA <- subset(resAll$by, group == "A")$estimate
  gB <- subset(resAll$by, group == "B")$estimate
  expect_lt(max(abs(gA - resA$mean)), 1e-10)
  expect_lt(max(abs(gB - resB$mean)), 1e-10)

  rmstA_true <- vapply(
    tau_test,
    function(tau) mean((1 - exp(-f$lambdaA * tau)) / f$lambdaA),
    numeric(1)
  )
  expect_lt(max(abs(resA$at_tau$estimate - rmstA_true)), 0.02)
})

test_that("rmst_delta returns valid curve and tau summaries", {
  f <- make_final_scope_fixture()
  tau_test <- c(1.2, 2.7, 4.6)

  d <- rmst_delta(f$xA, f$xB, tau = tau_test)

  expect_equal(length(d$delta), length(f$time))
  expect_equal(nrow(d$at_tau), length(tau_test))
  expect_gt(mean(d$delta), 0)
})
