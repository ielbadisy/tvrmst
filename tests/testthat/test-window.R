wtrue <- function(l, s, w) exp(-l * s) * (1 - exp(-l * w)) / l

test_that("rmstwindow matches closed-form exponential values", {
  f <- make_final_scope_fixture()
  s <- c(0, 0.73, 2, 3.4); w <- c(1, 1, 2.5, 1.6)
  r <- rmstwindow(f$xA, s, w)
  tr <- vapply(seq_along(s), function(j) wtrue(f$lambdaA, s[j], w[j]), numeric(f$nA))
  expect_lt(max(abs(r$individual - tr)), 1e-3)
  expect_equal(r$mean, colMeans(r$individual))
  expect_equal(r$by$estimate, r$mean)
})

test_that("on-grid windows are increments of rmstdynamic, s = 0 gives RMST", {
  f <- make_final_scope_fixture()
  R <- rmstdynamic(f$x_all, by = NULL)$individual
  r <- rmstwindow(f$x_all, s = c(0, 1, 2), w = c(2, 1, 3), by = NULL)
  ti <- function(t) which(abs(f$time - t) < 1e-9)
  expect_equal(r$individual[, 1], R[, ti(2)])
  expect_equal(r$individual[, 2], R[, ti(2)] - R[, ti(1)])
  expect_equal(r$individual[, 3], R[, ti(5)] - R[, ti(2)])
})

test_that("rmstcond: memoryless exponential and survivor-weighted mean", {
  f <- make_final_scope_fixture()
  s <- c(0, 1.25, 3); w <- 1.5
  r <- rmstcond(f$xA, s, w)
  tr <- outer(f$lambdaA, s, function(l, s) (1 - exp(-l * w)) / l)
  expect_lt(max(abs(r$individual - tr)), 1e-3)
  pop <- vapply(s, function(sj) mean(wtrue(f$lambdaA, sj, w)) / mean(exp(-f$lambdaA * sj)), 1)
  expect_lt(max(abs(r$mean - pop)), 1e-3)
  expect_gt(max(abs(r$mean - colMeans(r$individual))), 1e-3)
})

test_that("window and conditional deltas, bootstrap estimands", {
  f <- make_final_scope_fixture()
  s <- c(0, 1, 2); w <- 1
  dw <- rmstwindowdelta(f$xA, f$xB, s, w)
  dc <- rmstconddelta(f$xA, f$xB, s, w)
  expect_equal(dw$delta, rmstwindow(f$xB, s, w)$mean - rmstwindow(f$xA, s, w)$mean)
  expect_equal(dc$delta, rmstcond(f$xB, s, w)$mean - rmstcond(f$xA, s, w)$mean)
  expect_equal(dw$delta[1], rmstdelta(f$xA, f$xB, tau = 1)$at_tau$delta)
  b <- bootrmstdelta(f$xA, f$xB, R = 20, seed = 1, estimand = "conditional", s = s, w = w)
  expect_equal(b$estimate, dc$delta)
  expect_true(all(b$lo <= b$estimate & b$estimate <= b$hi))
  expect_equal(b$windows, data.frame(s = s, w = w))
  expect_equal(bootrmstdelta(f$xA, f$xB, R = 5, seed = 1)$estimate, rmstdelta(f$xA, f$xB)$delta)
})

test_that("windows beyond the grid or with bad s, w are errors", {
  f <- make_final_scope_fixture()
  expect_error(rmstwindow(f$xA, 4.5, 1), "last grid time")
  expect_error(rmstcond(f$xA, -1, 1), "s >= 0")
  expect_length(rmstwindow(f$xA, 4, 1)$mean, 1)
})

test_that("conditional RMST stays accurate when S(s) is tiny", {
  time <- seq(0, 5, by = 0.01); l <- 20
  r <- rmstcond(assurvmat(matrix(exp(-l * time), 1), time), c(1.5, 3), 1)
  expect_lt(max(abs(r$individual / ((1 - exp(-l)) / l) - 1)), 0.01)
})
