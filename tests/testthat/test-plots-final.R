library(testthat)

test_that("plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()
  tau_test <- c(1.2, 2.7, 4.6)

  resAll <- rmstdynamic(f$x_all)
  d <- rmstdelta(f$xA, f$xB, tau = tau_test)
  boot_d <- bootrmstdelta(f$xA, f$xB, R = 50, seed = 1)

  p_ind <- plotrmstindividualbygroup(resAll, f$x_all$group, n_show_per_group = 30)
  p_two <- plotrmsttwoarms(f$xA, f$xB, labels = c("A", "B"))
  p_del <- plotdeltacurve(d$time, d$delta, title = "ΔRMST(τ)=B-A", xlab = "τ", ylab = "ΔRMST")
  p_boot <- plotbootcurve(boot_d, title = "Bootstrap CI: ΔRMST(τ)", xlab = "τ", ylab = "ΔRMST")

  expect_s3_class(p_ind, "ggplot")
  expect_s3_class(p_two, "ggplot")
  expect_s3_class(p_del, "ggplot")
  expect_s3_class(p_boot, "ggplot")
})

test_that("plot helpers use theme_bw as the base theme", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()
  p_two <- plotrmsttwoarms(f$xA, f$xB)
  expected <- ggplot2::theme_bw(base_size = 12)

  expect_equal(
    ggplot2::calc_element("panel.background", p_two$theme),
    ggplot2::calc_element("panel.background", expected)
  )
  expect_equal(
    ggplot2::calc_element("panel.border", p_two$theme),
    ggplot2::calc_element("panel.border", expected)
  )
  expect_equal(
    ggplot2::calc_element("axis.line", p_two$theme),
    ggplot2::calc_element("axis.line", expected)
  )
})

test_that("plotrmsttwoarms supports custom labels, colors, and optional title", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()
  custom_colors <- c("#004E89", "#E07A1F")

  p_two <- plotrmsttwoarms(
    f$xA,
    f$xB,
    labels = c("Control", "Treatment"),
    title = NULL,
    xlab = "Follow-up time",
    ylab = "Restricted mean survival time",
    curve_colors = custom_colors
  )

  expect_null(p_two$labels$title)
  expect_identical(p_two$labels$x, "Follow-up time")
  expect_identical(p_two$labels$y, "Restricted mean survival time")

  color_scale <- p_two$scales$get_scales("colour")
  expect_equal(unname(color_scale$palette(2)), custom_colors)
})

test_that("plot helpers support axis rescaling and unit labels", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()
  d <- rmstdelta(f$xA, f$xB)
  boot_d <- bootrmstdelta(f$xA, f$xB, R = 30, seed = 1)

  p_two <- plotrmsttwoarms(
    f$xA,
    f$xB,
    x_scale = 2,
    y_scale = 0.5,
    x_unit = "months",
    y_unit = "months"
  )
  expect_equal(p_two$data$tau[1], f$xA$time[1] / 2)
  expect_equal(p_two$data$estimate[1], rmstdynamic(f$xA, by = NULL)$mean[1] / 0.5)
  expect_identical(p_two$labels$x, "Time (tau), months")
  expect_identical(p_two$labels$y, "RMST(tau), months")

  p_delta <- plotdeltacurve(
    d$time,
    d$delta,
    x_scale = 10,
    y_scale = 2,
    x_unit = "years",
    y_unit = "years"
  )
  expect_equal(p_delta$data$t[2], d$time[2] / 10)
  expect_equal(p_delta$data$delta[2], d$delta[2] / 2)
  expect_identical(p_delta$labels$x, "t (years)")
  expect_identical(p_delta$labels$y, "Delta (years)")

  p_boot <- plotbootcurve(
    boot_d,
    x_scale = 4,
    y_scale = 2,
    x_unit = "quarters",
    y_unit = "months"
  )
  expect_equal(p_boot$data$t[3], boot_d$time[3] / 4)
  expect_equal(p_boot$data$estimate[3], boot_d$estimate[3] / 2)
  expect_equal(p_boot$data$lo[3], boot_d$lo[3] / 2)
  expect_equal(p_boot$data$hi[3], boot_d$hi[3] / 2)
  expect_identical(p_boot$labels$x, "t (quarters)")
  expect_identical(p_boot$labels$y, "estimate (months)")
})

test_that("plot helpers reject invalid axis scales", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()

  expect_error(
    plotrmsttwoarms(f$xA, f$xB, x_scale = 0),
    "`x_scale` must be a single positive finite number."
  )
  expect_error(
    plotdeltacurve(1:3, c(0.1, 0.2, 0.3), y_scale = -1),
    "`y_scale` must be a single positive finite number."
  )
})
