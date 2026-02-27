library(testthat)

test_that("plot helpers return ggplot objects", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()
  tau_test <- c(1.2, 2.7, 4.6)

  resAll <- rmst_dynamic(f$x_all)
  d <- rmst_delta(f$xA, f$xB, tau = tau_test)
  boot_d <- boot_rmst_delta(f$xA, f$xB, R = 50, seed = 1)

  p_ind <- plot_rmst_individual_by_group(resAll, f$x_all$group, n_show_per_group = 30)
  p_two <- plot_rmst_two_arms(f$xA, f$xB, labels = c("A", "B"))
  p_del <- plot_delta_curve(d$time, d$delta, title = "ΔRMST(τ)=B-A", xlab = "τ", ylab = "ΔRMST")
  p_boot <- plot_boot_curve(boot_d, title = "Bootstrap CI: ΔRMST(τ)", xlab = "τ", ylab = "ΔRMST")

  expect_s3_class(p_ind, "ggplot")
  expect_s3_class(p_two, "ggplot")
  expect_s3_class(p_del, "ggplot")
  expect_s3_class(p_boot, "ggplot")
})

test_that("plot_rmst_two_arms supports custom labels, colors, and optional title", {
  skip_if_not_installed("ggplot2")

  f <- make_final_scope_fixture()
  custom_colors <- c("#004E89", "#E07A1F")

  p_two <- plot_rmst_two_arms(
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
