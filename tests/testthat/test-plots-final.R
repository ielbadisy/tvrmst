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
