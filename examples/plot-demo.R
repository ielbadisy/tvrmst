pkgload::load_all()

f <- make_final_scope_fixture()
tau_test <- c(1.2, 2.7, 4.6)

res_all <- rmst_dynamic(f$x_all)
d <- rmst_delta(f$xA, f$xB, tau = tau_test)
boot_d <- boot_rmst_delta(f$xA, f$xB, R = 50, seed = 1)

plots <- list(
  individual = plot_rmst_individual_by_group(
    res_all,
    f$x_all$group,
    n_show_per_group = 30,
    title = "Individual dynamic RMST by group"
  ),
  two_arms = plot_rmst_two_arms(
    f$xA,
    f$xB,
    labels = c("Control", "Treatment"),
    title = "Dynamic RMST",
    xlab = "Follow-up time",
    ylab = "Restricted mean survival time"
  ),
  delta = plot_delta_curve(
    d$time,
    d$delta,
    title = "Delta RMST",
    xlab = "Follow-up time",
    ylab = "RMST difference"
  ),
  bootstrap = plot_boot_curve(
    boot_d,
    title = "Bootstrap CI",
    xlab = "Follow-up time",
    ylab = "RMST difference"
  )
)

for (name in names(plots)) {
  message("Rendering plot: ", name)
  print(plots[[name]])
}
