# tvrmst

`tvrmst` computes RMST and time-varying RMST summaries from survival
predictions on a common time grid.

## Input contract

All compute/bootstrap functions use:

-   `time`: strictly increasing numeric vector, length `m >= 2`
-   `S`: numeric `matrix` or `data.frame`, shape `n_units x m` (rows =
    units, cols = time)

No vector input, no auto-transpose.

## Minimal API

<table>
<colgroup>
<col style="width: 33%" />
<col style="width: 33%" />
<col style="width: 33%" />
</colgroup>
<thead>
<tr class="header">
<th>Task</th>
<th>Function</th>
<th>Output</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>RMST at fixed horizon</td>
<td><code>rmst_tau()</code></td>
<td>numeric vector</td>
</tr>
<tr class="even">
<td>RMST curve</td>
<td><code>rmst_curve()</code></td>
<td>data.frame</td>
</tr>
<tr class="odd">
<td>RMST delta curve</td>
<td><code>rmst_delta_curve()</code></td>
<td>data.frame</td>
</tr>
<tr class="even">
<td>Window RMST contrast</td>
<td><code>rmst_window()</code></td>
<td>data.frame or matrix</td>
</tr>
<tr class="odd">
<td>Conditional tvRMST</td>
<td><code>tvrmst_cond()</code></td>
<td>data.frame or matrix</td>
</tr>
<tr class="even">
<td>Conditional difference</td>
<td><code>tvrmst_diff()</code></td>
<td>data.frame or matrix</td>
</tr>
<tr class="odd">
<td>Bootstrap CI (replicates)</td>
<td><code>bootstrap_tvrmst_diff_reps()</code></td>
<td>data.frame</td>
</tr>
</tbody>
</table>

Plot helpers:

-   `plot_survival_curves()`
-   `plot_rmst_curve()`
-   `plot_rmst_delta()`
-   `plot_tvrmst()`
-   `plot_tvrmst_diff()`
-   `plot_rmst_individual()`
-   `plot_rmst_mean()`

## Classical workflow (KM, simulated non-PH trial)

    library(tvrmst)
    library(survival)
    set.seed(1)

    n <- 400
    arm <- rbinom(n, 1, 0.5)
    # non-PH style simulation: treatment better early, attenuates later
    rate <- ifelse(arm == 1, 0.09 + 0.02 * runif(n), 0.12 + 0.01 * runif(n))
    Tevent <- rexp(n, rate = rate)
    C <- runif(n, 2, 12)
    time_obs <- pmin(Tevent, C)
    status <- as.integer(Tevent <= C)

    dat <- data.frame(time = time_obs, status = status, arm = factor(arm, labels = c("Control", "Treatment")))
    time_grid <- seq(0, 10, by = 0.25)

    km_to_grid <- function(df_arm, grid) {
      fit <- survfit(Surv(time, status) ~ 1, data = df_arm)
      sf <- summary(fit, times = grid, extend = TRUE)
      sf$surv
    }

    S0 <- km_to_grid(subset(dat, arm == "Control"), time_grid)
    S1 <- km_to_grid(subset(dat, arm == "Treatment"), time_grid)

    S0_mat <- rbind(Control = S0)
    S1_mat <- rbind(Treatment = S1)
    S_both <- rbind(Control = S0, Treatment = S1)

    rmst_tau(S_both, time_grid, tau = 6)

    rc0 <- rmst_curve(S0_mat, time_grid, statistic = "mean", probs = NULL)
    rc1 <- rmst_curve(S1_mat, time_grid, statistic = "mean", probs = NULL)
    rc_wide <- data.frame(tau = time_grid, Control = rc0$estimate, Treatment = rc1$estimate)
    dr <- rmst_delta_curve(rc_wide, arm1 = "Treatment", arm0 = "Control")

    s_grid <- seq(0, 7, by = 0.5)
    tau_win <- 2
    window_w <- rmst_window(S1_mat, S0_mat, time_grid, s_grid, tau = tau_win, statistic = "mean")
    mu0 <- tvrmst_cond(S0_mat, time_grid, s_grid, tau = tau_win, statistic = "mean")
    mu1 <- tvrmst_cond(S1_mat, time_grid, s_grid, tau = tau_win, statistic = "mean")
    delta_c <- tvrmst_diff(S1_mat, S0_mat, time_grid, s_grid, tau = tau_win, statistic = "mean")

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_survival_curves(S0_mat, S1_mat, time_grid, labels = c("Control", "Treatment"), show = "mean")
      plot_rmst_curve(rc0, title = "RMST curve (Control)")
      plot_rmst_delta(dr, title = "Delta RMST")
      plot_tvrmst(mu0, mu1, labels = c("Control", "Treatment"), title = "Conditional tvRMST")
      plot_tvrmst_diff(delta_c, title = "Delta conditional tvRMST")
    }

## Bootstrap CI from KM replicates

    set.seed(2)

    mk_rep <- function(df, grid) {
      idx <- sample(seq_len(nrow(df)), size = nrow(df), replace = TRUE)
      d <- df[idx, , drop = FALSE]
      S0_r <- km_to_grid(subset(d, arm == "Control"), grid)
      S1_r <- km_to_grid(subset(d, arm == "Treatment"), grid)
      list(S1 = rbind(Treatment = S1_r), S0 = rbind(Control = S0_r))
    }

    reps <- replicate(200, mk_rep(dat, time_grid), simplify = FALSE)

    boot_df <- bootstrap_tvrmst_diff_reps(
      reps = reps,
      time = time_grid,
      s_grid = s_grid,
      tau = tau_win,
      conf = 0.95,
      statistic = "mean"
    )

    if (requireNamespace("ggplot2", quietly = TRUE)) {
      plot_tvrmst_diff(boot_df, title = "Bootstrap CI for Delta conditional tvRMST")
    }

## Short survdnn demo (individualized RMST curves)

    if (requireNamespace("survdnn", quietly = TRUE)) {
      # Example sketch; adapt formula/data to your survdnn setup
      set.seed(3)
      id_train <- sample(seq_len(nrow(dat)), floor(0.7 * nrow(dat)))
      train_data <- dat[id_train, ]
      test_data <- dat[-id_train, ]

      fit <- survdnn::survdnn(
        formula = Surv(time, status) ~ arm,
        data = train_data
      )

      pred_df <- predict(fit, newdata = test_data, type = "survival", times = time_grid)
      RMST_dyn <- rmst_dynamic(pred_df, time_grid)

      if (requireNamespace("ggplot2", quietly = TRUE)) {
        plot_rmst_individual(RMST_dyn, group = test_data$arm)
      }
    }

## References

-   Uno H, Claggett B, Tian L, et al. Moving beyond the hazard ratio in
    quantifying the between-group difference in survival analysis. *J
    Clin Oncol.* 2014.
-   Royston P, Parmar MKB. Restricted mean survival time: an alternative
    to the hazard ratio for the design and analysis of randomized trials
    with a time-to-event outcome. *BMC Med Res Methodol.* 2013.
