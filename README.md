
# tvrmst

Time-varying RMST from survival curves on a time grid. The package is
model-agnostic: you provide a time vector and survival matrices, and it
returns RMST-based quantities and bootstrap summaries.

## Installation

From a local checkout:

``` r
# in R
install.packages("tvrmst", repos = NULL, type = "source")
```

## Quick start

``` r
library(tvrmst)

# time grid
t <- seq(0, 10, by = 1)

# survival matrices (n_time x n_series)
S0 <- cbind(exp(-0.20 * t))
S1 <- cbind(exp(-0.15 * t))

# RMST at a horizon
rmst_tau(t, S0, tau = 5)

# RMST curve
rmst_curve(t, S0)

# conditional tvRMST by landmark time
s_grid <- seq(0, 8, by = 1)
mu0 <- tvrmst_cond(t, S0, s_grid, tau = 2)
mu1 <- tvrmst_cond(t, S1, s_grid, tau = 2)

# difference in conditional tvRMST
delta <- tvrmst_diff(t, S1, S0, s_grid, tau = 2)

# bootstrap confidence intervals (columns are units)
S0_units <- cbind(exp(-0.20 * t), exp(-0.22 * t), exp(-0.18 * t), exp(-0.21 * t), exp(-0.19 * t))
S1_units <- cbind(exp(-0.15 * t), exp(-0.16 * t), exp(-0.14 * t), exp(-0.155 * t), exp(-0.145 * t))
boot <- bootstrap_tvrmst_diff_cols(t, S1_units, S0_units, s_grid, tau = 2, R = 200)
```

## Functions (progressive guide)

This section walks through the full API from core inputs to plotting.

### 1) Inputs

All functions use a time grid `t` and survival matrices `S` with
dimension `n_time x n_series`. Each column is a series (arm, group, or
unit).

``` r
t <- seq(0, 10, by = 1)
S0 <- cbind(Control = exp(-0.20 * t))
S1 <- cbind(Treatment = exp(-0.15 * t))
```

### 2) Classical RMST

`rmst_tau()` gives RMST at a single horizon.  
`rmst_curve()` evaluates RMST along the grid.  
`rmst_dynamic()` returns RMST at every grid time (including 0) for each
series.  
`rmst_delta_curve()` is the difference between two series in an
`rmst_curve()` output.

``` r
rmst_tau(t, S0, tau = 5)

rc <- rmst_curve(t, S0)
rd <- rmst_dynamic(t, S0)

rc2 <- rmst_curve(t, cbind(Control = S0[, 1], Treatment = S1[, 1]))
rmst_delta_curve(rc2, "Treatment", "Control")
```

### 3) Windowed and time‑varying RMST

`rmst_window()` computes a windowed RMST difference on landmarks.  
`tvrmst_cond()` computes conditional tvRMST for each series.  
`tvrmst_diff()` computes the difference in conditional tvRMST between
two arms.

``` r
s_grid <- seq(0, 8, by = 1)

rmst_window(t, S1, S0, s_grid, tau = 2)

mu0 <- tvrmst_cond(t, S0, s_grid, tau = 2)
mu1 <- tvrmst_cond(t, S1, s_grid, tau = 2)
delta <- tvrmst_diff(t, S1, S0, s_grid, tau = 2)
```

### 4) Bootstrap confidence intervals

`bootstrap_tvrmst_diff_cols()` bootstraps by resampling columns
(units).  
`bootstrap_tvrmst_diff_reps()` aggregates a list of replicate survival
matrices.

``` r
S0_units <- cbind(exp(-0.20 * t), exp(-0.22 * t), exp(-0.18 * t), exp(-0.21 * t), exp(-0.19 * t))
S1_units <- cbind(exp(-0.15 * t), exp(-0.16 * t), exp(-0.14 * t), exp(-0.155 * t), exp(-0.145 * t))
bootstrap_tvrmst_diff_cols(t, S1_units, S0_units, s_grid, tau = 2, R = 200)

reps <- lapply(1:10, function(i) {
  list(
    S1 = cbind(exp(-(0.15 + runif(1, -0.02, 0.02)) * t)),
    S0 = cbind(exp(-(0.20 + runif(1, -0.02, 0.02)) * t))
  )
})
bootstrap_tvrmst_diff_reps(t, reps, s_grid, tau = 2)
```

### 5) Identity checks

`check_identities()` validates key numerical identities for correctness.

``` r
check_identities(t, S0, s_grid, tau = 2)
```

### 6) Plot helpers (requires ggplot2)

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_survival_curves(t, S0, S1)
  plot_rmst_curve(rc)
  plot_rmst_delta(rmst_delta_curve(rc2, "Treatment", "Control"))
  plot_rmst_individual(rd)
  plot_rmst_mean(rd, group = rep("All", ncol(S0)))
  plot_tvrmst(mu0, mu1)
  plot_tvrmst_diff(delta)
}
```
