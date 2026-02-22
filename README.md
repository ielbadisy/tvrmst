---
output: github_document
---



# tvrmst

Time-varying RMST from survival curves on a time grid. The package is
model-agnostic: you provide a time vector and survival matrices, and it returns RMST-based quantities and bootstrap summaries.

## Installation

From a local checkout:

```r
# in R
install.packages("tvrmst", repos = NULL, type = "source")
```

## Quick start

```r
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

## Functions

Core functions include:

- `rmst_tau()`: RMST at a single horizon.
- `rmst_curve()`: RMST evaluated along the time grid.
- `rmst_dynamic()`: RMST at every grid time for each series (individual curves).
- `rmst_window()`: windowed RMST difference between two arms.
- `tvrmst_cond()`: conditional time-varying RMST.
- `tvrmst_diff()`: difference in conditional tvRMST between two arms.
- `bootstrap_tvrmst_diff_cols()` and `bootstrap_tvrmst_diff_reps()` for CIs.
- `check_identities()` to verify key RMST identities numerically.

Plot helpers (requires `ggplot2`):

- `plot_survival_curves()`
- `plot_rmst_curve()`
- `plot_rmst_delta()`
- `plot_rmst_individual()`
- `plot_rmst_mean()`
- `plot_tvrmst()`
- `plot_tvrmst_diff()`
