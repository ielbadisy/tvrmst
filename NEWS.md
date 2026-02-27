# tvrmst 0.0.6

- Added `rmst_dynamic()` for individual dynamic RMST curves over the grid.
- Added plot helpers: `plot_rmst_individual()` and `plot_rmst_mean()`.
- Added RMST curve plot helpers: `plot_rmst_curve()` and `plot_rmst_delta()`.
- Updated example script and technical report.
- Added `rmst_dynamic_curves()` and `plot_rmst_dynamic_by_cluster()` for
  model-agnostic individual dynamic RMST trajectories on `(t, S)` inputs.
- Added an end-to-end example script at
  `inst/examples/plot_individual_dynamic_rmst.R` that saves a faceted cluster
  plot PNG.
