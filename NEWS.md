# tvrmst 0.0.8

- `as_survmat()` now validates that each row of `S` is a survival curve: values in `[0, 1]` and non-increasing over time (numerical tolerance `sqrt(.Machine$double.eps)`). It errors when `t = 0` is present with `S(0) != 1` instead of silently resetting it to 1.

# tvrmst 0.0.7

- `as_survmat()` now enforces a time grid starting at 0 with `S = 1`: it prepends `t = 0`, `S = 1` (with a message) when the grid starts after 0, and resets `S(0)` to 1 with a message otherwise.
- Documented that `boot_rmst_delta()` resamples fixed prediction matrices and does not capture model uncertainty; a refit bootstrap is recommended.
- Documented the sign of `rmst_delta()` (B minus A).
- Fixed author name casing ("El Badisy") in LICENSE and LICENSE.md.

# tvrmst 0.0.6

- Added `rmst_dynamic()` for individual dynamic RMST curves over the grid.
- Added plot helpers: `plot_rmst_individual()` and `plot_rmst_mean()`.
- Added RMST curve plot helpers: `plot_rmst_curve()` and `plot_rmst_delta()`.
- Updated example script and technical report.
- Added `rmst_dynamic_curves()` and `plot_rmst_dynamic_by_cluster()` for
  model-agnostic individual dynamic RMST trajectories on `(t, S)` inputs.
