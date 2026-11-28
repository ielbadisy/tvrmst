# tvrmst 0.0.11

* New `rmstwindow()` and `rmstcond()`: window RMST `W(s, w)` (time alive within `[s, s + w]`) and conditional (landmark) RMST `C(s, w) = E[min(T - s, w) | T > s]`, per subject and for the population (mean curve). `s` and `w` need not lie on the grid.
* New `rmstwindowdelta()` and `rmstconddelta()`: B minus A contrasts of these estimands. The conditional contrast compares the survivors at `s` in each arm and is documented as prognostic, not a treatment effect.
* `bootrmstdelta()` gains `estimand`, `s` and `w` to bootstrap the window and conditional contrasts. The default (`"cumulative"`) keeps the previous results.

# tvrmst 0.0.10

* `rmstdynamic()`, `rmstdelta()`, `bootrmstdelta()`, `nobssurvmat()`, `bindsurvmat()` and `plotrmsttwoarms()` now also accept a numeric survival matrix carrying its grid in `attr(x, "times")`, such as the `survmat` matrix returned by survalis `predict()`. Before, this input failed because both packages use the class name `survmat`. `assurvmat()` now stores `S` as a plain numeric matrix.
* `bootrmstdelta()` gains `paired`. With `paired = TRUE` (same subjects in both arms, e.g. counterfactual predictions), subjects are resampled jointly, keeping the within-subject correlation. The default `paired = FALSE` keeps the independent-arm resampling and its results.

# tvrmst 0.0.9

Exported functions are renamed to glued lowercase (no underscores). The old names remain as deprecated wrappers for one release (see `?"tvrmst-deprecated"`). The `survmat` class name is unchanged.

| Old | New |
|---|---|
| `as_survmat()` | `assurvmat()` |
| `as_survprob_matrix()` | `assurvprobmatrix()` |
| `bind_survmat()` | `bindsurvmat()` |
| `boot_rmst_delta()` | `bootrmstdelta()` |
| `bootstrap_curve()` | `bootstrapcurve()` |
| `nobs_survmat()` | `nobssurvmat()` |
| `plot_boot_curve()` | `plotbootcurve()` |
| `plot_delta_curve()` | `plotdeltacurve()` |
| `plot_rmst_individual_by_group()` | `plotrmstindividualbygroup()` |
| `plot_rmst_two_arms()` | `plotrmsttwoarms()` |
| `rmst_delta()` | `rmstdelta()` |
| `rmst_dynamic()` | `rmstdynamic()` |

