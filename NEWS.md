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

