validate_surv_input <- function(S, time, name = "S") {
  if (!is.numeric(time) || length(time) < 2 || any(!is.finite(time)) || any(diff(time) <= 0)) {
    .stop("`time` must be numeric, finite, strictly increasing, and length >= 2.")
  }

  if (!is.matrix(S) && !is.data.frame(S)) {
    .stop("`S` must be a numeric matrix or data.frame with rows=units and columns=time points.")
  }

  S <- as.matrix(S)

  if (!is.numeric(S)) {
    .stop("`%s` must be numeric.", name)
  }

  if (ncol(S) != length(time)) {
    .stop("%s must have ncol(%s) == length(time). Got %d vs %d.",
          name, name, ncol(S), length(time))
  }

  if (any(!is.finite(S))) {
    .stop("`%s` contains non-finite values.", name)
  }

  tol <- 1e-8
  if (any(S < -tol | S > 1 + tol, na.rm = TRUE)) {
    warning(sprintf("`%s` has values outside [0,1].", name), call. = FALSE)
  }

  S
}
