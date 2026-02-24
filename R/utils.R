# internal helpers 

.stop <- function(...) stop(sprintf(...), call. = FALSE)

.check_time <- function(time) {
  if (!is.numeric(time) || length(time) < 2) .stop("`time` must be numeric, length >= 2.")
  if (any(!is.finite(time))) .stop("`time` contains non-finite values.")
  if (any(diff(time) <= 0)) .stop("`time` must be strictly increasing.")
  invisible(TRUE)
}

.coerce_unit_time <- function(S, time, name = "S") {
  .check_time(time)
  if (is.data.frame(S)) S <- as.matrix(S)
  if (!is.matrix(S)) .stop("`%s` must be a numeric matrix or data.frame.", name)
  if (!is.numeric(S)) .stop("`%s` must be numeric.", name)
  if (ncol(S) != length(time)) {
    .stop("`%s` must be n_units \u00d7 n_time (columns correspond to `time`). Ensure ncol(%s) == length(time).",
          name, name)
  }
  if (any(!is.finite(S))) .stop("`%s` contains non-finite values.", name)
  S
}

.check_survmat <- function(t, S, name = "S") {
  if (!is.numeric(t) || length(t) < 2) .stop("`t` must be numeric, length >= 2.")
  if (any(!is.finite(t))) .stop("`t` contains non-finite values.")
  if (!is.matrix(S)) .stop("`%s` must be a matrix (n_time x n_series).", name)
  if (nrow(S) != length(t)) .stop("nrow(%s) must equal length(t).", name)
  if (any(!is.finite(S))) .stop("`%s` contains non-finite values.", name)
  invisible(TRUE)
}

.sort_survmat <- function(t, S) {
  if (any(diff(t) <= 0)) .stop("`time` must be strictly increasing.")
  list(t = t, S = S)
}

.trapz_cols <- function(x, Y) {
  dt <- diff(x)
  YL <- Y[, -ncol(Y), drop = FALSE]
  YR <- Y[, -1, drop = FALSE]
  rowSums((YL + YR) / 2 * dt)
}

.cumtrapz_cols <- function(x, Y) {
  dt <- diff(x)
  YL <- Y[, -ncol(Y), drop = FALSE]
  YR <- Y[, -1, drop = FALSE]
  incr <- (YL + YR) / 2 * dt
  t(apply(incr, 1, cumsum))
}

# left-continuous step interpolation for unit-first matrices
# returns nrow(S) x length(xout)
.surv_at <- function(t, S, xout) {
  xout <- as.numeric(xout)
  if (!is.matrix(S)) S <- as.matrix(S)
  out <- vapply(
    seq_len(nrow(S)),
    function(i) {
      stats::approx(t, S[i, ], xout = xout, method = "constant", f = 0, rule = 2)$y
    },
    FUN.VALUE = numeric(length(xout))
  )
  if (is.vector(xout)) {
    dim(out) <- c(length(xout), nrow(S))
  }
  t(out)
}

# ensure grid includes time 0 (with S(0)=1) for integration from 0
# if time already has 0, leave as-is
.extend_to_zero <- function(t, S) {
  if (min(t) > 0) {
    t <- c(0, t)
    S <- cbind(rep(1, nrow(S)), S)
  }
  list(t = t, S = S)
}

# ensure grid includes tau endpoint for integration up to tau (constant extension)
.extend_to_tau <- function(t, S, tau) {
  if (max(t) < tau) {
    t <- c(t, tau)
    S <- cbind(S, S[, ncol(S), drop = FALSE])
  }
  list(t = t, S = S)
}
