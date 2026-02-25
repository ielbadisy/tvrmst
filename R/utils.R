# internal helpers

.stop <- function(...) stop(sprintf(...), call. = FALSE)

.check_time <- function(time) {
  if (!is.numeric(time) || length(time) < 2) .stop("`time` must be numeric, length >= 2.")
  if (any(!is.finite(time))) .stop("`time` contains non-finite values.")
  if (any(diff(time) <= 0)) .stop("`time` must be strictly increasing.")
  invisible(TRUE)
}

.coerce_unit_time <- function(S, time, name = "S") {
  validate_surv_input(S, time, name = name)
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

# linear interpolation for RMST-like trajectories
# returns nrow(Y) x length(xout)
.linear_at <- function(t, Y, xout) {
  xout <- as.numeric(xout)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  out <- vapply(
    seq_len(nrow(Y)),
    function(i) {
      stats::approx(t, Y[i, ], xout = xout, rule = 2)$y
    },
    FUN.VALUE = numeric(length(xout))
  )
  if (is.vector(xout)) {
    dim(out) <- c(length(xout), nrow(Y))
  }
  t(out)
}
