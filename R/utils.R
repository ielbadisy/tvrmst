# internal helpers 

.stop <- function(...) stop(sprintf(...), call. = FALSE)

.check_survmat <- function(t, S, name = "S") {
  if (!is.numeric(t) || length(t) < 2) .stop("`t` must be numeric, length >= 2.")
  if (any(!is.finite(t))) .stop("`t` contains non-finite values.")
  if (!is.matrix(S)) .stop("`%s` must be a matrix (n_time x n_series).", name)
  if (nrow(S) != length(t)) .stop("nrow(%s) must equal length(t).", name)
  if (any(!is.finite(S))) .stop("`%s` contains non-finite values.", name)
  invisible(TRUE)
}

.sort_survmat <- function(t, S) {
  ord <- order(t)
  list(t = t[ord], S = S[ord, , drop = FALSE])
}

.trapz_cols <- function(x, Y) {
  dt <- diff(x)
  YL <- Y[-nrow(Y), , drop = FALSE]
  YR <- Y[-1, , drop = FALSE]
  colSums((YL + YR) / 2 * dt)
}

.cumtrapz_cols <- function(x, Y) {
  dt <- diff(x)
  YL <- Y[-nrow(Y), , drop = FALSE]
  YR <- Y[-1, , drop = FALSE]
  incr <- (YL + YR) / 2 * dt
  apply(incr, 2, cumsum)
}

# left-continuous step interpolation for matrix-valued y
# returns length(xout) x ncol(S)
.surv_at <- function(t, S, xout) {
  xout <- as.numeric(xout)
  if (!is.matrix(S)) S <- as.matrix(S)
  out <- vapply(
    seq_len(ncol(S)),
    function(j) {
      stats::approx(t, S[, j], xout = xout, method = "constant", f = 0, rule = 2)$y
    },
    FUN.VALUE = numeric(length(xout))
  )
  # vapply returns length(xout) x ncol(S) *if* simplified like this:
  if (is.vector(xout)) {
    dim(out) <- c(length(xout), ncol(S))
  }
  out
}

# ensure grid includes time 0 (with S(0)=1) for integration from 0
# if t already has 0, leave as-is
.extend_to_zero <- function(t, S) {
  if (min(t) > 0) {
    t <- c(0, t)
    S <- rbind(rep(1, ncol(S)), S)
  }
  list(t = t, S = S)
}

# ensure grid includes tau endpoint for integration up to tau (constant extension)
.extend_to_tau <- function(t, S, tau) {
  if (max(t) < tau) {
    t <- c(t, tau)
    S <- rbind(S, S[nrow(S), , drop = FALSE])
  }
  list(t = t, S = S)
}
