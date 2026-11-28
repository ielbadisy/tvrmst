#' Window and conditional (landmark) RMST
#'
#' `rmstwindow()` computes the window RMST `W_i(s, w) = integral_s^(s+w) S_i(u) du`,
#' the expected time alive within `[s, s + w]`. `rmstcond()` computes the
#' conditional (landmark) RMST `C_i(s, w) = W_i(s, w) / S_i(s)
#' = E[min(T - s, w) | T > s]`, the expected time alive in the next `w` time
#' units for a subject alive at `s`.
#'
#' Survival is linearly interpolated between grid points, so the integrals are
#' exact for the trapezoidal rule used by [rmstdynamic()] and `s`, `s + w` need
#' not lie on the grid. `s` and `w` are recycled to a common length; each pair
#' defines one window.
#'
#' The `mean` element is the population quantity of the mean curve
#' `Sbar = colMeans(S)`. For the window RMST it equals the mean of the individual
#' values. For the conditional RMST it is `mean(W_i) / mean(S_i(s))`, the
#' conditional RMST of the population still alive at `s`, which is *not* the mean
#' of the individual `C_i`: survivors at `s` are weighted by their probability of
#' being alive. A subject with `S_i(s) = 0` has an undefined `C_i` (`NaN`).
#'
#' @param x A `survmat` object, or a survival matrix with a `"times"` attribute
#'   (see [rmstdynamic()]).
#' @param s Numeric landmark(s), `s >= 0`.
#' @param w Numeric window width(s), `w > 0`, with `s + w <= max(x$time)`.
#' @param by Optional grouping variable. Defaults to `x$group`.
#'
#' @return A list with `windows` (data frame of `s`, `w`), `individual` (matrix,
#'   subjects x windows), `mean` (population value per window) and, with `by`,
#'   `by` (data frame of `group`, `s`, `w`, `estimate`).
#' @examples
#' time <- seq(0, 5, by = 0.1)
#' S <- rbind(exp(-0.2 * time), exp(-0.5 * time))
#' x <- assurvmat(S, time, group = c("A", "B"))
#' rmstwindow(x, s = 0:3, w = 1)$by
#' rmstcond(x, s = 0:3, w = 1)$by
#' @name rmstwindow
NULL

#' @rdname rmstwindow
#' @export
rmstwindow <- function(x, s, w, by = x$group) {
  x <- tosm(x)
  win(x, s, w, by, cond = FALSE)
}

#' @rdname rmstwindow
#' @export
rmstcond <- function(x, s, w, by = x$group) {
  x <- tosm(x)
  win(x, s, w, by, cond = TRUE)
}

#' Window and conditional RMST difference between two arms
#'
#' `rmstwindowdelta()` returns `W_B(s, w) - W_A(s, w)` and `rmstconddelta()`
#' returns `C_B(s, w) - C_A(s, w)`, both from the mean curves of each arm (see
#' [rmstwindow()]). The sign is B minus A, as in [rmstdelta()].
#'
#' For `s > 0` the conditional difference compares the subjects still alive at
#' `s` in each arm. These groups differ when the arms have different survival
#' before `s`, even under randomisation, so the conditional difference is a
#' prognostic contrast, not a treatment effect. Use the window difference for a
#' time-localised treatment effect.
#'
#' @inheritParams rmstdelta
#' @inheritParams rmstwindow
#'
#' @return A data frame with `s`, `w`, `est_A`, `est_B` and `delta`.
#' @examples
#' time <- seq(0, 5, by = 0.1)
#' xA <- assurvmat(matrix(exp(-0.5 * time), 1), time)
#' xB <- assurvmat(matrix(exp(-0.3 * time), 1), time)
#' rmstwindowdelta(xA, xB, s = 0:3, w = 1)
#' rmstconddelta(xA, xB, s = 0:3, w = 1)
#' @name rmstwindowdelta
NULL

#' @rdname rmstwindowdelta
#' @export
rmstwindowdelta <- function(xA, xB, s, w) windelta(xA, xB, s, w, cond = FALSE)

#' @rdname rmstwindowdelta
#' @export
rmstconddelta <- function(xA, xB, s, w) windelta(xA, xB, s, w, cond = TRUE)

win <- function(x, s, w, by, cond) {
  stopifnot(inherits(x, "survmat"))
  g <- wingrid(x$time, s, w)
  Ss <- survat(x$time, x$S, g$s)
  W <- cumat(x$time, x$S, g$s + g$w) - cumat(x$time, x$S, g$s)
  ind <- if (cond) W / Ss else W
  m <- function(idx) {
    Wm <- colMeans(W[idx, , drop = FALSE])
    if (cond) Wm / colMeans(Ss[idx, , drop = FALSE]) else Wm
  }
  out <- list(estimand = if (cond) "rmstcond" else "rmstwindow", windows = g,
              individual = ind, mean = m(seq_len(nrow(W))))
  if (!is.null(by)) {
    by <- as.factor(by)
    out$by <- do.call(rbind, lapply(levels(by), function(l)
      data.frame(group = l, g, estimate = m(which(by == l)))))
  }
  out
}

windelta <- function(xA, xB, s, w, cond) {
  xA <- tosm(xA); xB <- tosm(xB)
  stopifnot(inherits(xA, "survmat"), inherits(xB, "survmat"))
  a <- win(xA, s, w, NULL, cond)
  b <- win(xB, s, w, NULL, cond)
  data.frame(a$windows, est_A = a$mean, est_B = b$mean, delta = b$mean - a$mean)
}

wingrid <- function(time, s, w) {
  g <- data.frame(s = as.numeric(s), w = as.numeric(w))
  if (any(g$s < 0) || any(g$w <= 0)) stop("Need `s >= 0` and `w > 0`.", call. = FALSE)
  if (any(g$s + g$w > max(time) + 1e-12)) {
    stop("`s + w` must not exceed the last grid time (", max(time), ").", call. = FALSE)
  }
  g
}

# S linearly interpolated at t (rows = subjects, columns = t)
survat <- function(time, S, t) {
  k <- findInterval(t, time, rightmost.closed = TRUE)
  f <- (t - time[k]) / (time[k + 1L] - time[k])
  f[!is.finite(f)] <- 0
  k1 <- pmin(k + 1L, length(time))
  S[, k, drop = FALSE] + sweep(S[, k1, drop = FALSE] - S[, k, drop = FALSE], 2, f, `*`)
}

# integral_0^t of linearly interpolated S, exact for the trapezoidal rule
cumat <- function(time, S, t) {
  k <- findInterval(t, time, rightmost.closed = TRUE)
  R <- trapz_cum(time, S)
  St <- survat(time, S, t)
  R[, k, drop = FALSE] + sweep(S[, k, drop = FALSE] + St, 2, (t - time[k]) / 2, `*`)
}
