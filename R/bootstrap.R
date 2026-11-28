#' Generic bootstrap for curve-valued estimators
#'
#' @param est_fun Function of one argument `r`; called with `NULL` for the
#'   point estimate and with integers `1:R` for bootstrap replicates.
#' @param R Number of replicates.
#' @param conf Confidence level in `(0, 1)`.
#' @param seed Optional RNG seed.
#' @param keep_reps If `TRUE`, include replicate matrix in output.
#'
#' @return A list with `estimate`, `lo`, `hi`, and `conf`.
#' @examples
#' vals <- c(1, 2, 4, 8)
#' boot <- bootstrapcurve(
#'   function(r) {
#'     if (is.null(r)) {
#'       mean(vals)
#'     } else {
#'       mean(sample(vals, replace = TRUE))
#'     }
#'   },
#'   R = 20,
#'   seed = 1
#' )
#' boot$estimate
#' @export
bootstrapcurve <- function(est_fun, R = 300, conf = 0.95, seed = NULL, keep_reps = FALSE) {
  stopifnot(is.function(est_fun), R >= 2, conf > 0 && conf < 1)
  if (!is.null(seed)) set.seed(seed)

  est0 <- as.numeric(est_fun(NULL))
  p <- length(est0)
  reps <- matrix(NA_real_, R, p)

  for (r in seq_len(R)) reps[r, ] <- as.numeric(est_fun(r))

  alpha <- (1 - conf) / 2
  lo <- apply(reps, 2, stats::quantile, probs = alpha, na.rm = TRUE, names = FALSE)
  hi <- apply(reps, 2, stats::quantile, probs = 1 - alpha, na.rm = TRUE, names = FALSE)

  out <- list(estimate = est0, lo = lo, hi = hi, conf = conf)
  if (keep_reps) out$reps <- reps
  out
}

#' Bootstrap CI for dynamic RMST delta curve
#'
#' Resamples subjects (rows) of the fixed prediction matrices `xA` and `xB`,
#' conditional on the fitted model(s) that produced them. The resulting bands
#' therefore reflect sampling variability of the predicted curves only, and do
#' NOT capture model uncertainty (estimation of the model parameters). For
#' valid inference, use a refit bootstrap: resample the original data, refit
#' the model, recompute the survival matrices and [rmstdelta()] in each
#' replicate, then take percentile limits of the replicated delta curves.
#'
#' @param xA A `survmat` object for arm A, or a survival matrix with a
#'   `"times"` attribute (see [rmstdynamic()]).
#' @param xB As `xA`, for arm B.
#' @param R Number of bootstrap replicates.
#' @param conf Confidence level in `(0, 1)`.
#' @param seed Optional RNG seed.
#' @param paired If `TRUE`, row `i` of `xA` and row `i` of `xB` are the same
#'   subject (e.g. counterfactual predictions of one sample under two
#'   treatments): subjects are resampled jointly, the same rows in both arms,
#'   which keeps the within-subject correlation of the two curves. Requires
#'   `nrow(xA) == nrow(xB)`; strata are taken from `xA$group`. If `FALSE` (the
#'   default, for arms of genuinely different subjects), the two arms are
#'   resampled independently. Using `FALSE` on paired curves ignores the
#'   covariance and gives too wide bands.
#' @param estimand `"cumulative"` (default, the curve of [rmstdelta()]),
#'   `"window"` ([rmstwindowdelta()]) or `"conditional"` ([rmstconddelta()]).
#' @param s,w Landmarks and window widths for the window and conditional
#'   estimands (see [rmstwindow()]); ignored for `"cumulative"`.
#'
#' @return A list with point estimate and percentile confidence bands, and
#'   `time` (cumulative) or `windows` (window, conditional).
#' @examples
#' time <- c(0, 1, 2)
#' xA <- assurvmat(rbind(c(1.0, 0.9, 0.7), c(1.0, 0.8, 0.6)), time, group = c("A", "A"))
#' xB <- assurvmat(rbind(c(1.0, 0.95, 0.8), c(1.0, 0.85, 0.7)), time, group = c("B", "B"))
#' boot <- bootrmstdelta(xA, xB, R = 10, seed = 1)
#' boot$estimate
#' bootrmstdelta(xA, xB, R = 10, seed = 1, paired = TRUE)$lo
#' bootrmstdelta(xA, xB, R = 10, seed = 1, estimand = "window", s = c(0, 1), w = 1)$estimate
#' @export
bootrmstdelta <- function(xA, xB, R = 300, conf = 0.95, seed = NULL, paired = FALSE,
                          estimand = c("cumulative", "window", "conditional"), s = NULL, w = NULL) {
  estimand <- match.arg(estimand)
  xA <- tosm(xA); xB <- tosm(xB)
  stopifnot(inherits(xA, "survmat"), inherits(xB, "survmat"))
  if (paired && nobssurvmat(xA) != nobssurvmat(xB)) {
    stop("`paired = TRUE` needs the same number of rows in `xA` and `xB`.", call. = FALSE)
  }
  time <- xA$time
  dfun <- switch(estimand,
    cumulative = function(a, b) rmstdelta(a, b)$delta,
    window = function(a, b) rmstwindowdelta(a, b, s, w)$delta,
    conditional = function(a, b) rmstconddelta(a, b, s, w)$delta)

  est_fun <- function(r) {
    if (is.null(r)) {
      dfun(xA, xB)
    } else {
      iA <- resample_idx(nobssurvmat(xA), xA$group)
      iB <- if (paired) iA else resample_idx(nobssurvmat(xB), xB$group)
      xA_b <- assurvmat(xA$S[iA, , drop = FALSE], time, group = xA$group[iA])
      xB_b <- assurvmat(xB$S[iB, , drop = FALSE], time, group = xB$group[iB])
      dfun(xA_b, xB_b)
    }
  }

  out <- bootstrapcurve(est_fun, R = R, conf = conf, seed = seed)
  if (estimand == "cumulative") out$time <- time else out$windows <- wingrid(time, s, w)
  out$estimand <- "bootrmstdelta"
  out
}
