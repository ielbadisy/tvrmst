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
#' boot <- bootstrap_curve(
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
bootstrap_curve <- function(est_fun, R = 300, conf = 0.95, seed = NULL, keep_reps = FALSE) {
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
#' @param xA A `survmat` object for arm A.
#' @param xB A `survmat` object for arm B.
#' @param R Number of bootstrap replicates.
#' @param conf Confidence level in `(0, 1)`.
#' @param seed Optional RNG seed.
#'
#' @return A list with point estimate and percentile confidence bands.
#' @examples
#' time <- c(0, 1, 2)
#' xA <- as_survmat(rbind(c(1.0, 0.9, 0.7), c(1.0, 0.8, 0.6)), time, group = c("A", "A"))
#' xB <- as_survmat(rbind(c(1.0, 0.95, 0.8), c(1.0, 0.85, 0.7)), time, group = c("B", "B"))
#' boot <- boot_rmst_delta(xA, xB, R = 10, seed = 1)
#' boot$estimate
#' @export
boot_rmst_delta <- function(xA, xB, R = 300, conf = 0.95, seed = NULL) {
  stopifnot(inherits(xA, "survmat"), inherits(xB, "survmat"))
  time <- xA$time

  est_fun <- function(r) {
    if (is.null(r)) {
      rmst_delta(xA, xB)$delta
    } else {
      iA <- resample_idx(nobs_survmat(xA), xA$group)
      iB <- resample_idx(nobs_survmat(xB), xB$group)
      xA_b <- as_survmat(xA$S[iA, , drop = FALSE], time, group = xA$group[iA])
      xB_b <- as_survmat(xB$S[iB, , drop = FALSE], time, group = xB$group[iB])
      rmst_delta(xA_b, xB_b)$delta
    }
  }

  out <- bootstrap_curve(est_fun, R = R, conf = conf, seed = seed)
  out$time <- time
  out$estimand <- "boot_rmst_delta"
  out
}
