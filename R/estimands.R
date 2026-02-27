#' Dynamic RMST trajectories
#'
#' Computes `RMST_i(tau) = integral_0^tau S_i(u) du` over the provided grid.
#'
#' @param x A `survmat` object.
#' @param tau Optional numeric horizons for interpolation.
#' @param by Optional grouping variable. Defaults to `x$group`.
#'
#' @return A list with individual and mean RMST curves.
#' @export
rmst_dynamic <- function(x, tau = NULL, by = x$group) {
  stopifnot(inherits(x, "survmat"))
  time <- x$time
  ind <- trapz_cum(time, x$S)
  mean_curve <- colMeans(ind)

  out <- list(
    estimand = "rmst_dynamic",
    time = time,
    individual = ind,
    mean = mean_curve
  )

  if (!is.null(by)) {
    by <- as.factor(by)
    out$by <- do.call(rbind, lapply(levels(by), function(g) {
      idx <- which(by == g)
      data.frame(group = g, tau = time, estimate = colMeans(ind[idx, , drop = FALSE]))
    }))
  }

  if (!is.null(tau)) {
    tau <- as.numeric(tau)
    out$at_tau <- data.frame(tau = tau, estimate = lininterp_vec(time, mean_curve, tau))

    if (!is.null(by)) {
      out$by_at_tau <- do.call(rbind, lapply(levels(by), function(g) {
        idx <- which(by == g)
        mc <- colMeans(ind[idx, , drop = FALSE])
        data.frame(group = g, tau = tau, estimate = lininterp_vec(time, mc, tau))
      }))
    }
  }

  out
}

#' Dynamic RMST difference between two arms
#'
#' Computes `Delta(tau) = RMST_B(tau) - RMST_A(tau)`.
#'
#' @param xA A `survmat` object for arm A.
#' @param xB A `survmat` object for arm B.
#' @param tau Optional numeric horizons for interpolation.
#'
#' @return A list with RMST curves for both arms and the delta curve.
#' @export
rmst_delta <- function(xA, xB, tau = NULL) {
  stopifnot(inherits(xA, "survmat"), inherits(xB, "survmat"))
  stopifnot(length(xA$time) == length(xB$time), all(abs(xA$time - xB$time) < 1e-12))

  rA <- rmst_dynamic(xA, by = NULL)
  rB <- rmst_dynamic(xB, by = NULL)

  out <- list(
    estimand = "rmst_delta",
    time = rA$time,
    rmst_A = rA$mean,
    rmst_B = rB$mean,
    delta = rB$mean - rA$mean
  )

  if (!is.null(tau)) {
    tau <- as.numeric(tau)
    out$at_tau <- data.frame(
      tau = tau,
      rmst_A = lininterp_vec(out$time, out$rmst_A, tau),
      rmst_B = lininterp_vec(out$time, out$rmst_B, tau)
    )
    out$at_tau$delta <- out$at_tau$rmst_B - out$at_tau$rmst_A
  }

  out
}
