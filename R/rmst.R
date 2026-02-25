#' Dynamic RMST on the provided time grid
#'
#' @param S Survival predictions as matrix/data.frame with rows = units and
#'   columns = time points.
#' @param time Numeric strictly increasing time vector.
#' @return Numeric matrix of shape n_units x length(time).
#' @export
rmst_dynamic <- function(S, time) {
  S <- validate_surv_input(S, time, name = "S")
  time <- as.numeric(time)

  n_units <- nrow(S)
  n_time <- ncol(S)
  out <- matrix(0, nrow = n_units, ncol = n_time)

  if (n_time > 1) {
    dt <- diff(time)
    mids <- (S[, -n_time, drop = FALSE] + S[, -1, drop = FALSE]) / 2
    incr <- sweep(mids, 2, dt, "*")
    out[, -1] <- t(apply(incr, 1, cumsum))
  }

  if (!is.null(rownames(S))) rownames(out) <- rownames(S)
  colnames(out) <- as.character(time)
  out
}

#' RMST at a fixed horizon
#'
#' @param S Survival predictions as matrix/data.frame with rows = units and
#'   columns = time points.
#' @param time Numeric strictly increasing time vector.
#' @param tau Numeric scalar horizon.
#' @return Numeric vector length n_units. Attribute `tau_used` records the grid
#'   time used after nearest-grid snapping.
#' @export
rmst_tau <- function(S, time, tau) {
  S <- validate_surv_input(S, time, name = "S")
  time <- as.numeric(time)

  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau)) {
    .stop("`tau` must be a single finite numeric value.")
  }

  rmst_mat <- rmst_dynamic(S, time)
  idx <- which.min(abs(time - tau))
  out <- rmst_mat[, idx]

  if (!is.null(rownames(S))) {
    names(out) <- rownames(S)
  } else {
    names(out) <- paste0("unit", seq_len(nrow(S)))
  }
  attr(out, "tau_used") <- time[idx]
  out
}

#' RMST curve summarized across units
#'
#' @param S Survival predictions as matrix/data.frame with rows = units and
#'   columns = time points.
#' @param time Numeric strictly increasing time vector.
#' @param statistic Summary statistic across units: `"mean"` or `"median"`.
#' @param probs Quantiles for interval summary when `conf = NULL`.
#' @param conf Optional confidence level in (0,1) for normal-approx intervals.
#' @return data.frame with columns `tau`, `estimate`, and optional `lower`,
#'   `upper`.
#' @export
rmst_curve <- function(S, time,
                       statistic = c("mean", "median"),
                       probs = c(0.025, 0.975),
                       conf = NULL) {
  S <- validate_surv_input(S, time, name = "S")
  statistic <- match.arg(statistic)
  time <- as.numeric(time)

  rmst_mat <- rmst_dynamic(S, time)

  if (statistic == "mean") {
    estimate <- colMeans(rmst_mat, na.rm = TRUE)
  } else {
    estimate <- apply(rmst_mat, 2, stats::median, na.rm = TRUE)
  }

  out <- data.frame(tau = time, estimate = estimate)

  if (!is.null(conf)) {
    if (!is.numeric(conf) || length(conf) != 1 || !is.finite(conf) || conf <= 0 || conf >= 1) {
      .stop("`conf` must be a single number in (0,1).")
    }
    z <- stats::qnorm(1 - (1 - conf) / 2)
    if (statistic == "mean") {
      se <- apply(rmst_mat, 2, stats::sd, na.rm = TRUE) / sqrt(nrow(rmst_mat))
    } else {
      se <- apply(rmst_mat, 2, stats::mad, constant = 1.4826, na.rm = TRUE) / sqrt(nrow(rmst_mat))
    }
    out$lower <- estimate - z * se
    out$upper <- estimate + z * se
    return(out)
  }

  if (!is.null(probs)) {
    if (!is.numeric(probs) || length(probs) != 2 || any(!is.finite(probs)) ||
        any(probs < 0) || any(probs > 1) || probs[1] >= probs[2]) {
      .stop("`probs` must be two increasing probabilities in [0,1].")
    }
    qmat <- apply(
      rmst_mat,
      2,
      stats::quantile,
      probs = probs,
      na.rm = TRUE,
      names = FALSE
    )
    out$lower <- as.numeric(qmat[1, ])
    out$upper <- as.numeric(qmat[2, ])
  }

  out
}

#' Delta RMST curve from a wide rmst-curve table
#'
#' @param rmst_curve_df Wide data.frame with `tau` and one column per arm.
#' @param arm1 Column name for arm 1.
#' @param arm0 Column name for arm 0.
#' @return data.frame with columns `tau`, `estimate`.
#' @export
rmst_delta_curve <- function(rmst_curve_df, arm1, arm0) {
  if (!is.data.frame(rmst_curve_df) || !"tau" %in% names(rmst_curve_df)) {
    .stop("`rmst_curve_df` must be a data.frame with column `tau`.")
  }
  if (!arm1 %in% names(rmst_curve_df) || !arm0 %in% names(rmst_curve_df)) {
    .stop("`arm1` and `arm0` must be columns in `rmst_curve_df`.")
  }
  data.frame(
    tau = rmst_curve_df$tau,
    estimate = rmst_curve_df[[arm1]] - rmst_curve_df[[arm0]]
  )
}
