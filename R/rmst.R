#' RMST at a single horizon tau
#'
#' RMST(tau) = \eqn{\int_0^\tau S(u)\,du} computed by trapezoid rule on the grid.
#'
#' @param S Survival matrix (n_units x n_time).
#' @param time Numeric time grid (length n_time).
#' @param tau Nonnegative scalar horizon.
#' @return Numeric vector of length nrow(S).
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' rmst_tau(S, t, tau = 3)
#' @export
rmst_tau <- function(S, time, tau) {
  S <- .coerce_unit_time(S, time, "S")
  if (!is.numeric(tau) || length(tau) != 1 || !is.finite(tau) || tau < 0) {
    .stop("`tau` must be a single nonnegative number.")
  }
  time <- as.numeric(time)
  ez <- .extend_to_zero(time, S); time <- ez$t; S <- ez$S
  dt <- diff(time)
  S_left <- S[, -ncol(S), drop = FALSE]
  S_right <- S[, -1, drop = FALSE]
  S_mid <- (S_left + S_right) / 2
  area_incr <- sweep(S_mid, 2, dt, "*")
  rmst_no0 <- t(apply(area_incr, 1, cumsum))
  rmst_mat <- cbind(0, rmst_no0)

  if (tau <= time[1]) {
    out <- rep(tau, nrow(S))
    if (!is.null(rownames(S))) names(out) <- rownames(S)
    return(out)
  }

  if (tau >= time[length(time)]) {
    tau_idx <- length(time)
  } else {
    tau_idx <- max(which(time <= tau))
  }

  out <- rmst_mat[, tau_idx]
  if (tau > time[tau_idx]) {
    dt_last <- tau - time[tau_idx]
    out <- out + S[, tau_idx] * dt_last
  }
  if (!is.null(rownames(S))) {
    names(out) <- rownames(S)
  } else {
    names(out) <- paste0("unit", seq_len(nrow(S)))
  }
  out
}

#' RMST curve evaluated at tau = `time[-1]`
#'
#' Returns a data.frame with tau and RMST(tau) for each unit column.
#'
#' @param S Survival matrix (n_units x n_time).
#' @param time Numeric time grid (length n_time).
#' @return data.frame with column tau and one column per unit.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' rmst_curve(S, t)
#' @export
rmst_curve <- function(S, time) {
  S <- .coerce_unit_time(S, time, "S")
  time <- as.numeric(time)
  ez <- .extend_to_zero(time, S); time <- ez$t; S <- ez$S
  dt <- diff(time)
  S_left <- S[, -ncol(S), drop = FALSE]
  S_right <- S[, -1, drop = FALSE]
  S_mid <- (S_left + S_right) / 2
  area_incr <- sweep(S_mid, 2, dt, "*")
  rmst_mat <- t(apply(area_incr, 1, cumsum))
  out <- data.frame(tau = time[-1], t(rmst_mat), check.names = FALSE)
  if (!is.null(rownames(S))) {
    names(out)[-1] <- rownames(S)
  } else {
    names(out)[-1] <- paste0("unit", seq_len(nrow(S)))
  }
  out
}

#' Individual (dynamic) RMST curves aligned to the time grid
#'
#' Returns RMST values for each unit at every grid time, including tau = 0.
#'
#' @param S Survival matrix (n_units x n_time).
#' @param time Numeric time grid (length n_time).
#' @return data.frame with column tau (including 0) and one column per unit.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- rbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' rmst_dynamic(S, t)
#' @export
rmst_dynamic <- function(S, time) {
  S <- .coerce_unit_time(S, time, "S")
  time <- as.numeric(time)
  ez <- .extend_to_zero(time, S); time <- ez$t; S <- ez$S
  dt <- diff(time)
  S_left <- S[, -ncol(S), drop = FALSE]
  S_right <- S[, -1, drop = FALSE]
  S_mid <- (S_left + S_right) / 2
  area_incr <- sweep(S_mid, 2, dt, "*")
  rmst_no0 <- t(apply(area_incr, 1, cumsum))
  rmst_mat <- cbind(0, rmst_no0)
  out <- data.frame(tau = time, t(rmst_mat), check.names = FALSE)
  if (!is.null(rownames(S))) {
    names(out)[-1] <- rownames(S)
  } else {
    names(out)[-1] <- paste0("unit", seq_len(nrow(S)))
  }
  out
}

#' Delta RMST curve between two arms in an rmst_curve() output
#'
#' @param rmst_curve_df Output of rmst_curve().
#' @param arm1 Column name for arm 1.
#' @param arm0 Column name for arm 0.
#' @return data.frame(tau, delta)
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- exp(-0.2 * t)
#' S1 <- exp(-0.15 * t)
#' S <- rbind(Control = S0, Treatment = S1)
#' rc <- rmst_curve(S, t)
#' rmst_delta_curve(rc, "Treatment", "Control")
#' @export
rmst_delta_curve <- function(rmst_curve_df, arm1, arm0) {
  if (!is.data.frame(rmst_curve_df) || !"tau" %in% names(rmst_curve_df)) {
    .stop("`rmst_curve_df` must come from rmst_curve().")
  }
  if (!arm1 %in% names(rmst_curve_df) || !arm0 %in% names(rmst_curve_df)) {
    .stop("`arm1` and `arm0` must be columns in rmst_curve_df.")
  }
  data.frame(tau = rmst_curve_df$tau, delta = rmst_curve_df[[arm1]] - rmst_curve_df[[arm0]])
}
