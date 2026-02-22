#' Plot two survival curves (first column only) as ggplot
#'
#' @param t Time grid.
#' @param S0 Survival matrix for control arm.
#' @param S1 Survival matrix for treatment arm.
#' @param labels Character vector of length 2 for legend labels.
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- cbind(A = exp(-0.2 * t))
#' S1 <- cbind(A = exp(-0.15 * t))
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_survival_curves(t, S0, S1)
#' }
#' @export
plot_survival_curves <- function(t, S0, S1, labels = c("Control","Treatment"),
                                 title = "Survival curves") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  .check_survmat(t, as.matrix(S0), "S0")
  .check_survmat(t, as.matrix(S1), "S1")

  df <- data.frame(
    time = rep(t, 2),
    survival = c(S0[, 1], S1[, 1]),
    arm = rep(labels, each = length(t))
  )
  ggplot2::ggplot(df, ggplot2::aes(time, survival, color = arm)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Time", y = "Survival")
}

#' Plot conditional tvRMST curves for two arms (first series by default)
#'
#' @param mu0_df Output of `tvrmst_cond()` for arm 0.
#' @param mu1_df Output of `tvrmst_cond()` for arm 1.
#' @param labels Character vector of length 2 for legend labels.
#' @param series_col Optional column name to plot.
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' S1 <- cbind(A = exp(-0.15 * t), B = exp(-0.25 * t))
#' s_grid <- c(0, 1, 2)
#' mu0 <- tvrmst_cond(t, S0, s_grid, tau = 2)
#' mu1 <- tvrmst_cond(t, S1, s_grid, tau = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_tvrmst(mu0, mu1, series_col = "A")
#' }
#' @export
plot_tvrmst <- function(mu0_df, mu1_df, labels = c("Control","Treatment"),
                        series_col = NULL,
                        title = "Conditional tvRMST mu_c(s,tau) by arm") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  if (!is.data.frame(mu0_df) || !is.data.frame(mu1_df) || !"s" %in% names(mu0_df) || !"s" %in% names(mu1_df)) {
    .stop("mu0_df and mu1_df must be outputs of tvrmst_cond() (with column 's').")
  }

  cols0 <- setdiff(names(mu0_df), "s")
  cols1 <- setdiff(names(mu1_df), "s")
  if (length(cols0) < 1 || length(cols1) < 1) .stop("No series columns in tvrmst_cond outputs.")

  if (is.null(series_col)) {
    c0 <- cols0[1]; c1 <- cols1[1]
  } else {
    if (!series_col %in% cols0 || !series_col %in% cols1) .stop("series_col must exist in both inputs.")
    c0 <- series_col; c1 <- series_col
  }

  df <- rbind(
    data.frame(s = mu0_df$s, value = mu0_df[[c0]], arm = labels[1]),
    data.frame(s = mu1_df$s, value = mu1_df[[c1]], arm = labels[2])
  )

  ggplot2::ggplot(df, ggplot2::aes(s, value, color = arm)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Landmark time s", y = "mu_c(s,tau)")
}

#' Plot Delta_c(s,tau) with optional CI (estimate/lower/upper)
#'
#' @param delta_df Output of `tvrmst_diff()` or a data.frame with columns
#'   `s`, `estimate`, `lower`, `upper`.
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' S1 <- cbind(A = exp(-0.15 * t), B = exp(-0.25 * t))
#' s_grid <- c(0, 1, 2)
#' delta <- tvrmst_diff(t, S1, S0, s_grid, tau = 2)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_tvrmst_diff(delta)
#' }
#' @export
plot_tvrmst_diff <- function(delta_df, title = "Delta_c(s,tau)") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  if (!is.data.frame(delta_df) || !"s" %in% names(delta_df)) .stop("delta_df must have column 's'.")

  if (all(c("estimate","lower","upper") %in% names(delta_df))) {
    ggplot2::ggplot(delta_df, ggplot2::aes(s, estimate)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), alpha = 0.2) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = title, x = "Landmark time s", y = "Delta_c(s,tau)")
  } else {
    ycol <- if ("delta_c" %in% names(delta_df)) "delta_c" else setdiff(names(delta_df), "s")[1]
    ggplot2::ggplot(delta_df, ggplot2::aes(s, .data[[ycol]])) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::theme_minimal() +
      ggplot2::labs(title = title, x = "Landmark time s", y = "Delta_c(s,tau)")
  }
}

#' Plot RMST curves over tau for multiple series
#'
#' @param rmst_curve_df Output of `rmst_curve()`.
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t))
#' rc <- rmst_curve(t, S)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_rmst_curve(rc)
#' }
#' @export
plot_rmst_curve <- function(rmst_curve_df, title = "RMST curve") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  if (!is.data.frame(rmst_curve_df) || !"tau" %in% names(rmst_curve_df)) {
    .stop("`rmst_curve_df` must come from rmst_curve().")
  }
  series <- setdiff(names(rmst_curve_df), "tau")
  if (length(series) < 1) .stop("No series columns in rmst_curve_df.")

  df <- do.call(
    rbind,
    lapply(series, function(s) {
      data.frame(tau = rmst_curve_df$tau, value = rmst_curve_df[[s]], series = s)
    })
  )

  ggplot2::ggplot(df, ggplot2::aes(.data$tau, .data$value, color = .data$series)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "RMST(tau)")
}

#' Plot Delta RMST curve over tau
#'
#' @param rmst_delta_df Output of `rmst_delta_curve()`.
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S0 <- exp(-0.2 * t)
#' S1 <- exp(-0.15 * t)
#' rc <- rmst_curve(t, cbind(Control = S0, Treatment = S1))
#' dr <- rmst_delta_curve(rc, "Treatment", "Control")
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_rmst_delta(dr)
#' }
#' @export
plot_rmst_delta <- function(rmst_delta_df, title = "Delta RMST curve") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  if (!is.data.frame(rmst_delta_df) || !"tau" %in% names(rmst_delta_df) || !"delta" %in% names(rmst_delta_df)) {
    .stop("`rmst_delta_df` must come from rmst_delta_curve().")
  }

  ggplot2::ggplot(rmst_delta_df, ggplot2::aes(.data$tau, .data$delta)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "Delta RMST(tau)")
}

#' Plot individual dynamic RMST curves with optional grouping
#'
#' @param rmst_dynamic_df Output of `rmst_dynamic()`.
#' @param group Optional vector of group labels (length = number of series).
#' @param alpha Line transparency for individual curves.
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t), C = exp(-0.25 * t))
#' rd <- rmst_dynamic(t, S)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_rmst_individual(rd, group = c("G1", "G2", "G1"))
#' }
#' @export
plot_rmst_individual <- function(rmst_dynamic_df, group = NULL, alpha = 0.15,
                                 title = "Individual RMST curves") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  if (!is.data.frame(rmst_dynamic_df) || !"tau" %in% names(rmst_dynamic_df)) {
    .stop("`rmst_dynamic_df` must come from rmst_dynamic().")
  }
  series <- setdiff(names(rmst_dynamic_df), "tau")
  if (length(series) < 1) .stop("No series columns in rmst_dynamic_df.")
  if (!is.null(group) && length(group) != length(series)) {
    .stop("`group` must be NULL or have length equal to number of series.")
  }

  df <- do.call(
    rbind,
    lapply(seq_along(series), function(i) {
      data.frame(
        tau = rmst_dynamic_df$tau,
        rmst = rmst_dynamic_df[[series[i]]],
        id = series[i],
        group = if (is.null(group)) "All" else as.character(group[i])
      )
    })
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$tau, .data$rmst, group = .data$id)) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "RMST(tau)")

  if (!is.null(group)) p <- p + ggplot2::facet_wrap(~ group)
  p
}

#' Plot mean dynamic RMST curves by group
#'
#' @param rmst_dynamic_df Output of `rmst_dynamic()`.
#' @param group Group labels (length = number of series).
#' @param title Plot title.
#' @examples
#' t <- seq(0, 5, by = 1)
#' S <- cbind(A = exp(-0.2 * t), B = exp(-0.3 * t), C = exp(-0.25 * t))
#' rd <- rmst_dynamic(t, S)
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plot_rmst_mean(rd, group = c("G1", "G2", "G1"))
#' }
#' @export
plot_rmst_mean <- function(rmst_dynamic_df, group, title = "Mean RMST by group") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting. Install it or avoid plot_*().")
  }
  if (!is.data.frame(rmst_dynamic_df) || !"tau" %in% names(rmst_dynamic_df)) {
    .stop("`rmst_dynamic_df` must come from rmst_dynamic().")
  }
  series <- setdiff(names(rmst_dynamic_df), "tau")
  if (length(series) < 1) .stop("No series columns in rmst_dynamic_df.")
  if (missing(group) || length(group) != length(series)) {
    .stop("`group` must be provided and have length equal to number of series.")
  }

  df <- do.call(
    rbind,
    lapply(seq_along(series), function(i) {
      data.frame(
        tau = rmst_dynamic_df$tau,
        rmst = rmst_dynamic_df[[series[i]]],
        group = as.character(group[i])
      )
    })
  )

  ggplot2::ggplot(df, ggplot2::aes(.data$tau, .data$rmst, group = .data$group)) +
    ggplot2::stat_summary(fun = mean, geom = "line") +
    ggplot2::facet_wrap(~ group) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "RMST(tau)")
}
