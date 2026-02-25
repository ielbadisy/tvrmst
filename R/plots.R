#' Plot survival curves from unit x time matrices
#'
#' @param S0 Control survival matrix/data.frame (rows = units, cols = time).
#' @param S1 Treatment survival matrix/data.frame (rows = units, cols = time).
#' @param time Numeric time vector.
#' @param labels Legend labels.
#' @param title Optional title.
#' @param show `"mean"` uses row means; `"first"` uses first row.
#' @return A ggplot object.
#' @export
plot_survival_curves <- function(S0, S1, time,
                                 labels = c("Control", "Treatment"),
                                 title = NULL,
                                 show = c("mean", "first")) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  show <- match.arg(show)

  if (!is.matrix(S0) && !is.data.frame(S0)) .stop("`S0` must be a matrix or data.frame.")
  if (!is.matrix(S1) && !is.data.frame(S1)) .stop("`S1` must be a matrix or data.frame.")
  S0 <- as.matrix(S0)
  S1 <- as.matrix(S1)
  if (!is.numeric(S0) || !is.numeric(S1)) .stop("`S0` and `S1` must be numeric.")
  if (!is.numeric(time) || length(time) < 2 || any(diff(time) <= 0)) {
    .stop("`time` must be numeric, strictly increasing, and length >= 2.")
  }
  if (ncol(S0) != length(time) || ncol(S1) != length(time)) {
    .stop("`S0` and `S1` must have ncol(.) == length(time).")
  }

  y0 <- if (show == "mean") colMeans(S0, na.rm = TRUE) else S0[1, ]
  y1 <- if (show == "mean") colMeans(S1, na.rm = TRUE) else S1[1, ]

  df <- data.frame(
    time = rep(time, 2),
    estimate = c(y0, y1),
    arm = rep(labels, each = length(time))
  )

  ggplot2::ggplot(df, ggplot2::aes(time, estimate, color = arm)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Time", y = "Survival")
}

#' Plot RMST curve object
#'
#' @param rmst_curve_df Data frame with `tau`, `estimate`, optional `lower`,
#'   `upper`, and optional group column.
#' @param title Optional title.
#' @param group_col Optional grouping column name.
#' @return A ggplot object.
#' @export
plot_rmst_curve <- function(rmst_curve_df, title = NULL, group_col = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  if (!is.data.frame(rmst_curve_df) || !all(c("tau", "estimate") %in% names(rmst_curve_df))) {
    .stop("`rmst_curve_df` must have columns `tau` and `estimate`.")
  }

  plot_df <- rmst_curve_df
  if (!is.null(group_col)) {
    if (!group_col %in% names(plot_df)) .stop("`group_col` is not a column in `rmst_curve_df`.")
    plot_df$group <- as.character(plot_df[[group_col]])
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = tau, y = estimate, color = group, fill = group, group = group))
  } else {
    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = tau, y = estimate))
  }

  if (all(c("lower", "upper") %in% names(rmst_curve_df))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      color = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "RMST")
}

#' Plot RMST delta curve
#'
#' @param rmst_delta_df Data frame with `tau`, `estimate`, and optional
#'   `lower`, `upper`.
#' @param title Optional title.
#' @return A ggplot object.
#' @export
plot_rmst_delta <- function(rmst_delta_df, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  if (!is.data.frame(rmst_delta_df) || !all(c("tau", "estimate") %in% names(rmst_delta_df))) {
    .stop("`rmst_delta_df` must have columns `tau` and `estimate`.")
  }

  p <- ggplot2::ggplot(rmst_delta_df, ggplot2::aes(tau, estimate))

  if (all(c("lower", "upper") %in% names(rmst_delta_df))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "Delta RMST")
}

#' Plot conditional tvRMST curves for two arms
#'
#' @param mu0_df Data frame with columns `s`, `estimate`, optional `lower`,
#'   `upper` for control arm.
#' @param mu1_df Data frame with columns `s`, `estimate`, optional `lower`,
#'   `upper` for treatment arm.
#' @param labels Legend labels.
#' @param title Optional title.
#' @return A ggplot object.
#' @export
plot_tvrmst <- function(mu0_df, mu1_df,
                        labels = c("Control", "Treatment"),
                        title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  req <- c("s", "estimate")
  if (!is.data.frame(mu0_df) || !all(req %in% names(mu0_df))) {
    .stop("`mu0_df` must have columns `s` and `estimate`.")
  }
  if (!is.data.frame(mu1_df) || !all(req %in% names(mu1_df))) {
    .stop("`mu1_df` must have columns `s` and `estimate`.")
  }

  d0 <- mu0_df
  d1 <- mu1_df
  d0$arm <- labels[1]
  d1$arm <- labels[2]
  df <- rbind(d0, d1)

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = s, y = estimate, color = arm, fill = arm)
  )

  if (all(c("lower", "upper") %in% names(df))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2,
      color = NA
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Landmark s", y = "Conditional tvRMST")
}

#' Plot tvRMST difference (and bootstrap CI when available)
#'
#' @param delta_df Data frame with columns `s`, `estimate`, optional `lower`,
#'   `upper`.
#' @param title Optional title.
#' @return A ggplot object.
#' @export
plot_tvrmst_diff <- function(delta_df, title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  if (!is.data.frame(delta_df) || !all(c("s", "estimate") %in% names(delta_df))) {
    .stop("`delta_df` must have columns `s` and `estimate`.")
  }

  p <- ggplot2::ggplot(delta_df, ggplot2::aes(s, estimate))

  if (all(c("lower", "upper") %in% names(delta_df))) {
    p <- p + ggplot2::geom_ribbon(
      ggplot2::aes(ymin = lower, ymax = upper),
      alpha = 0.2
    )
  }

  p +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Landmark s", y = "Delta conditional tvRMST")
}

#' Plot individual RMST trajectories
#'
#' @param RMST_mat Numeric matrix from `rmst_dynamic()` (n_units x n_time).
#' @param group Optional grouping vector of length n_units.
#' @param max_units Maximum number of units to plot.
#' @param alpha Line alpha.
#' @param title Optional title.
#' @return A ggplot object.
#' @export
plot_rmst_individual <- function(RMST_mat, group = NULL,
                                 max_units = 150,
                                 alpha = 0.1,
                                 title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  if (!is.matrix(RMST_mat) || !is.numeric(RMST_mat)) {
    .stop("`RMST_mat` must be a numeric matrix from rmst_dynamic().")
  }

  n_units <- nrow(RMST_mat)
  if (!is.null(group) && length(group) != n_units) {
    .stop("`group` must be NULL or length nrow(RMST_mat).")
  }

  ids <- seq_len(n_units)
  if (n_units > max_units) {
    set.seed(1)
    ids <- sort(sample(ids, max_units))
  }

  tau <- suppressWarnings(as.numeric(colnames(RMST_mat)))
  if (any(!is.finite(tau))) tau <- seq_len(ncol(RMST_mat))

  df <- data.frame(
    id = rep(ids, each = length(tau)),
    tau = rep(tau, times = length(ids)),
    rmst = as.vector(t(RMST_mat[ids, , drop = FALSE]))
  )

  if (!is.null(group)) {
    df$group <- as.character(group[df$id])
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(tau, rmst, group = id)) +
    ggplot2::geom_line(alpha = alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "RMST")

  if (!is.null(group)) {
    p <- p + ggplot2::facet_wrap(~group)
  }

  p
}

#' Plot group-level RMST summaries
#'
#' @param RMST_mat Numeric matrix from `rmst_dynamic()` (n_units x n_time).
#' @param group Group labels of length n_units.
#' @param statistic `"mean"` or `"median"`.
#' @param title Optional title.
#' @return A ggplot object.
#' @export
plot_rmst_mean <- function(RMST_mat, group,
                           statistic = c("mean", "median"),
                           title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    .stop("ggplot2 is required for plotting.")
  }
  statistic <- match.arg(statistic)

  if (!is.matrix(RMST_mat) || !is.numeric(RMST_mat)) {
    .stop("`RMST_mat` must be a numeric matrix from rmst_dynamic().")
  }
  if (missing(group) || length(group) != nrow(RMST_mat)) {
    .stop("`group` must have length nrow(RMST_mat).")
  }

  tau <- suppressWarnings(as.numeric(colnames(RMST_mat)))
  if (any(!is.finite(tau))) tau <- seq_len(ncol(RMST_mat))

  g <- as.character(group)
  groups <- unique(g)

  df <- do.call(
    rbind,
    lapply(groups, function(gr) {
      idx <- which(g == gr)
      vals <- if (statistic == "mean") {
        colMeans(RMST_mat[idx, , drop = FALSE], na.rm = TRUE)
      } else {
        apply(RMST_mat[idx, , drop = FALSE], 2, stats::median, na.rm = TRUE)
      }
      data.frame(group = gr, tau = tau, estimate = as.numeric(vals))
    })
  )

  ggplot2::ggplot(df, ggplot2::aes(tau, estimate, color = group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = title, x = "Horizon tau", y = "RMST")
}
