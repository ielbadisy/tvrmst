#' Plot two survival curves (first column only) as ggplot
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
