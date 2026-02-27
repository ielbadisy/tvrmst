#' Plot individual dynamic RMST curves by group
#'
#' @param res Result from [rmst_dynamic()] containing `individual` and `time`.
#' @param group Group vector aligned with rows of `res$individual`.
#' @param n_show_per_group Maximum number of individual curves shown per group.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_rmst_individual_by_group <- function(res, group, n_show_per_group = 30,
                                          title = "Individual dynamic RMST by group") {
  .require_ggplot2()
  stopifnot(!is.null(res$individual))

  ind <- res$individual
  stopifnot(nrow(ind) == length(group))
  group <- as.factor(group)

  idx <- unlist(lapply(levels(group), function(g) {
    ids <- which(group == g)
    if (length(ids) <= n_show_per_group) ids else sample(ids, n_show_per_group)
  }), use.names = FALSE)

  ind_sub <- ind[idx, , drop = FALSE]
  g_sub <- group[idx]
  k <- nrow(ind_sub)
  m <- ncol(ind_sub)

  df <- data.frame(
    id = rep(seq_len(k), each = m),
    group = rep(g_sub, each = m),
    tau = rep(res$time, times = k),
    estimate = as.vector(t(ind_sub))
  )

  df_mean <- do.call(rbind, lapply(levels(group), function(g) {
    ids <- which(group == g)
    data.frame(group = g, tau = res$time, mean = colMeans(ind[ids, , drop = FALSE]))
  }))

  ggplot2::ggplot(df, ggplot2::aes(x = tau, y = estimate, group = id)) +
    ggplot2::geom_line(alpha = 0.25) +
    ggplot2::geom_line(
      data = df_mean,
      ggplot2::aes(x = tau, y = mean),
      inherit.aes = FALSE,
      linewidth = 1
    ) +
    ggplot2::facet_wrap(~group) +
    ggplot2::labs(title = title, x = "tau", y = "RMST_i(tau)")
}

#' Plot mean RMST curves for two arms
#'
#' @param xA A `survmat` object for arm A.
#' @param xB A `survmat` object for arm B.
#' @param labels Two legend labels.
#' @param title Plot title.
#'
#' @return A `ggplot` object.
#' @export
plot_rmst_two_arms <- function(xA, xB,
                               labels = c("Arm A", "Arm B"),
                               title = "Dynamic RMST (Two Arms)") {
  .require_ggplot2()
  stopifnot(inherits(xA, "survmat"), inherits(xB, "survmat"))
  stopifnot(length(xA$time) == length(xB$time), all(abs(xA$time - xB$time) < 1e-12))

  rA <- rmst_dynamic(xA, by = NULL)
  rB <- rmst_dynamic(xB, by = NULL)

  df <- rbind(
    data.frame(group = labels[1], tau = rA$time, estimate = rA$mean),
    data.frame(group = labels[2], tau = rB$time, estimate = rB$mean)
  )

  ggplot2::ggplot(df, ggplot2::aes(x = tau, y = estimate, group = group)) +
    ggplot2::geom_line(ggplot2::aes(linetype = group), linewidth = 1) +
    ggplot2::labs(title = title, x = "Time (tau)", y = "RMST(tau)") +
    ggplot2::theme_minimal()
}

#' Plot a delta curve
#'
#' @param grid X-axis grid.
#' @param delta Y values.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#'
#' @return A `ggplot` object.
#' @export
plot_delta_curve <- function(grid, delta, title = "Delta curve", xlab = "t", ylab = "Delta") {
  .require_ggplot2()
  df <- data.frame(t = grid, delta = delta)
  ggplot2::ggplot(df, ggplot2::aes(x = t, y = delta)) +
    ggplot2::geom_line() +
    ggplot2::labs(title = title, x = xlab, y = ylab)
}

#' Plot bootstrap estimate with confidence ribbon
#'
#' @param boot List returned by [bootstrap_curve()] or [boot_rmst_delta()].
#' @param grid Optional x-axis grid. If omitted, uses `boot$time`.
#' @param title Plot title.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#'
#' @return A `ggplot` object.
#' @export
plot_boot_curve <- function(boot, grid = NULL, title = "Bootstrap curve", xlab = "t", ylab = "estimate") {
  .require_ggplot2()
  if (is.null(grid)) {
    if (!is.null(boot$time)) {
      grid <- boot$time
    } else {
      stop("Provide `grid=` or ensure boot has $time.", call. = FALSE)
    }
  }

  df <- data.frame(t = grid, estimate = boot$estimate, lo = boot$lo, hi = boot$hi)
  ggplot2::ggplot(df, ggplot2::aes(x = t, y = estimate)) +
    ggplot2::geom_line() +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), alpha = 0.2) +
    ggplot2::labs(title = title, x = xlab, y = ylab)
}
