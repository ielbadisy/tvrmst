#' Construct a survival matrix object
#'
#' The time grid must start at 0 with `S = 1` for every subject, since RMST is
#' integrated from 0. If the grid starts after 0, a column `t = 0`, `S = 1` is
#' prepended automatically and a message is emitted. If the grid starts at 0
#' and some `S[, 1]` differ from 1, an error is raised: the input is not a
#' valid survival curve and is not silently altered.
#'
#' Each row of `S` must be a valid survival curve: values in `[0, 1]` and
#' non-increasing over time. Violations larger than a numerical tolerance of
#' `sqrt(.Machine$double.eps)` are errors.
#'
#' @param S Numeric matrix of survival probabilities (rows = subjects, columns = time).
#' @param time Numeric strictly increasing time grid with `length(time) == ncol(S)`.
#'   Must be non-negative; see Details on the grid at 0.
#' @param id Optional subject identifiers. Defaults to `seq_len(nrow(S))`.
#' @param group Optional group labels with length `nrow(S)`.
#'
#' @return An object of class `"survmat"`. Its time grid always starts at 0.
#' @examples
#' time <- c(0, 1, 2)
#' S <- rbind(
#'   c(1.0, 0.8, 0.6),
#'   c(1.0, 0.7, 0.5)
#' )
#' x <- as_survmat(S, time, group = c("A", "B"))
#' x
#' @export
as_survmat <- function(S, time, id = NULL, group = NULL) {
  stopifnot(is.matrix(S))
  stopifnot(is.numeric(time), length(time) == ncol(S))
  stopifnot(all(is.finite(time)), all(diff(time) > 0))
  if (!all(is.finite(S))) stop("`S` must contain only finite values.", call. = FALSE)
  if (any(time < 0)) stop("`time` must be non-negative.", call. = FALSE)

  tol <- sqrt(.Machine$double.eps)
  if (any(S < -tol | S > 1 + tol)) {
    stop("`S` must lie in [0, 1].", call. = FALSE)
  }
  S <- pmin(pmax(S, 0), 1)

  if (ncol(S) > 1) {
    incr <- S[, -1, drop = FALSE] - S[, -ncol(S), drop = FALSE]
    bad <- which(apply(incr > tol, 1, any))
    if (length(bad)) {
      stop(sprintf(
        "`S` must be non-increasing over time; %d row(s) increase (first: row %d).",
        length(bad), bad[1]
      ), call. = FALSE)
    }
  }

  if (time[1] > 0) {
    message("Time grid does not start at 0: prepending t = 0 with S = 1.")
    time <- c(0, time)
    S <- cbind(1, S)
  } else {
    bad0 <- which(abs(S[, 1] - 1) > tol)
    if (length(bad0)) {
      stop(sprintf(
        "S(0) must equal 1 at t = 0; %d row(s) differ (first: row %d, S = %g).",
        length(bad0), bad0[1], S[bad0[1], 1]
      ), call. = FALSE)
    }
    S[, 1] <- 1
  }

  n <- nrow(S)

  if (is.null(id)) id <- seq_len(n)
  stopifnot(length(id) == n)

  if (!is.null(group)) {
    stopifnot(length(group) == n)
    group <- as.factor(group)
  }

  structure(
    list(S = S, time = as.numeric(time), id = id, group = group),
    class = "survmat"
  )
}

#' Number of rows in a survmat
#'
#' @param x A `survmat` object.
#'
#' @return Number of subjects.
#' @examples
#' time <- c(0, 1, 2)
#' S <- rbind(
#'   c(1.0, 0.8, 0.6),
#'   c(1.0, 0.7, 0.5)
#' )
#' x <- as_survmat(S, time)
#' nobs_survmat(x)
#' @export
nobs_survmat <- function(x) {
  stopifnot(inherits(x, "survmat"))
  nrow(x$S)
}

#' Row-bind survmat objects on a common time grid
#'
#' @param ... One or more `survmat` objects.
#'
#' @return A combined `survmat`.
#' @examples
#' time <- c(0, 1, 2)
#' xA <- as_survmat(matrix(c(1.0, 0.9, 0.7), nrow = 1), time, group = "A")
#' xB <- as_survmat(matrix(c(1.0, 0.8, 0.6), nrow = 1), time, group = "B")
#' x <- bind_survmat(xA, xB)
#' nobs_survmat(x)
#' @export
bind_survmat <- function(...) {
  xs <- list(...)
  stopifnot(length(xs) >= 1)
  stopifnot(all(vapply(xs, inherits, logical(1), "survmat")))

  time <- xs[[1]]$time
  ok_time <- vapply(
    xs,
    function(z) length(z$time) == length(time) && all(abs(z$time - time) < 1e-12),
    logical(1)
  )
  stopifnot(all(ok_time))

  S <- do.call(rbind, lapply(xs, `[[`, "S"))
  group <- do.call(c, lapply(xs, function(z) {
    if (!is.null(z$group)) as.character(z$group) else rep("G", nrow(z$S))
  }))

  as_survmat(S, time, id = seq_len(nrow(S)), group = group)
}
