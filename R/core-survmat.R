#' Construct a survival matrix object
#'
#' @param S Numeric matrix of survival probabilities (rows = subjects, columns = time).
#' @param time Numeric strictly increasing time grid with `length(time) == ncol(S)`.
#' @param id Optional subject identifiers. Defaults to `seq_len(nrow(S))`.
#' @param group Optional group labels with length `nrow(S)`.
#'
#' @return An object of class `"survmat"`.
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
  stopifnot(all(is.finite(S)))
  stopifnot(all(S >= 0 & S <= 1))

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
