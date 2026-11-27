trapz_cum <- function(time, Y) {
  stopifnot(is.matrix(Y), length(time) == ncol(Y))
  dt <- diff(time)
  mids <- (Y[, -1, drop = FALSE] + Y[, -ncol(Y), drop = FALSE]) / 2
  inc <- mids * matrix(dt, nrow(Y), length(dt), byrow = TRUE)
  cbind(0, t(apply(inc, 1, cumsum)))
}

lininterp_vec <- function(x, y, xout) {
  stats::approx(x = x, y = y, xout = xout, rule = 2, ties = "ordered")$y
}

resample_idx <- function(n, group = NULL) {
  if (is.null(group)) {
    sample.int(n, n, replace = TRUE)
  } else {
    group <- as.factor(group)
    unlist(lapply(levels(group), function(g) {
      idx <- which(group == g)
      sample(idx, length(idx), replace = TRUE)
    }), use.names = FALSE)
  }
}

# Accept a survival matrix carrying its grid in attr "times" (e.g. survalis
# survmat, which shares the class name but is a matrix, not a list).
tosm <- function(x) {
  if (is.matrix(x) && !is.null(attr(x, "times"))) assurvmat(x, attr(x, "times")) else x
}
