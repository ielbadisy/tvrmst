#' Coerce survival predictions to unit-first orientation
#'
#' @param S Survival matrix or data.frame.
#' @param time Numeric time grid (length n_time).
#' @param by One of "cols", "rows", or "auto".
#' @return Numeric matrix with rows = units and columns = time.
#' @examples
#' time <- seq(0, 5, by = 1)
#' S_time_rows <- cbind(exp(-0.2 * time), exp(-0.3 * time))
#' S_units <- as_unit_time(S_time_rows, time, by = "rows")
#' rmst_dynamic(S_units, time)
#' @export
as_unit_time <- function(S, time, by = c("cols", "rows", "auto")) {
  by <- match.arg(by)
  .check_time(time)

  if (is.data.frame(S)) S <- as.matrix(S)
  if (!is.matrix(S)) .stop("`S` must be a numeric matrix or data.frame.")
  if (!is.numeric(S)) .stop("`S` must be numeric.")

  n_rows <- nrow(S)
  n_cols <- ncol(S)
  t_len <- length(time)

  if (by == "cols") {
    if (n_cols != t_len) {
      .stop("`S` must be n_units x n_time (columns correspond to `time`). Ensure ncol(S) == length(time).")
    }
    return(S)
  }

  if (by == "rows") {
    if (n_rows != t_len) {
      .stop("`S` must be n_time x n_units when by = \"rows\" (nrow(S) == length(time)).")
    }
    return(t(S))
  }

  match_rows <- n_rows == t_len
  match_cols <- n_cols == t_len

  if (match_cols && !match_rows) return(S)
  if (match_rows && !match_cols) return(t(S))
  if (match_rows && match_cols) {
    .stop("Ambiguous orientation: both nrow(S) and ncol(S) match length(time). Use by = \"rows\" or \"cols\".")
  }

  .stop("`S` must be n_units x n_time (columns correspond to `time`). Ensure ncol(S) == length(time).")
}
