#' Coerce survival predictions to a numeric matrix
#'
#' @param pred Prediction output from a model.
#' @param times Unused placeholder for API compatibility.
#'
#' @return Numeric matrix with rows as subjects and columns as times.
#' @examples
#' pred <- data.frame(
#'   t0 = c(1.0, 1.0),
#'   t1 = c(0.8, 0.7),
#'   t2 = c(0.6, 0.5)
#' )
#' as_survprob_matrix(pred)
#' @export
as_survprob_matrix <- function(pred, times = NULL) {
  if (is.matrix(pred)) return(pred)

  if (is.data.frame(pred)) {
    M <- as.matrix(pred)
    storage.mode(M) <- "double"
    return(M)
  }

  if (is.list(pred)) {
    if (!is.null(pred$survival)) {
      M <- as.matrix(pred$survival)
      storage.mode(M) <- "double"
      return(M)
    }
    if (!is.null(pred$pred) && (is.matrix(pred$pred) || is.data.frame(pred$pred))) {
      M <- as.matrix(pred$pred)
      storage.mode(M) <- "double"
      return(M)
    }
  }

  stop("Unsupported predict() output. Try str(pred) to inspect.", call. = FALSE)
}
