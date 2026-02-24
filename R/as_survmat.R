#' Construct a survmat object
#'
#' @param t Numeric time grid (strictly increasing).
#' @param S Numeric matrix (n_time x n_series).
#' @return An object of class "survmat" with components t and S.
#' @export
survmat <- function(t, S) {
  if (!is.numeric(t) || length(t) < 2) {
    .stop("`t` must be numeric, length >= 2.")
  }
  if (any(!is.finite(t))) {
    .stop("`t` contains non-finite values.")
  }
  if (!all(diff(t) > 0)) {
    .stop("`t` must be strictly increasing.")
  }
  if (!is.matrix(S)) {
    .stop("`S` must be a numeric matrix (n_time x n_series).")
  }
  if (!is.numeric(S)) {
    .stop("`S` must be numeric.")
  }
  if (nrow(S) != length(t)) {
    .stop("nrow(S) must equal length(t).")
  }
  if (any(!is.finite(S))) {
    .stop("`S` contains non-finite values.")
  }
  structure(list(t = t, S = S), class = "survmat")
}

#' Coerce inputs to a survmat object
#'
#' @param x Survival predictions: matrix, data.frame, list-like object, or survmat.
#' @param t Optional numeric time grid. If omitted, `as_survmat()` tries to
#'   extract time from the object (e.g., `t=<time>` names, `attr(x, "time")`,
#'   `attr(x, "times")`, `x$time`, `x$times`, numeric dimnames).
#' @param orient One of "auto", "time_rows", "time_cols".
#' @param strict Logical. If TRUE, enforce strict checks on time grids.
#' @return A "survmat" object.
#' @export
as_survmat <- function(x, t = NULL,
                       orient = c("auto", "time_rows", "time_cols"),
                       strict = TRUE) {
  orient <- match.arg(orient)
  t_explicit <- t

  if (inherits(x, "survmat")) {
    return(x)
  }

  pick_named <- function(obj, keys) {
    for (k in keys) {
      if (!is.null(obj[[k]])) return(obj[[k]])
    }
    NULL
  }

  parse_numeric_names <- function(nm) {
    if (is.null(nm)) return(NULL)
    v <- suppressWarnings(as.numeric(nm))
    if (anyNA(v)) return(NULL)
    v
  }

  parse_t_style_names <- function(nm) {
    if (is.null(nm)) return(NULL)
    hit <- grepl("^t=", nm)
    if (!any(hit)) return(NULL)
    if (!all(hit)) {
      .stop("If time is encoded in names, all names must use the 't=<time>' format.")
    }
    tt <- suppressWarnings(as.numeric(sub("^t=", "", nm)))
    if (anyNA(tt) || any(!is.finite(tt)) || length(tt) < 2) {
      .stop("Could not parse a valid time grid from names formatted as 't=<time>'.")
    }
    tt
  }

  get_time_attr <- function(obj) {
    tt <- attr(obj, "time", exact = TRUE)
    if (is.null(tt)) tt <- attr(obj, "times", exact = TRUE)
    tt
  }

  ensure_numeric_time <- function(tt, label) {
    if (is.null(tt)) return(NULL)
    if (!is.numeric(tt)) .stop(sprintf("`%s` must be numeric.", label))
    if (any(!is.finite(tt))) .stop(sprintf("`%s` contains non-finite values.", label))
    as.numeric(tt)
  }

  maybe_check_same_time <- function(t_ref, t_other, msg) {
    if (!strict || is.null(t_ref) || is.null(t_other)) return(invisible(NULL))
    tol <- 1e-12
    if (length(t_ref) != length(t_other) || max(abs(t_ref - t_other)) > tol) {
      .stop(msg)
    }
    invisible(NULL)
  }

  as_survmat_matrix <- function(S_raw, t_from_source = NULL, source_label = "x") {
    if (!is.matrix(S_raw)) .stop(sprintf("`%s` must be a matrix.", source_label))
    if (!is.numeric(S_raw)) .stop(sprintf("`%s` must be numeric.", source_label))

    t_use <- ensure_numeric_time(t_explicit, "t")
    t_src <- ensure_numeric_time(t_from_source, "time grid")

    if (is.null(t_use)) {
      t_use <- t_src
    } else {
      maybe_check_same_time(
        t_use, t_src,
        "Provided `t` does not match time grid encoded in the input object."
      )
    }

    if (orient == "time_rows") {
      if (is.null(t_use) || length(t_use) != nrow(S_raw)) {
        .stop("`t` must match nrow(x) when orient = \"time_rows\".")
      }
      return(survmat(t_use, S_raw))
    }

    if (orient == "time_cols") {
      if (is.null(t_use) || length(t_use) != ncol(S_raw)) {
        .stop("`t` must match ncol(x) when orient = \"time_cols\".")
      }
      return(survmat(t_use, t(S_raw)))
    }

    # orient == "auto"
    if (!is.null(t_use)) {
      match_rows <- length(t_use) == nrow(S_raw)
      match_cols <- length(t_use) == ncol(S_raw)
      if (match_rows && !match_cols) return(survmat(t_use, S_raw))
      if (match_cols && !match_rows) return(survmat(t_use, t(S_raw)))
      if (match_rows && match_cols) {
        # deterministic default for matrix/list objects: time in rows
        return(survmat(t_use, S_raw))
      }
      .stop("`t` length must match nrow(x) or ncol(x) when orient = \"auto\".")
    }

    .stop("Could not determine a valid time grid. Provide `t`, `attr(x, \"time\"|\"times\")`, or numeric dimnames.")
  }

  if (is.matrix(x)) {
    t_attr <- get_time_attr(x)
    if (is.null(t_attr)) {
      if (orient == "time_rows") t_attr <- parse_numeric_names(rownames(x))
      if (orient == "time_cols") t_attr <- parse_numeric_names(colnames(x))
      if (orient == "auto") {
        t_attr <- parse_numeric_names(rownames(x))
        if (is.null(t_attr)) t_attr <- parse_numeric_names(colnames(x))
      }
    }
    return(as_survmat_matrix(x, t_from_source = t_attr, source_label = "x"))
  }

  if (is.data.frame(x)) {
    S_raw <- as.matrix(x)
    if (!is.numeric(S_raw)) .stop("`x` must contain numeric values.")

    t_cols <- parse_t_style_names(colnames(x))
    t_attr <- get_time_attr(x)

    if (!is.null(t_cols)) {
      maybe_check_same_time(
        ensure_numeric_time(t_explicit, "t"),
        t_cols,
        "Provided `t` does not match times encoded in data.frame column names."
      )
      return(survmat(t_cols, t(S_raw)))
    }

    t_pref <- ensure_numeric_time(t_explicit, "t")
    if (!is.null(t_pref)) {
      maybe_check_same_time(
        t_pref, ensure_numeric_time(t_attr, "time grid"),
        "Provided `t` does not match time grid encoded in the input object."
      )
    }
    if (is.null(t_pref)) t_pref <- ensure_numeric_time(t_attr, "time grid")

    if (is.null(t_pref)) {
      .stop("For data.frame inputs, provide `t`, 't=<time>' column names, or `attr(x, \"time\"|\"times\")`.")
    }

    if (orient == "time_rows") {
      if (length(t_pref) != nrow(S_raw)) .stop("`t` must match nrow(x) when orient = \"time_rows\".")
      return(survmat(t_pref, S_raw))
    }

    if (orient == "time_cols") {
      if (length(t_pref) != ncol(S_raw)) .stop("`t` must match ncol(x) when orient = \"time_cols\".")
      return(survmat(t_pref, t(S_raw)))
    }

    # orient == "auto"
    match_rows <- length(t_pref) == nrow(S_raw)
    match_cols <- length(t_pref) == ncol(S_raw)
    if (match_cols && !match_rows) return(survmat(t_pref, t(S_raw)))
    if (match_rows && !match_cols) return(survmat(t_pref, S_raw))
    if (match_rows && match_cols) {
      # deterministic default for data.frame objects: rows are units, cols are time
      return(survmat(t_pref, t(S_raw)))
    }
    .stop("`t` length must match nrow(x) or ncol(x) when orient = \"auto\".")
  }

  if (is.list(x)) {
    S_list <- pick_named(x, c("surv", "survival", "S"))
    t_list <- pick_named(x, c("time", "times"))
    if (is.null(S_list) || is.null(t_list)) {
      .stop("If `x` is a list, it must contain one of {surv, survival, S} and one of {time, times}.")
    }

    S_mat <- as.matrix(S_list)
    if (!is.numeric(S_mat)) .stop("List survival component must be numeric.")
    return(as_survmat_matrix(S_mat, t_from_source = t_list, source_label = "x[[surv|survival|S]]"))
  }

  .stop("`x` must be a survmat, matrix, data.frame, or list-like survival prediction object.")
}

align_survmat_pair <- function(sm1, sm0, tol = 1e-12) {
  if (length(sm1$t) != length(sm0$t) || max(abs(sm1$t - sm0$t)) > tol) {
    .stop(paste0(
      "Time grids differ. Evaluate both survival objects on the same grid ",
      "(e.g., pass the same `times=`) before comparing."
    ))
  }
  list(t = sm1$t, S1 = sm1$S, S0 = sm0$S)
}
