#' Stop with a named synthaio error class
#'
#' @param class Error class name, e.g. `"synthaio_bad_time"`.
#' @param ... Message pieces pasted together.
#' @noRd
synthaio_stop <- function(class, ...) {
  stop(errorCondition(paste0(...), class = class))
}

#' Optional `{unit}_name` column next to a numeric unit id
#' @noRd
sc_spec_name_col <- function(unit) {
  paste0(unit, "_name")
}

#' Map a character treat name to the unit id via `{unit}_name`
#' @noRd
sc_spec_resolve_treat <- function(treat, data, unit) {
  present <- unique(data[[unit]])
  if (length(treat) == 1L && treat %in% present) {
    return(treat)
  }
  name_col <- sc_spec_name_col(unit)
  if (name_col %in% names(data)) {
    hit <- unique(data[[unit]][as.character(data[[name_col]]) == as.character(treat)])
    hit <- hit[!is.na(hit)]
    if (length(hit) == 1L) {
      return(hit[[1L]])
    }
  }
  synthaio_stop(
    "synthaio_missing_treat",
    "treated unit ", treat, " is not in `data`"
  )
}

#' Named id → label map from `{unit}_name`, or NULL
#' @noRd
sc_spec_unit_labels <- function(data, unit) {
  name_col <- sc_spec_name_col(unit)
  if (!name_col %in% names(data)) {
    return(NULL)
  }
  ids <- data[[unit]]
  labs <- as.character(data[[name_col]])
  keep <- !duplicated(ids) & !is.na(ids)
  stats::setNames(labs[keep], as.character(ids[keep]))
}

#' Format unit ids as "California (3)" when labels exist
#' @noRd
sc_format_unit <- function(spec, ids) {
  ids_chr <- as.character(ids)
  labels <- spec$unit_labels
  if (is.null(labels) || !length(labels)) {
    return(ids_chr)
  }
  lab <- unname(labels[ids_chr])
  ifelse(
    is.na(lab) | !nzchar(lab) | lab == ids_chr,
    ids_chr,
    paste0(lab, " (", ids_chr, ")")
  )
}

#' Coerce a `unit` / `time` argument to a column name
#' @noRd
sc_spec_colname <- function(expr) {
  if (is.character(expr) && length(expr) == 1L && !is.na(expr) && nzchar(expr)) {
    return(expr)
  }
  if (is.symbol(expr)) {
    return(as.character(expr))
  }
  stop("`unit` and `time` must be column names (strings or unquoted names).",
       call. = FALSE)
}

#' Split a formula RHS on `+`
#' @noRd
sc_spec_split_terms <- function(expr) {
  if (is.call(expr) && identical(expr[[1L]], quote(`(`))) {
    return(sc_spec_split_terms(expr[[2L]]))
  }
  if (is.call(expr) && identical(expr[[1L]], quote(`+`))) {
    return(c(sc_spec_split_terms(expr[[2L]]), sc_spec_split_terms(expr[[3L]])))
  }
  list(expr)
}

#' Parse one predictor term into variable, times, and matrix name
#' @noRd
sc_spec_parse_term <- function(expr) {
  if (is.symbol(expr)) {
    var <- as.character(expr)
    return(list(var = var, times = NULL, name = var))
  }
  if (!is.call(expr) || length(expr) < 2L) {
    stop("cannot parse predictor term: ", deparse1(expr), call. = FALSE)
  }
  fun <- expr[[1L]]
  if (!is.symbol(fun)) {
    stop("cannot parse predictor term: ", deparse1(expr), call. = FALSE)
  }
  var <- as.character(fun)
  arg <- expr[[2L]]
  if (is.numeric(arg) && length(arg) == 1L) {
    times <- arg
    name <- paste(var, times, sep = "_")
  } else if (is.call(arg) && identical(arg[[1L]], quote(`:`))) {
    times <- eval(arg, envir = baseenv())
    name <- paste(var, times[1L], times[length(times)], sep = "_")
  } else if (is.call(arg) && identical(arg[[1L]], quote(c))) {
    times <- eval(arg, envir = baseenv())
    name <- paste(c(var, times), collapse = "_")
  } else {
    times <- eval(arg, envir = baseenv())
    if (length(times) == 1L) {
      name <- paste(var, times, sep = "_")
    } else if (length(times) > 1L && all(diff(as.numeric(times)) == 1)) {
      name <- paste(var, times[1L], times[length(times)], sep = "_")
    } else {
      name <- paste(c(var, times), collapse = "_")
    }
  }
  if (!is.numeric(times) || !length(times) || anyNA(times)) {
    stop("predictor times must be numeric: ", deparse1(expr), call. = FALSE)
  }
  list(var = var, times = times, name = name)
}

#' Error if any requested times are absent from the panel
#' @noRd
sc_spec_check_times <- function(values, panel_times, what) {
  missing <- setdiff(values, panel_times)
  if (length(missing)) {
    synthaio_stop(
      "synthaio_bad_time",
      what, " includes times not in the panel: ",
      paste(missing, collapse = ", ")
    )
  }
}

#' Per-unit mean of `var` over `times`
#'
#' Averages use `na.rm = TRUE` so missing beer / income / age years
#' (Synth-style incomplete covariates) drop out of the mean.
#' @noRd
sc_spec_unit_means <- function(panel, unit, time, var, ids, times) {
  out <- numeric(length(ids))
  names(out) <- as.character(ids)
  u <- panel[[unit]]
  tms <- panel[[time]]
  v <- panel[[var]]
  for (i in seq_along(ids)) {
    out[[i]] <- mean(v[u == ids[[i]] & tms %in% times], na.rm = TRUE)
  }
  out
}

#' Outcome path for one unit, aligned to `times`
#' @noRd
sc_spec_unit_path <- function(panel, unit, time, var, id, times) {
  rows <- panel[[unit]] == id
  panel[[var]][rows][match(times, panel[[time]][rows])]
}

#' Parse a long panel plus a synth2-like formula into SCM matrices
#'
#' Turns a long unit-time data frame and a formula into the predictor and
#' outcome matrices an estimator can consume. Bare covariates are averaged
#' over `xperiod`. Calls such as `cigsale(1988)`, `cigsale(1980:1988)`, and
#' `cigsale(c(1982, 1986, 1988))` select those times and override `xperiod`
#' for that predictor only.
#'
#' Window arguments default independently: `preperiod`, `xperiod`, and
#' `mspeperiod` each default to all times `< trperiod`; `postperiod` defaults
#' to all times `>= trperiod`. Setting `xperiod` does **not** change
#' `mspeperiod`.
#'
#' @param formula A formula. The left-hand side is the outcome. Right-hand
#'   side terms are bare covariate names or `var(times)` selectors.
#' @param data A long panel data frame.
#' @param unit Column name of the unit id (string or unquoted name).
#' @param time Column name of the time index (string or unquoted name).
#' @param treat Treated unit id, matching `data[[unit]]`.
#' @param trperiod Treatment-period time point.
#' @param counit Optional vector of donor unit ids. The treated unit, if
#'   listed, is dropped with a warning.
#' @param xperiod Times used to average bare covariates. Default: all times
#'   `< trperiod`.
#' @param mspeperiod Times used for the MSPE / nested objective. Default:
#'   all times `< trperiod`, independent of `xperiod`.
#' @param preperiod Pre-treatment window (RMSE, \(R^2\)). Default: all times
#'   `< trperiod`.
#' @param postperiod Post-treatment window (ATT, post MSPE). Default: all
#'   times `>= trperiod`.
#'
#' @return An object of class `sc_spec`: a list with matrices `X0`, `X1`,
#'   `Y0`, `Y1`, `X0_scaled`, `X1_scaled`; window vectors `donors`, `times`,
#'   `pre_times`, `post_times`, `mspe_times`, `x_times`; slots `outcome`,
#'   `treat`, `trperiod`; and the constructor inputs so later steps can
#'   rebuild a spec.
#'
#' @examples
#' spec <- sc_spec(
#'   cigsale ~ lnincome + cigsale(1988),
#'   data = smoking, unit = "state", time = "year",
#'   treat = 3, trperiod = 1989, xperiod = 1980:1988
#' )
#' dim(spec$X0)
#'
#' @export
sc_spec <- function(formula,
                    data,
                    unit,
                    time,
                    treat,
                    trperiod,
                    counit = NULL,
                    xperiod = NULL,
                    mspeperiod = NULL,
                    preperiod = NULL,
                    postperiod = NULL) {
  unit <- sc_spec_colname(substitute(unit))
  time <- sc_spec_colname(substitute(time))

  if (!is.data.frame(data) || nrow(data) == 0L) {
    synthaio_stop("synthaio_not_panel", "`data` is empty or not a panel")
  }
  if (!all(c(unit, time) %in% names(data))) {
    stop("`unit` and `time` must name columns in `data`.", call. = FALSE)
  }
  if (anyDuplicated(data[c(unit, time)])) {
    synthaio_stop("synthaio_not_panel", "duplicate unit-time pairs")
  }
  if (!is.numeric(trperiod) || length(trperiod) != 1L || is.na(trperiod)) {
    stop("`trperiod` must be a single time point.", call. = FALSE)
  }
  if (!inherits(formula, "formula") || length(formula) < 3L) {
    stop("`formula` must be of the form outcome ~ predictors.", call. = FALSE)
  }
  lhs <- formula[[2L]]
  if (!is.symbol(lhs)) {
    stop("left-hand side of `formula` must be a single outcome name.",
         call. = FALSE)
  }
  outcome <- as.character(lhs)
  preds <- lapply(sc_spec_split_terms(formula[[3L]]), sc_spec_parse_term)
  needed <- unique(c(outcome, vapply(preds, `[[`, character(1L), "var")))
  missing_cols <- setdiff(needed, names(data))
  if (length(missing_cols)) {
    stop("missing columns: ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }

  treat <- sc_spec_resolve_treat(treat, data, unit)
  present <- unique(data[[unit]])
  unit_labels <- sc_spec_unit_labels(data, unit)
  if (is.null(counit)) {
    donors <- sort(present[present != treat])
  } else {
    if (any(counit == treat)) {
      warning(
        "treated unit listed in `counit` is dropped from the donor pool",
        call. = FALSE
      )
    }
    donors <- unique(counit)
    donors <- donors[donors != treat & donors %in% present]
  }
  if (length(donors) < 2L) {
    synthaio_stop(
      "synthaio_single_donor",
      "need at least 2 donor units, got ", length(donors)
    )
  }

  ids <- c(treat, donors)
  panel <- data[data[[unit]] %in% ids, , drop = FALSE]
  times <- sort(unique(panel[[time]]))
  if (nrow(panel) != length(ids) * length(times)) {
    synthaio_stop(
      "synthaio_unbalanced",
      "panel is unbalanced: a required unit-time is missing"
    )
  }

  panel_times <- unique(data[[time]])
  default_pre <- times[times < trperiod]
  default_post <- times[times >= trperiod]
  x_times <- if (is.null(xperiod)) default_pre else xperiod
  mspe_times <- if (is.null(mspeperiod)) default_pre else mspeperiod
  pre_times <- if (is.null(preperiod)) default_pre else preperiod
  post_times <- if (is.null(postperiod)) default_post else postperiod
  sc_spec_check_times(x_times, panel_times, "`xperiod`")
  sc_spec_check_times(mspe_times, panel_times, "`mspeperiod`")
  sc_spec_check_times(pre_times, panel_times, "`preperiod`")
  sc_spec_check_times(post_times, panel_times, "`postperiod`")
  for (pred in preds) {
    if (!is.null(pred$times)) {
      sc_spec_check_times(pred$times, panel_times, pred$name)
    }
  }

  pred_names <- vapply(preds, `[[`, character(1L), "name")
  X <- matrix(
    NA_real_,
    nrow = length(preds),
    ncol = length(ids),
    dimnames = list(pred_names, as.character(ids))
  )
  for (k in seq_along(preds)) {
    pred_times <- preds[[k]]$times
    if (is.null(pred_times)) {
      pred_times <- x_times
    }
    X[k, ] <- sc_spec_unit_means(
      panel, unit, time, preds[[k]]$var, ids, pred_times
    )
  }
  # A 1-row X[, 1] drops to an unnamed scalar; keep predictor names.
  X1 <- stats::setNames(as.numeric(X[, 1L]), rownames(X))
  X0 <- X[, -1L, drop = FALSE]

  Y <- matrix(
    NA_real_,
    nrow = length(times),
    ncol = length(ids),
    dimnames = list(as.character(times), as.character(ids))
  )
  for (j in seq_along(ids)) {
    Y[, j] <- sc_spec_unit_path(
      panel, unit, time, outcome, ids[[j]], times
    )
  }
  Y1 <- stats::setNames(as.numeric(Y[, 1L]), rownames(Y))
  Y0 <- Y[, -1L, drop = FALSE]

  X_all <- cbind(X0, X1)
  sds <- apply(X_all, 1L, stats::sd)
  scale <- sds
  scale[!is.finite(scale) | scale == 0] <- 1
  X0_scaled <- X0 / scale
  X1_scaled <- X1 / scale

  structure(
    list(
      formula = formula,
      data = data,
      unit = unit,
      time = time,
      treat = treat,
      trperiod = trperiod,
      counit = counit,
      xperiod = xperiod,
      mspeperiod = mspeperiod,
      preperiod = preperiod,
      postperiod = postperiod,
      donors = donors,
      times = times,
      pre_times = pre_times,
      post_times = post_times,
      mspe_times = mspe_times,
      x_times = x_times,
      outcome = outcome,
      unit_labels = unit_labels,
      X0 = X0,
      X1 = X1,
      Y0 = Y0,
      Y1 = Y1,
      X0_scaled = X0_scaled,
      X1_scaled = X1_scaled
    ),
    class = "sc_spec"
  )
}
