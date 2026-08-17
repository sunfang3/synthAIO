#' Stop unless `x` is an `scm_fit`
#' @noRd
sc_result_require <- function(x) {
  if (!inherits(x, "scm_fit")) {
    stop("`x` must be an scm_fit object.", call. = FALSE)
  }
  invisible(x)
}

#' Nonzero donor weights as a two-column table
#' @noRd
sc_result_nonzero_w <- function(W) {
  w <- as.numeric(W)
  names(w) <- names(W)
  nz <- w[w > 0]
  data.frame(
    unit = names(nz),
    weight = as.numeric(nz),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Print an `scm_fit`
#'
#' Reports the treated unit, treatment time, method, pre-treatment RMSE,
#' pre-treatment \(R^2\), ATT, and donor-pool size. Does not recompute
#' estimates.
#'
#' @param x An `scm_fit`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.scm_fit <- function(x, ...) {
  cat("Synthetic control\n")
  cat("  Treated:    ", x$spec$treat, "\n", sep = "")
  cat("  Time:       ", x$spec$trperiod, "\n", sep = "")
  cat("  Method:     ", x$method, "\n", sep = "")
  cat("  RMSE:       ", format(x$rmse, digits = 4), "\n", sep = "")
  cat("  R-squared:  ", format(x$r2, digits = 4), "\n", sep = "")
  cat("  ATT:        ", format(x$att, digits = 4), "\n", sep = "")
  cat("  Donors:     ", length(x$spec$donors), "\n", sep = "")
  invisible(x)
}

#' Summarize an `scm_fit`
#'
#' Prints the `print()` header plus the predictor balance table and
#' nonzero donor weights. If an in-space placebo is attached, also
#' prints the MSPE table and Fisher p-values. Without a placebo,
#' p-values are not mentioned. Does not recompute estimates.
#'
#' @param object An `scm_fit`.
#' @param ... Passed to `print()`.
#' @return `object`, invisibly.
#' @export
summary.scm_fit <- function(object, ...) {
  print(object, ...)

  cat("\nPredictor balance\n")
  print(object$balance, row.names = FALSE)

  cat("\nDonor weights (nonzero)\n")
  wtab <- sc_result_nonzero_w(object$W)
  if (!nrow(wtab)) {
    cat("(none)\n")
  } else {
    print(wtab, row.names = FALSE)
  }

  plc <- object$placebo_space
  if (!is.null(plc)) {
    if (!is.null(plc$table)) {
      cat("\nPlacebo MSPE\n")
      print(plc$table, row.names = FALSE)
    }
    if (!is.null(plc$p_unfiltered)) {
      cat(
        "\nFisher p-value (unfiltered): ",
        format(plc$p_unfiltered, digits = 4),
        "\n",
        sep = ""
      )
    }
    if (!is.null(plc$p_filtered)) {
      cat(
        "Fisher p-value (filtered):   ",
        format(plc$p_filtered, digits = 4),
        "\n",
        sep = ""
      )
    }
  }

  invisible(object)
}

#' Extract pieces of an `scm_fit`
#'
#' These helpers return stored fields. They do not refit or recompute
#' estimates. `glance()` is defined here (no **broom** dependency).
#'
#' @param x An `scm_fit`.
#' @param object An `scm_fit`.
#' @param ... Unused.
#'
#' @return
#' * `coef()`: donor weights \(W\) (names are donor ids).
#' * `v_weights()`: predictor weights \(V\).
#' * `balance()`: the predictor balance data frame.
#' * `effects()`: the treated-minus-synthetic gap series.
#' * `glance()`: one-row data frame of scalars (`treated`, `trperiod`,
#'   `method`, `n_donors`, `att`, `rmse`, `r2`, `pre_mspe`, `post_mspe`,
#'   `mspe_ratio`).
#'
#' @examples
#' spec <- sc_spec(
#'   cigsale ~ lnincome + cigsale(1988),
#'   data = smoking, unit = "state", time = "year",
#'   treat = 3, trperiod = 1989, xperiod = 1980:1988
#' )
#' fit <- sc_fit(spec, method = "regression")
#' coef(fit)
#' v_weights(fit)
#' glance(fit)
#'
#' @name scm_fit-extractors
NULL

#' @export
#' @rdname scm_fit-extractors
coef.scm_fit <- function(object, ...) {
  object$W
}

#' @export
#' @rdname scm_fit-extractors
v_weights <- function(x, ...) {
  sc_result_require(x)
  x$V
}

#' @export
#' @rdname scm_fit-extractors
balance <- function(x, ...) {
  sc_result_require(x)
  x$balance
}

#' @export
#' @rdname scm_fit-extractors
effects <- function(x, ...) {
  sc_result_require(x)
  x$effect
}

#' @export
#' @rdname scm_fit-extractors
glance <- function(x, ...) {
  sc_result_require(x)
  data.frame(
    treated = x$spec$treat,
    trperiod = x$spec$trperiod,
    method = x$method,
    n_donors = length(x$spec$donors),
    att = unname(x$att),
    rmse = unname(x$rmse),
    r2 = unname(x$r2),
    pre_mspe = unname(x$pre_mspe),
    post_mspe = unname(x$post_mspe),
    mspe_ratio = unname(x$mspe_ratio),
    stringsAsFactors = FALSE
  )
}
