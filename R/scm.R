#' Map scm flags onto an [sc_fit()] method
#'
#' `allopt = TRUE` wins. A contradictory `nested = FALSE` is ignored
#' with a warning. `custom_v` is used only when neither nested flag
#' selected a method.
#' @noRd
scm_method <- function(nested, allopt, custom_v, nested_explicit = FALSE) {
  if (isTRUE(allopt)) {
    if (nested_explicit && identical(nested, FALSE)) {
      warning(
        "allopt = TRUE implies nested; ignoring nested = FALSE",
        call. = FALSE
      )
    }
    return("allopt")
  }
  if (isTRUE(nested)) {
    return("nested")
  }
  if (!is.null(custom_v)) {
    return("custom")
  }
  "regression"
}

#' Normalize the `placebo` argument into dispatcher inputs
#'
#' `NULL` and `list()` mean no placebo. Missing `cut` defaults to `Inf`.
#' Space runs when `unit` is `TRUE` or a vector of ids; time runs when
#' `period` is set.
#' @noRd
scm_placebo_args <- function(placebo) {
  if (is.null(placebo)) {
    return(NULL)
  }
  if (!is.list(placebo)) {
    stop("`placebo` must be NULL or a list.", call. = FALSE)
  }
  unit <- placebo$unit
  period <- placebo$period
  cut <- if (is.null(placebo$cut)) Inf else placebo$cut
  run_space <- !is.null(unit) && !identical(unit, FALSE)
  run_time <- !is.null(period)
  if (!run_space && !run_time) {
    return(NULL)
  }
  list(
    unit = if (run_space) unit else FALSE,
    period = period,
    cut = cut
  )
}

#' Estimate a synthetic control and optional diagnostics
#'
#' Builds an [sc_spec()], fits it with [sc_fit()], then optionally
#' attaches in-time / in-space placebos and leave-one-out robustness
#' to the same `scm_fit` object. Contains no estimation math of its
#' own.
#'
#' Method flags: `allopt = TRUE` selects `"allopt"` and implies nested.
#' Otherwise `nested = TRUE` selects `"nested"`. Otherwise a non-`NULL`
#' `custom_v` selects `"custom"`. Otherwise `"regression"`.
#'
#' `placebo` is `NULL` or a list with optional `unit` (`TRUE` or donor
#' ids), `period` (fake treatment time), and `cut` (pre-MSPE filter).
#' `placebo = list()` runs no placebo. If both `unit` and `period` are
#' set, in-time is run first, then in-space on the shifted spec.
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
#' @param nested Use nested MSPE minimization for \(V\).
#' @param allopt Multi-start nested optimization. Implies `nested`.
#' @param custom_v Length-\(K\) non-negative predictor weights.
#' @param placebo `NULL` or a list with `unit`, `period`, and/or `cut`.
#' @param loo If `TRUE`, attach leave-one-out donor robustness.
#' @param margin Extra primal-feasibility slack for the inner QP.
#' @param maxiter Inner QP iteration cap (also outer BFGS `maxit`).
#' @param sigf Solver significance digits; also the LOO weight threshold.
#' @param bound Upper bound on each donor weight.
#' @param start_seed RNG seed for the jittered allopt start.
#'
#' @return An object of class `scm_fit`. Optional slots `$placebo_time`,
#'   `$placebo_space`, and `$loo` are attached when requested.
#'
#' @examples
#' fit <- scm(
#'   cigsale ~ lnincome + cigsale(1988),
#'   data = smoking, unit = "state", time = "year",
#'   treat = 3, trperiod = 1989, xperiod = 1980:1988
#' )
#' fit$att
#'
#' @export
scm <- function(formula,
                data,
                unit,
                time,
                treat,
                trperiod,
                counit = NULL,
                xperiod = NULL,
                mspeperiod = NULL,
                preperiod = NULL,
                postperiod = NULL,
                nested = FALSE,
                allopt = FALSE,
                custom_v = NULL,
                placebo = NULL,
                loo = FALSE,
                margin = 0.05,
                maxiter = 1000,
                sigf = 7,
                bound = 10,
                start_seed = 1) {
  unit <- sc_spec_colname(substitute(unit))
  time <- sc_spec_colname(substitute(time))

  spec <- do.call(
    sc_spec,
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
      postperiod = postperiod
    )
  )

  method <- scm_method(
    nested = nested,
    allopt = allopt,
    custom_v = custom_v,
    nested_explicit = !missing(nested)
  )

  solver <- list(
    margin = margin,
    maxiter = maxiter,
    sigf = sigf,
    bound = bound,
    start_seed = start_seed
  )
  fit <- do.call(
    sc_fit,
    c(
      list(spec = spec, method = method, custom_v = custom_v),
      solver
    )
  )

  plc <- scm_placebo_args(placebo)
  if (!is.null(plc)) {
    plc_out <- do.call(
      sc_placebo,
      c(
        list(
          spec = spec,
          fit = fit,
          unit = plc$unit,
          period = plc$period,
          cut = plc$cut
        ),
        solver
      )
    )
    fit$placebo_time <- plc_out$placebo_time
    fit$placebo_space <- plc_out$placebo_space
  }

  if (isTRUE(loo)) {
    fit$loo <- do.call(sc_loo, c(list(spec = spec, fit = fit), solver))
  }

  fit
}
