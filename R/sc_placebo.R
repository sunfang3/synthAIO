#' Rebuild an `sc_spec` by overriding constructor inputs
#'
#' `sc_spec()` reads `unit` / `time` via `substitute()`, so callers must
#' rebuild with `do.call` and character column names — never
#' `sc_spec(..., unit = spec$unit, time = spec$time)`.
#' @noRd
sc_rebuild_spec <- function(spec, ...) {
  args <- list(
    formula = spec$formula,
    data = spec$data,
    unit = spec$unit,
    time = spec$time,
    treat = spec$treat,
    trperiod = spec$trperiod,
    counit = spec$counit,
    xperiod = spec$xperiod,
    mspeperiod = spec$mspeperiod,
    preperiod = spec$preperiod,
    postperiod = spec$postperiod
  )
  override <- list(...)
  args[names(override)] <- override
  do.call(sc_spec, args)
}

#' Reconstruct a formula term for `var` at `times`
#' @noRd
sc_placebo_term_expr <- function(var, times) {
  if (length(times) == 1L) {
    return(call(var, times))
  }
  if (length(times) > 1L && all(diff(as.numeric(times)) == 1)) {
    return(call(var, call(":", times[[1L]], times[[length(times)]])))
  }
  call(var, as.call(c(quote(c), as.list(times))))
}

#' Drop / truncate predictors that use times at or after fake `t0`
#' @noRd
sc_placebo_shift_formula <- function(formula, t0, drop_bare = FALSE) {
  terms <- sc_spec_split_terms(formula[[3L]])
  kept <- list()
  for (term in terms) {
    parsed <- sc_spec_parse_term(term)
    if (is.null(parsed$times)) {
      if (!drop_bare) {
        kept[[length(kept) + 1L]] <- term
      }
      next
    }
    keep_times <- parsed$times[parsed$times < t0]
    if (!length(keep_times)) {
      next
    }
    if (length(keep_times) == length(parsed$times)) {
      kept[[length(kept) + 1L]] <- term
    } else {
      kept[[length(kept) + 1L]] <- sc_placebo_term_expr(parsed$var, keep_times)
    }
  }
  if (!length(kept)) {
    return(NULL)
  }
  rhs <- kept[[1L]]
  if (length(kept) > 1L) {
    for (i in seq.int(2L, length(kept))) {
      rhs <- call("+", rhs, kept[[i]])
    }
  }
  out <- formula
  out[[3L]] <- rhs
  out
}

#' Fit one placebo spec with the same method as `fit`
#' @noRd
sc_placebo_fit <- function(spec, fit, ...) {
  dots <- list(...)
  args <- c(list(spec = spec, method = fit$method), dots)
  if (identical(fit$method, "custom") && !("custom_v" %in% names(dots))) {
    args$custom_v <- fit$V
  }
  do.call(sc_fit, args)
}

#' Year-wise left and two-sided placebo p-values
#'
#' Includes the treated unit. Empty kept-placebo set returns `NA` paths.
#' @noRd
sc_placebo_year_p <- function(effect, treat, kept, times) {
  time_names <- as.character(times)
  na_p <- stats::setNames(rep(NA_real_, length(times)), time_names)
  units <- colnames(effect)
  treat_name <- as.character(treat)
  keep_names <- units[kept]
  if (!treat_name %in% keep_names || length(keep_names) < 2L) {
    return(list(p_left = na_p, p_two = na_p))
  }
  E <- effect[, keep_names, drop = FALSE]
  e_tr <- E[, treat_name]
  list(
    p_left = stats::setNames(rowMeans(E <= e_tr), time_names),
    p_two = stats::setNames(rowMeans(abs(E) >= abs(e_tr)), time_names)
  )
}

#' In-space placebo: each donor as treated, remaining units as donors
#'
#' Locked Fisher p is `#{ratio >= ratio_tr} / N_units` including the
#' treated unit (no mid-p, no `+1` in the denominator). `cut` drops
#' placebos with pre MSPE `>` `cut` times the treated pre MSPE before
#' filtered ranking and year-wise p.
#'
#' Reserved mixed args: if `period` is set, run [sc_placebo_time()] first
#' and then this in-space loop on the shifted spec.
#'
#' @param spec An `sc_spec`.
#' @param fit The original `scm_fit` (same method is reused).
#' @param cut Pre-MSPE filter relative to the treated unit. Default `Inf`.
#' @param unit Reserved mixed flag; `TRUE` (default) runs the donor loop.
#' @param period Reserved mixed fake treatment time, or `NULL`.
#' @param ... Passed to [sc_fit()].
#' @return List with tidy `table` (`unit`, `pre_mspe`, `post_mspe`,
#'   `ratio`, `kept`) plus `p_unfiltered`, `p_filtered`, `p_left`, `p_two`.
#' @noRd
sc_placebo_space <- function(spec,
                             fit,
                             cut = Inf,
                             unit = TRUE,
                             period = NULL,
                             ...) {
  if (!inherits(spec, "sc_spec")) {
    stop("`spec` must be an sc_spec object.", call. = FALSE)
  }
  if (!inherits(fit, "scm_fit")) {
    stop("`fit` must be an scm_fit object.", call. = FALSE)
  }
  if (!is.numeric(cut) || length(cut) != 1L || is.na(cut) || cut < 0) {
    stop("`cut` must be a single non-negative number.", call. = FALSE)
  }

  if (!is.null(period)) {
    time_res <- sc_placebo_time(spec, period, method = fit$method, ...)
    spec <- time_res$spec
    fit <- time_res$fit
  }
  if (identical(unit, FALSE)) {
    stop("`unit` is FALSE; nothing to do for in-space placebos.", call. = FALSE)
  }

  units <- c(spec$treat, spec$donors)
  n_units <- length(units)
  pre_mspe <- numeric(n_units)
  post_mspe <- numeric(n_units)
  ratio <- numeric(n_units)
  effect <- matrix(
    NA_real_,
    nrow = length(spec$times),
    ncol = n_units,
    dimnames = list(as.character(spec$times), as.character(units))
  )

  pre_mspe[[1L]] <- fit$pre_mspe
  post_mspe[[1L]] <- fit$post_mspe
  ratio[[1L]] <- fit$mspe_ratio
  effect[, 1L] <- fit$effect[as.character(spec$times)]

  pool <- units
  for (i in seq.int(2L, n_units)) {
    j <- units[[i]]
    spec_j <- sc_rebuild_spec(
      spec,
      treat = j,
      counit = pool[pool != j]
    )
    fit_j <- sc_placebo_fit(spec_j, fit, ...)
    pre_mspe[[i]] <- fit_j$pre_mspe
    post_mspe[[i]] <- fit_j$post_mspe
    ratio[[i]] <- fit_j$mspe_ratio
    effect[, i] <- fit_j$effect[as.character(spec$times)]
  }

  threshold <- cut * pre_mspe[[1L]]
  kept <- c(TRUE, pre_mspe[-1L] <= threshold)
  n_placebo_kept <- sum(kept[-1L])
  if (n_placebo_kept == 0L) {
    warning(
      warningCondition(
        "cut dropped every placebo; year-wise p-values are NA",
        class = "synthaio_cut_empty"
      )
    )
  }

  table <- data.frame(
    unit = units,
    pre_mspe = pre_mspe,
    post_mspe = post_mspe,
    ratio = ratio,
    kept = kept,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  p_unfiltered <- sum(ratio >= ratio[[1L]]) / n_units
  if (n_placebo_kept == 0L) {
    p_filtered <- NA_real_
  } else {
    p_filtered <- sum(ratio[kept] >= ratio[[1L]]) / sum(kept)
  }
  year_p <- sc_placebo_year_p(effect, spec$treat, kept, spec$times)

  list(
    table = table,
    p_unfiltered = p_unfiltered,
    p_filtered = p_filtered,
    p_left = year_p$p_left,
    p_two = year_p$p_two,
    effect = effect,
    spec = spec,
    fit = fit
  )
}

#' In-time placebo: fake treatment at `period` < `trperiod`
#'
#' Point predictors with time `>= period` are dropped. Bare covariates
#' are re-averaged on `xperiod` times strictly before `period`. If no
#' predictors remain, errors with `synthaio_no_predictors`.
#'
#' @param spec An `sc_spec`.
#' @param period Fake treatment time, strictly before `spec$trperiod`.
#' @param ... Passed to [sc_fit()] (`method`, QP args).
#' @return List with the shifted `spec` / `fit`, `effect_fake` on
#'   `[period, original trperiod)`, and `effect_post` on the original
#'   post window at the new weights.
#' @noRd
sc_placebo_time <- function(spec, period, ...) {
  if (!inherits(spec, "sc_spec")) {
    stop("`spec` must be an sc_spec object.", call. = FALSE)
  }
  if (!is.numeric(period) || length(period) != 1L || is.na(period)) {
    stop("`period` must be a single time point.", call. = FALSE)
  }
  if (period >= spec$trperiod) {
    stop("`period` must be strictly before `trperiod`.", call. = FALSE)
  }

  xperiod <- spec$xperiod
  drop_bare <- FALSE
  if (!is.null(xperiod)) {
    xperiod <- xperiod[xperiod < period]
    drop_bare <- !length(xperiod)
    if (drop_bare) {
      xperiod <- NULL
    }
  }

  formula <- sc_placebo_shift_formula(
    spec$formula, period, drop_bare = drop_bare
  )
  if (is.null(formula)) {
    synthaio_stop(
      "synthaio_no_predictors",
      "no predictors remain after shifting trperiod to ", period
    )
  }

  new_spec <- sc_rebuild_spec(
    spec,
    formula = formula,
    trperiod = period,
    xperiod = xperiod
  )
  new_fit <- sc_fit(new_spec, ...)

  times <- spec$times
  fake_times <- times[times >= period & times < spec$trperiod]
  effect <- new_fit$effect
  effect_fake <- effect[as.character(fake_times)]
  effect_post <- effect[as.character(spec$post_times)]

  list(
    period = period,
    spec = new_spec,
    fit = new_fit,
    effect = effect,
    effect_fake = effect_fake,
    effect_post = effect_post
  )
}

#' Placebo dispatcher: in-time, in-space, or mixed (time then space)
#'
#' If both `unit` and `period` are set, shift `trperiod` first and run
#' in-space placebos on that shifted spec.
#'
#' @param spec An `sc_spec`.
#' @param fit An `scm_fit`.
#' @param unit `TRUE` to run in-space placebos.
#' @param period Fake treatment time, or `NULL`.
#' @param cut Passed to [sc_placebo_space()].
#' @param ... Passed to the placebo workers.
#' @noRd
sc_placebo <- function(spec, fit, unit = FALSE, period = NULL, cut = Inf, ...) {
  time_res <- NULL
  space_spec <- spec
  space_fit <- fit
  if (!is.null(period)) {
    time_res <- sc_placebo_time(spec, period, method = fit$method, ...)
    space_spec <- time_res$spec
    space_fit <- time_res$fit
  }
  space_res <- NULL
  if (!identical(unit, FALSE) && !is.null(unit)) {
    space_res <- sc_placebo_space(space_spec, space_fit, cut = cut, ...)
  }
  list(placebo_time = time_res, placebo_space = space_res)
}
