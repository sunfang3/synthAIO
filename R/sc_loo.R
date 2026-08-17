#' Leave-one-out donor robustness
#'
#' Drops each donor with \(W_j > 10^{-sigf}\) in turn, rebuilds the spec
#' on the remaining donors, and refits with the same method. Returns the
#' original effect plus a period-wise min/max band over the LOO paths.
#'
#' `sigf` is taken from `...` when supplied, otherwise `7`.
#'
#' @param spec An `sc_spec`.
#' @param fit An `scm_fit`.
#' @param ... Passed to [sc_fit()] (`method` is taken from `fit`). `sigf`
#'   also sets the nonzero-weight threshold.
#' @return List with `band` (`time`, `effect`, `loo_min`, `loo_max`) and
#'   `fits`, a named list of per-dropped-unit `W`.
#' @noRd
sc_loo <- function(spec, fit, ...) {
  if (!inherits(spec, "sc_spec")) {
    stop("`spec` must be an sc_spec object.", call. = FALSE)
  }
  if (!inherits(fit, "scm_fit")) {
    stop("`fit` must be an scm_fit object.", call. = FALSE)
  }

  dots <- list(...)
  sigf <- if ("sigf" %in% names(dots)) dots$sigf else 7
  dropped <- names(fit$W)[as.numeric(fit$W) > 10^(-sigf)]
  if (!length(dropped)) {
    synthaio_stop(
      "synthaio_loo_empty",
      "no donor weight exceeds 10^(-sigf); nothing to drop"
    )
  }

  times <- spec$times
  time_names <- as.character(times)
  loo_effect <- matrix(
    NA_real_,
    nrow = length(times),
    ncol = length(dropped),
    dimnames = list(time_names, dropped)
  )
  fits <- vector("list", length(dropped))
  names(fits) <- dropped

  donors <- spec$donors
  for (j in dropped) {
    spec_j <- sc_rebuild_spec(
      spec,
      counit = donors[as.character(donors) != j]
    )
    fit_j <- sc_placebo_fit(spec_j, fit, ...)
    fits[[j]] <- fit_j$W
    loo_effect[, j] <- fit_j$effect[time_names]
  }

  list(
    band = data.frame(
      time = times,
      effect = as.numeric(fit$effect[time_names]),
      loo_min = apply(loo_effect, 1L, min),
      loo_max = apply(loo_effect, 1L, max),
      row.names = NULL,
      stringsAsFactors = FALSE
    ),
    fits = fits
  )
}
