#' Predictor percent bias; `NA` when the treated mean is zero
#' @noRd
sc_bias_pct <- function(value, treated) {
  out <- rep(NA_real_, length(treated))
  ok <- treated != 0
  out[ok] <- 100 * (value[ok] - treated[ok]) / treated[ok]
  out
}

#' Fit a synthetic control for one specification
#'
#' Solves predictor weights \eqn{V} and donor weights \eqn{W}, then
#' builds treated/synthetic paths, ATT, pre-treatment fit, MSPE
#' ratio, and a predictor balance table.
#'
#' @param spec An `sc_spec`.
#' @param method One of `"regression"` (default), `"nested"`, `"allopt"`,
#'   `"custom"`.
#' @param custom_v Length-`K` non-negative weights. Required for
#'   `method = "custom"`.
#' @param ... Passed to the \(V\) solver (`start_seed`, and the inner
#'   QP arguments `margin`, `maxiter`, `sigf`, `bound`).
#'
#' @return An object of class `scm_fit`: a list with `spec`, `method`,
#'   `V`, `W`, `y_treated`, `y_synth`, `effect`, `att`, `rmse`, `r2`,
#'   `pre_mspe`, `post_mspe`, `mspe_ratio`, and `balance`.
#'
#' @examples
#' spec <- sc_spec(
#'   cigsale ~ lnincome + cigsale(1988),
#'   data = smoking, unit = "state", time = "year",
#'   treat = 3, trperiod = 1989, xperiod = 1980:1988
#' )
#' fit <- sc_fit(spec, method = "regression")
#' fit$att
#'
#' @export
sc_fit <- function(spec, method = "regression", custom_v = NULL, ...) {
  if (!inherits(spec, "sc_spec")) {
    stop("`spec` must be an sc_spec object.", call. = FALSE)
  }
  method <- match.arg(method, c("regression", "nested", "allopt", "custom"))

  solved <- sc_solve_v(spec, method = method, custom_v = custom_v, ...)
  V <- solved$V
  W <- solved$W

  y_treated <- as.numeric(spec$Y1)
  names(y_treated) <- names(spec$Y1)
  y_synth <- as.numeric(as.matrix(spec$Y0) %*% as.numeric(W))
  names(y_synth) <- names(spec$Y1)
  effect <- y_treated - y_synth

  pre_idx <- match(spec$pre_times, spec$times)
  post_idx <- match(spec$post_times, spec$times)
  mspe_idx <- match(spec$mspe_times, spec$times)

  att <- mean(effect[post_idx])
  rmse <- sqrt(mean(effect[pre_idx]^2))
  y_pre <- y_treated[pre_idx]
  ss_res <- sum((y_pre - y_synth[pre_idx])^2)
  ss_tot <- sum((y_pre - mean(y_pre))^2)
  r2 <- if (!is.finite(ss_tot) || ss_tot == 0) {
    NA_real_
  } else {
    1 - ss_res / ss_tot
  }
  pre_mspe <- mean(effect[mspe_idx]^2)
  post_mspe <- mean(effect[post_idx]^2)
  mspe_ratio <- post_mspe / pre_mspe

  treated <- as.numeric(spec$X1)
  names(treated) <- names(spec$X1)
  synthetic <- as.numeric(as.matrix(spec$X0) %*% as.numeric(W))
  avg_control <- as.numeric(rowMeans(as.matrix(spec$X0)))
  balance <- data.frame(
    predictor = names(spec$X1),
    v_weight = as.numeric(V[names(spec$X1)]),
    treated = treated,
    synthetic = synthetic,
    bias_pct = sc_bias_pct(synthetic, treated),
    avg_control = avg_control,
    avg_bias_pct = sc_bias_pct(avg_control, treated),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  structure(
    list(
      spec = spec,
      method = method,
      V = V,
      W = W,
      y_treated = y_treated,
      y_synth = y_synth,
      effect = effect,
      att = att,
      rmse = rmse,
      r2 = r2,
      pre_mspe = pre_mspe,
      post_mspe = post_mspe,
      mspe_ratio = mspe_ratio,
      balance = balance
    ),
    class = "scm_fit"
  )
}
