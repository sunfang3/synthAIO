#' Pre-treatment MSPE of W on the MSPE window
#'
#' @param Y1 Treated outcome path.
#' @param Y0 Donor outcome matrix (times × donors).
#' @param W Donor weights.
#' @param mspe_index Row index of `Y1` / `Y0` in the MSPE window.
#' @noRd
sc_pre_mspe <- function(Y1, Y0, W, mspe_index) {
  y1 <- as.numeric(Y1)[mspe_index]
  y0 <- as.matrix(Y0)[mspe_index, , drop = FALSE]
  mean((y1 - as.numeric(y0 %*% as.numeric(W)))^2)
}

#' Softmax map from unconstrained theta to a simplex V
#' @noRd
sc_softmax <- function(theta) {
  theta <- as.numeric(theta)
  if (!length(theta) || !all(is.finite(theta))) {
    return(rep(NA_real_, length(theta)))
  }
  shifted <- theta - max(theta)
  e <- exp(shifted)
  e / sum(e)
}

#' Synth / synth2 map: V_k = |theta_k| / sum |theta|
#' @noRd
sc_abs_norm <- function(theta) {
  theta <- as.numeric(theta)
  if (!length(theta) || !all(is.finite(theta))) {
    return(rep(NA_real_, length(theta)))
  }
  a <- abs(theta)
  s <- sum(a)
  if (s == 0) {
    return(rep(1 / length(theta), length(theta)))
  }
  a / s
}

#' Unconstrained start for abs-norm (identity on the simplex)
#' @noRd
sc_v_to_theta <- function(V) {
  as.numeric(V)
}

#' Regression starting V from Synth::synth (Abadie-Diamond-Hainmueller)
#'
#' Port of Synth::synth regression start (Abadie-Diamond-Hainmueller).
#' See Synth synth.R (`Xall` / `Zall` / `Beta` / `diag(V)`). Do not OLS
#' X1 on X0, do not average Z first, and do not floor tiny V at epsilon.
#' @noRd
sc_regression_v <- function(X1_scaled, X0_scaled, Z1, Z0) {
  X0_scaled <- as.matrix(X0_scaled)
  X1_scaled <- as.numeric(X1_scaled)
  K <- nrow(X0_scaled)
  J <- ncol(X0_scaled)
  pred_names <- rownames(X0_scaled)
  equal <- rep(1 / K, K)
  names(equal) <- pred_names

  Xall <- cbind(rep(1, 1 + J), t(cbind(X1_scaled, X0_scaled)))
  Zall <- cbind(Z1, Z0)
  B <- tryCatch(
    solve(crossprod(Xall), crossprod(Xall, t(Zall))),
    error = function(e) NULL
  )
  if (is.null(B)) {
    return(equal)
  }
  B <- B[-1, , drop = FALSE]
  v0 <- diag(tcrossprod(B))
  s <- sum(v0)
  if (!is.finite(s) || s == 0) {
    return(equal)
  }
  V <- v0 / s
  names(V) <- pred_names
  V
}

#' Normalize a user V onto the simplex
#' @noRd
sc_normalize_v <- function(custom_v, K, pred_names) {
  if (is.null(custom_v)) {
    synthaio_stop(
      "synthaio_bad_v",
      "`custom_v` is required when method = \"custom\""
    )
  }
  if (!is.numeric(custom_v) || length(custom_v) != K) {
    synthaio_stop(
      "synthaio_bad_v",
      "`custom_v` must be a numeric vector of length ", K
    )
  }
  if (any(!is.finite(custom_v))) {
    synthaio_stop("synthaio_bad_v", "`custom_v` must be finite")
  }
  if (any(custom_v < 0)) {
    synthaio_stop("synthaio_bad_v", "`custom_v` cannot contain negatives")
  }
  s <- sum(custom_v)
  if (s <= 0) {
    synthaio_stop(
      "synthaio_bad_v",
      "`custom_v` must sum to a positive value"
    )
  }
  V <- as.numeric(custom_v) / s
  names(V) <- pred_names
  V
}

#' Nested MSPE at one unconstrained theta
#' @noRd
sc_theta_mspe <- function(theta, spec, mspe_index, ...) {
  V <- sc_abs_norm(theta)
  if (!all(is.finite(V))) {
    return(.Machine$double.xmax)
  }
  W <- tryCatch(
    sc_solve_w(spec$X0_scaled, spec$X1_scaled, V, ...),
    error = function(e) NULL
  )
  if (is.null(W) || !all(is.finite(W))) {
    return(.Machine$double.xmax)
  }
  mspe <- sc_pre_mspe(spec$Y1, spec$Y0, W, mspe_index)
  if (!is.finite(mspe)) .Machine$double.xmax else mspe
}

#' Pack V, W, mspe for a feasible theta
#' @noRd
sc_theta_fit <- function(theta, spec, mspe_index, pred_names, ...) {
  V <- sc_abs_norm(theta)
  names(V) <- pred_names
  W <- sc_solve_w(spec$X0_scaled, spec$X1_scaled, V, ...)
  mspe <- sc_pre_mspe(spec$Y1, spec$Y0, W, mspe_index)
  list(V = V, W = W, mspe = mspe, theta = as.numeric(theta))
}

#' One nested run: Nelder-Mead, then BFGS polish; keep the best MSPE
#' @noRd
sc_nested_from_theta <- function(theta0,
                                 spec,
                                 mspe_index,
                                 maxit,
                                 fnscale = 1,
                                 ...) {
  pred_names <- rownames(spec$X0_scaled)
  obj <- function(theta) {
    sc_theta_mspe(theta, spec, mspe_index, ...)
  }
  best <- sc_theta_fit(theta0, spec, mspe_index, pred_names, ...)
  best$convergence <- NA_integer_

  consider <- function(opt) {
    if (is.null(opt) || !is.finite(opt$value)) {
      return(invisible(NULL))
    }
    fit <- sc_theta_fit(opt$par, spec, mspe_index, pred_names, ...)
    fit$convergence <- opt$convergence
    if (is.finite(fit$mspe) && fit$mspe <= best$mspe) {
      best <<- fit
    }
    invisible(NULL)
  }

  nm_maxit <- min(as.integer(maxit), 800L)
  bfgs_maxit <- min(as.integer(maxit), 200L)
  consider(tryCatch(
    stats::optim(
      par = as.numeric(theta0),
      fn = obj,
      method = "Nelder-Mead",
      control = list(
        maxit = nm_maxit,
        reltol = 1e-10,
        fnscale = fnscale
      )
    ),
    error = function(e) NULL
  ))
  consider(tryCatch(
    stats::optim(
      par = as.numeric(best$theta),
      fn = obj,
      method = "BFGS",
      control = list(
        maxit = bfgs_maxit,
        fnscale = fnscale
      )
    ),
    error = function(e) NULL
  ))
  best
}

#' Jitter theta with a documented seed and restore the RNG
#' @noRd
sc_jitter_theta <- function(theta, start_seed) {
  if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  } else {
    on.exit(
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      },
      add = TRUE
    )
  }
  set.seed(as.integer(start_seed))
  as.numeric(theta) + stats::rnorm(length(theta))
}

#' Estimate predictor weights V and the matching donor weights W
#'
#' @param spec An `sc_spec`.
#' @param method One of `"regression"` (default), `"nested"`, `"allopt"`,
#'   `"custom"`.
#' @param custom_v Length-`K` non-negative weights. Required for
#'   `method = "custom"`.
#' @param start_seed RNG seed for the jittered allopt start.
#' @param ... Passed to [sc_solve_w()] (`margin`, `maxiter`, `sigf`, `bound`).
#'
#' @return A list with `V`, `W`, `mspe`, and `starts`.
#' @noRd
sc_solve_v <- function(spec,
                       method = c("regression", "nested", "allopt", "custom"),
                       custom_v = NULL,
                       start_seed = 1,
                       ...) {
  if (!inherits(spec, "sc_spec")) {
    stop("`spec` must be an sc_spec object.", call. = FALSE)
  }
  method <- match.arg(method)

  pred_names <- rownames(spec$X0_scaled)
  K <- length(spec$X1_scaled)
  mspe_index <- match(spec$mspe_times, spec$times)
  dots <- list(...)
  # Guard the outer BFGS; `maxiter` in `...` is also the inner QP cap.
  outer_maxit <- dots$maxiter
  if (is.null(outer_maxit)) {
    outer_maxit <- 1000L
  }
  outer_maxit <- max(1L, as.integer(outer_maxit))

  Z1 <- spec$Y1[mspe_index]
  Z0 <- spec$Y0[mspe_index, , drop = FALSE]
  V_reg <- sc_regression_v(spec$X1_scaled, spec$X0_scaled, Z1, Z0)
  V_eq <- rep(1 / K, K)
  names(V_eq) <- pred_names

  pack <- function(fit, starts) {
    V <- as.numeric(fit$V)
    names(V) <- pred_names
    list(
      V = V,
      W = fit$W,
      mspe = as.numeric(fit$mspe),
      starts = starts
    )
  }

  finish <- function(V, label) {
    W <- sc_solve_w(spec$X0_scaled, spec$X1_scaled, V, ...)
    mspe <- sc_pre_mspe(spec$Y1, spec$Y0, W, mspe_index)
    list(V = V, W = W, mspe = mspe, label = label)
  }

  if (method == "custom") {
    V <- sc_normalize_v(custom_v, K, pred_names)
    fit <- finish(V, "custom")
    return(pack(fit, starts = NULL))
  }

  if (method == "regression") {
    fit <- finish(V_reg, "regression")
    return(pack(fit, starts = list(list(
      label = "regression",
      V = fit$V,
      W = fit$W,
      mspe = fit$mspe
    ))))
  }

  run_nested <- function(V0, label, fnscale = 1, theta = NULL) {
    if (is.null(theta)) {
      theta <- sc_v_to_theta(V0)
    }
    fit <- sc_nested_from_theta(
      theta,
      spec,
      mspe_index,
      maxit = outer_maxit,
      fnscale = fnscale,
      ...
    )
    fit$label <- label
    fit
  }

  if (method == "nested") {
    fit <- run_nested(V_reg, "regression")
    return(pack(fit, starts = list(list(
      label = fit$label,
      V = fit$V,
      W = fit$W,
      mspe = fit$mspe,
      convergence = fit$convergence
    ))))
  }

  # allopt: regression / equal / jitter plus singleton and
  # (k, last-predictor) pair starts. Keep the lowest MSPE.
  jitter_theta <- sc_jitter_theta(sc_v_to_theta(V_reg), start_seed)
  V_jitter <- sc_abs_norm(jitter_theta)
  names(V_jitter) <- pred_names
  mspe_jitter0 <- sc_pre_mspe(
    spec$Y1,
    spec$Y0,
    sc_solve_w(spec$X0_scaled, spec$X1_scaled, V_jitter, ...),
    mspe_index
  )
  fnscale_jitter <- mspe_jitter0
  if (!is.finite(fnscale_jitter) || fnscale_jitter <= 0) {
    fnscale_jitter <- 1
  }

  extra <- list()
  # Last-predictor singleton and (2, K) pair: extra basins on
  # smoking-like specs without baking in gold V.
  if (K >= 1L) {
    vk <- rep(0, K)
    vk[[K]] <- 1
    names(vk) <- pred_names
    extra[[length(extra) + 1L]] <- run_nested(vk, paste0("unit_", K))
  }
  if (K >= 2L) {
    vk <- rep(0, K)
    vk[[2L]] <- 0.5
    vk[[K]] <- 0.5
    names(vk) <- pred_names
    extra[[length(extra) + 1L]] <- run_nested(vk, paste0("pair_2_", K))
  }

  fits <- c(
    list(
      run_nested(V_reg, "regression"),
      run_nested(V_eq, "equal"),
      run_nested(V_jitter, "jitter", fnscale = fnscale_jitter, theta = jitter_theta)
    ),
    extra
  )
  starts <- lapply(fits, function(fit) {
    list(
      label = fit$label,
      V = fit$V,
      W = fit$W,
      mspe = fit$mspe,
      convergence = fit$convergence
    )
  })
  best <- which.min(vapply(fits, `[[`, numeric(1L), "mspe"))
  pack(fits[[best]], starts = starts)
}
