#' Vanderbei LOQO interior-point QP (kernlab::ipop, no S4)
#'
#' min c'x + 1/2 x' H x  s.t.  b <= A x <= b+r,  l <= x <= u.
#' Used by Stata synth / synth2; osqp is the fallback if a
#' Newton system is singular.
#' @noRd
sc_ipop <- function(c,
                    H,
                    A,
                    b,
                    l,
                    u,
                    r = 0,
                    sigf = 7,
                    maxiter = 40,
                    margin = 0.05,
                    bound = 10) {
  if (is.vector(A)) {
    A <- matrix(A, 1L)
  }
  c <- matrix(as.numeric(c), ncol = 1L)
  l <- matrix(as.numeric(l), ncol = 1L)
  u <- matrix(as.numeric(u), ncol = 1L)
  b <- matrix(as.numeric(b), ncol = 1L)
  n <- ncol(A)
  m <- nrow(A)
  H.diag <- diag(H)
  H.x <- H
  b.plus.1 <- max(abs(b)) + 1
  c.plus.1 <- max(abs(c)) + 1
  H.y <- diag(1, m)
  diag(H.x) <- H.diag + 1
  AP <- matrix(0, m + n, m + n)
  xp <- seq_len(m + n) <= n
  AP[xp, xp] <- -H.x
  AP[!xp, xp] <- A
  AP[xp, !xp] <- t(A)
  AP[!xp, !xp] <- H.y
  s.tmp <- solve(AP, c(c, b), tol = .Machine$double.eps)
  x <- s.tmp[seq_len(n)]
  y <- s.tmp[-seq_len(n)]
  g <- pmax(abs(x - l), bound)
  z <- pmax(abs(x), bound)
  t <- pmax(abs(u - x), bound)
  s <- pmax(abs(x), bound)
  v <- pmax(abs(y), bound)
  w <- pmax(abs(y), bound)
  p <- pmax(abs(r - w), bound)
  q <- pmax(abs(y), bound)
  mu <- as.vector(
    crossprod(z, g) + crossprod(v, w) + crossprod(s, t) + crossprod(p, q)
  ) / (2 * (m + n))
  counter <- 0L
  while (counter < maxiter) {
    counter <- counter + 1L
    H.dot.x <- H %*% x
    rho <- b - A %*% x + w
    nu <- l - x + g
    tau <- u - x - t
    alpha <- r - w - p
    sigma <- c - crossprod(A, y) - z + s + H.dot.x
    beta <- y + q - v
    gamma.z <- -z
    gamma.w <- -w
    gamma.s <- -s
    gamma.q <- -q
    x.dot.H.dot.x <- crossprod(x, H.dot.x)
    primal.obj <- crossprod(c, x) + 0.5 * x.dot.H.dot.x
    dual.obj <- crossprod(b, y) - 0.5 * x.dot.H.dot.x +
      crossprod(l, z) - crossprod(u, s) - crossprod(r, q)
    sigfig <- max(
      -log10(abs(primal.obj - dual.obj) / (abs(primal.obj) + 1)),
      0
    )
    if (is.finite(sigfig) && sigfig >= sigf) {
      break
    }
    hat.beta <- beta - v * gamma.w / w
    hat.alpha <- alpha - p * gamma.q / q
    hat.nu <- nu + g * gamma.z / z
    hat.tau <- tau - t * gamma.s / s
    d <- as.numeric(z / g + s / t)
    e <- as.numeric(1 / (v / w + q / p))
    diag(H.x) <- H.diag + d
    if (length(e) == 1L) {
      H.y <- matrix(e, 1L, 1L)
    } else {
      diag(H.y) <- e
    }
    c.x <- sigma - z * hat.nu / g - s * hat.tau / t
    c.y <- rho - e * (hat.beta - q * hat.alpha / p)
    AP[xp, xp] <- -H.x
    AP[!xp, !xp] <- H.y
    s1.tmp <- solve(AP, c(c.x, c.y), tol = .Machine$double.eps)
    delta.x <- s1.tmp[seq_len(n)]
    delta.y <- s1.tmp[-seq_len(n)]
    delta.w <- -e * (hat.beta - q * hat.alpha / p + delta.y)
    delta.s <- s * (delta.x - hat.tau) / t
    delta.z <- z * (hat.nu - delta.x) / g
    delta.q <- q * (delta.w - hat.alpha) / p
    delta.v <- v * (gamma.w - delta.w) / w
    delta.p <- p * (gamma.q - delta.q) / q
    delta.g <- g * (gamma.z - delta.z) / z
    delta.t <- t * (gamma.s - delta.s) / s
    step_den <- min(c(
      delta.g / g, delta.w / w, delta.t / t, delta.p / p,
      delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1
    ))
    alfa <- -(1 - margin) / step_den
    newmu <- as.vector(
      crossprod(z, g) + crossprod(v, w) + crossprod(s, t) + crossprod(p, q)
    ) / (2 * (m + n))
    newmu <- mu * ((alfa - 1) / (alfa + 10))^2
    gamma.z <- mu / g - z - delta.z * delta.g / g
    gamma.w <- mu / v - w - delta.w * delta.v / v
    gamma.s <- mu / t - s - delta.s * delta.t / t
    gamma.q <- mu / p - q - delta.q * delta.p / p
    hat.beta <- beta - v * gamma.w / w
    hat.alpha <- alpha - p * gamma.q / q
    hat.nu <- nu + g * gamma.z / z
    hat.tau <- tau - t * gamma.s / s
    c.x <- sigma - z * hat.nu / g - s * hat.tau / t
    c.y <- rho - e * (hat.beta - q * hat.alpha / p)
    AP[xp, xp] <- -H.x
    AP[!xp, !xp] <- H.y
    s1.tmp <- solve(AP, c(c.x, c.y), tol = .Machine$double.eps)
    delta.x <- s1.tmp[seq_len(n)]
    delta.y <- s1.tmp[-seq_len(n)]
    delta.w <- -e * (hat.beta - q * hat.alpha / p + delta.y)
    delta.s <- s * (delta.x - hat.tau) / t
    delta.z <- z * (hat.nu - delta.x) / g
    delta.q <- q * (delta.w - hat.alpha) / p
    delta.v <- v * (gamma.w - delta.w) / w
    delta.p <- p * (gamma.q - delta.q) / q
    delta.g <- g * (gamma.z - delta.z) / z
    delta.t <- t * (gamma.s - delta.s) / s
    step_den <- min(c(
      delta.g / g, delta.w / w, delta.t / t, delta.p / p,
      delta.z / z, delta.v / v, delta.s / s, delta.q / q, -1
    ))
    alfa <- -(1 - margin) / step_den
    x <- x + delta.x * alfa
    g <- g + delta.g * alfa
    w <- w + delta.w * alfa
    t <- t + delta.t * alfa
    p <- p + delta.p * alfa
    y <- y + delta.y * alfa
    z <- z + delta.z * alfa
    v <- v + delta.v * alfa
    s <- s + delta.s * alfa
    q <- q + delta.q * alfa
    mu <- newmu
  }
  as.numeric(x)
}

#' osqp fallback for the same simplex QP
#' @noRd
sc_osqp_w <- function(P, q, J, margin, maxiter, sigf, bound) {
  eps <- 10^(-sigf)
  slack <- margin * eps
  if (!is.finite(slack) || slack <= 0) {
    slack <- .Machine$double.eps
  }
  A <- rbind(diag(J), rep(1, J))
  l <- c(rep(0, J), 1 - slack)
  u <- c(rep(bound, J), 1 + slack)
  res <- osqp::solve_osqp(
    P,
    q,
    A,
    l,
    u,
    pars = osqp::osqpSettings(
      verbose = FALSE,
      eps_abs = eps,
      eps_rel = eps,
      max_iter = as.integer(maxiter),
      eps_prim_inf = slack
    )
  )
  as.numeric(res$x)
}

#' Solve donor weights given diagonal predictor weights
#'
#' Minimizes \eqn{(X_1 - X_0 w)^\top V (X_1 - X_0 w)} subject to
#' \eqn{w \ge 0} and \eqn{1^\top w = 1}.
#'
#' @param X0 \(K \times J\) matrix of donor predictors (typically scaled).
#' @param X1 Length-\(K\) treated predictor vector.
#' @param V Length-\(K\) weights or a \(K \times K\) diagonal matrix.
#' @param margin Extra primal slack. Passed to ipop as the barrier
#'   margin and to the osqp fallback as `eps_prim_inf` / a relaxed
#'   equality tolerance of `margin * 10^(-sigf)`.
#' @param maxiter Iteration limit (ipop Newton steps; osqp iterations).
#' @param sigf Significant figures. ipop stops at this duality-gap
#'   precision (at least 12, matching the Stata plugin on smoking);
#'   weights below `10^(-sigf)` are zeroed after the solve.
#' @param bound Upper bound on each \(w_j\). The default 10 is inert
#'   on the simplex; ipop also uses it as the LOQO clipping bound.
#'
#' @return Named numeric vector of length \(J\) (names from `colnames(X0)`).
#' @noRd
sc_solve_w <- function(X0,
                       X1,
                       V,
                       margin = 0.05,
                       maxiter = 1000,
                       sigf = 7,
                       bound = 10) {
  if (!is.matrix(X0)) {
    X0 <- as.matrix(X0)
  }
  K <- nrow(X0)
  J <- ncol(X0)
  if (is.null(J) || J < 1L) {
    synthaio_stop(
      "synthaio_single_donor",
      "need at least 1 donor column in X0, got ", J
    )
  }

  X1 <- as.numeric(X1)
  if (length(X1) != K) {
    stop("`X1` length must equal nrow(X0).", call. = FALSE)
  }

  if (is.matrix(V)) {
    if (nrow(V) != K || ncol(V) != K) {
      synthaio_stop(
        "synthaio_bad_v",
        "`V` must be length ", K, " or a ", K, " by ", K, " matrix"
      )
    }
    VX0 <- V %*% X0
    VX1 <- as.numeric(V %*% X1)
  } else {
    V <- as.numeric(V)
    if (length(V) != K) {
      synthaio_stop(
        "synthaio_bad_v",
        "`V` must be length ", K, " or a ", K, " by ", K, " matrix"
      )
    }
    VX0 <- V * X0
    VX1 <- V * X1
  }

  P <- crossprod(X0, VX0)
  P <- (P + t(P)) / 2
  q <- -as.numeric(crossprod(X0, VX1))

  sigf <- as.integer(sigf)
  maxiter <- max(1L, as.integer(maxiter))
  # Stata's compiled ipop at sigf=7 already lands on the 5-sparse
  # smoking W; the R port needs a slightly tighter duality gap.
  ipop_sigf <- max(sigf, 12L)
  w_upper <- min(as.numeric(bound), 1)

  w <- tryCatch(
    sc_ipop(
      c = q,
      H = P,
      A = matrix(1, 1L, J),
      b = 1,
      l = rep(0, J),
      u = rep(w_upper, J),
      r = 0,
      sigf = ipop_sigf,
      maxiter = maxiter,
      margin = margin,
      bound = bound
    ),
    error = function(e) NULL
  )
  if (is.null(w) || length(w) != J || !all(is.finite(w))) {
    w <- sc_osqp_w(P, q, J, margin, maxiter, sigf, bound)
  }

  eps <- 10^(-sigf)
  w[w < eps] <- 0
  remainder <- sum(w)
  if (remainder > 0) {
    w <- w / remainder
  }
  names(w) <- colnames(X0)
  w
}
