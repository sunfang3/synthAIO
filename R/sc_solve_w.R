#' Solve donor weights given diagonal predictor weights
#'
#' Minimizes \eqn{(X_1 - X_0 w)^\top V (X_1 - X_0 w)} subject to
#' \eqn{w \ge 0} and \eqn{1^\top w = 1}.
#'
#' @param X0 \(K \times J\) matrix of donor predictors (typically scaled).
#' @param X1 Length-\(K\) treated predictor vector.
#' @param V Length-\(K\) weights or a \(K \times K\) diagonal matrix.
#' @param margin Extra primal slack. Passed to osqp as `eps_prim_inf` and
#'   as a relaxed equality tolerance of `margin * 10^(-sigf)`.
#' @param maxiter osqp iteration limit.
#' @param sigf Significant figures. Solver `eps_abs` / `eps_rel` are
#'   `10^(-sigf)`; weights below that threshold are zeroed after the solve.
#' @param bound Upper bound on each \(w_j\). The default 10 is inert because
#'   \(w_j \le 1\) on the simplex.
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

  eps <- 10^(-sigf)
  slack <- margin * eps
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

  w <- as.numeric(res$x)
  w[w < eps] <- 0
  remainder <- sum(w)
  if (remainder > 0) {
    w <- w / remainder
  }
  names(w) <- colnames(X0)
  w
}
