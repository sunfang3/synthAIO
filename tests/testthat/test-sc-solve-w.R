test_that("X1 equal to donor 1 puts all weight on donor 1", {
  X0 <- cbind(d1 = c(1, 2), d2 = c(5, 6))
  rownames(X0) <- c("a", "b")
  X1 <- c(1, 2)
  w <- synthaio:::sc_solve_w(X0, X1, V = c(1, 1))
  expect_equal(unname(w), c(1, 0), tolerance = 1e-6)
  expect_equal(names(w), c("d1", "d2"))
})

test_that("X1 midpoint of two donors splits weight equally", {
  X0 <- cbind(d1 = c(0, 0), d2 = c(2, 2))
  rownames(X0) <- c("a", "b")
  X1 <- c(1, 1)
  w <- synthaio:::sc_solve_w(X0, X1, V = c(1, 1))
  expect_equal(unname(w), c(0.5, 0.5), tolerance = 1e-6)
})

test_that("X0 with 0 columns is synthaio_single_donor", {
  X0 <- matrix(numeric(), nrow = 2L, ncol = 0L)
  expect_error(
    synthaio:::sc_solve_w(X0, X1 = c(1, 2), V = c(1, 1)),
    class = "synthaio_single_donor"
  )
})

test_that("V length not equal to nrow(X0) errors", {
  X0 <- cbind(d1 = c(1, 2), d2 = c(3, 4))
  X1 <- c(1, 2)
  expect_error(
    synthaio:::sc_solve_w(X0, X1, V = 1),
    class = "synthaio_bad_v"
  )
  expect_error(
    synthaio:::sc_solve_w(X0, X1, V = diag(3)),
    class = "synthaio_bad_v"
  )
})

test_that("weights lie on the simplex", {
  X0 <- cbind(d1 = c(1, 0, 2), d2 = c(0, 1, 3), d3 = c(2, 2, 1))
  X1 <- c(1, 1, 2)
  sigf <- 7
  w <- synthaio:::sc_solve_w(X0, X1, V = c(1, 1, 1), sigf = sigf)
  expect_equal(sum(w), 1)
  expect_true(all(w >= -10^(-sigf)))
})

test_that("custom V that zeros a predictor does not use that row", {
  X0 <- cbind(d1 = c(1, 100), d2 = c(9, 2))
  X1 <- c(1, 2)
  w_first <- synthaio:::sc_solve_w(X0, X1, V = c(1, 0))
  expect_equal(unname(w_first), c(1, 0), tolerance = 1e-5)
  w_second <- synthaio:::sc_solve_w(X0, X1, V = c(0, 1))
  expect_equal(unname(w_second), c(0, 1), tolerance = 1e-5)
  w_diag <- synthaio:::sc_solve_w(X0, X1, V = diag(c(1, 0)))
  expect_equal(unname(w_diag), c(1, 0), tolerance = 1e-5)
})

test_that("smoking Prop 99 equal V returns a length-38 simplex", {
  spec <- sc_spec(
    cigsale ~ lnincome + age15to24 + retprice + beer +
      cigsale(1988) + cigsale(1980) + cigsale(1975),
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    xperiod = 1980:1988
  )
  K <- nrow(spec$X0_scaled)
  V <- rep(1 / K, K)
  w <- synthaio:::sc_solve_w(spec$X0_scaled, spec$X1_scaled, V)
  expect_length(w, 38L)
  expect_equal(sum(w), 1)
  expect_true(all(w >= -1e-7))
  expect_equal(names(w), colnames(spec$X0_scaled))
})
