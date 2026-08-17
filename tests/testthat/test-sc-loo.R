# 4-unit / 3-donor toys live in this file only (not package data).

loo_panel <- function() {
  data.frame(
    id = rep(1:4, each = 6L),
    t = rep(1:6, 4L),
    y = c(
      10.0, 11.0, 12.0, 13.0, 30.0, 31.0,
      11.0, 12.0, 13.0, 14.0, 15.0, 16.0,
      9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
      50.0, 1.0, 50.0, 1.0, 20.0, 20.0
    ),
    x = c(
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.10, 1.10, 1.10, 1.10, 1.10, 1.10,
      0.90, 0.90, 0.90, 0.90, 0.90, 0.90,
      8.00, 0.00, 8.00, 0.00, 3.00, 3.00
    )
  )
}

loo_spec <- function() {
  sc_spec(
    y ~ x + y(4),
    data = loo_panel(),
    unit = "id",
    time = "t",
    treat = 1,
    trperiod = 5,
    xperiod = 1:4
  )
}

test_that("3-donor toy LOO drops each nonzero donor; band contains original effect", {
  spec <- loo_spec()
  fit <- sc_fit(spec, method = "regression")
  nz <- names(fit$W)[fit$W > 10^(-7)]
  expect_gte(length(nz), 2L)
  expect_lte(length(nz), 3L)

  out <- synthaio:::sc_loo(spec, fit)

  expect_s3_class(out$band, "data.frame")
  expect_equal(names(out$band), c("time", "effect", "loo_min", "loo_max"))
  expect_equal(out$band$time, spec$times)
  expect_equal(out$band$effect, unname(fit$effect[as.character(spec$times)]))
  expect_length(out$fits, length(nz))
  expect_equal(sort(names(out$fits)), sort(nz))
  expect_true(all(out$band$loo_min <= out$band$effect))
  expect_true(all(out$band$effect <= out$band$loo_max))
})

test_that("W all below threshold errors with synthaio_loo_empty", {
  spec <- loo_spec()
  fit <- sc_fit(spec, method = "regression")
  fit$W[] <- 0

  expect_error(
    synthaio:::sc_loo(spec, fit),
    class = "synthaio_loo_empty"
  )
})

test_that("non-scm_fit input errors", {
  spec <- loo_spec()
  expect_error(synthaio:::sc_loo(spec, list()), regexp = "scm_fit")
  expect_error(synthaio:::sc_loo(spec, spec), regexp = "scm_fit")
})

test_that("dropping the largest donor still produces sum(W) = 1", {
  spec <- loo_spec()
  fit <- sc_fit(spec, method = "regression")
  out <- synthaio:::sc_loo(spec, fit)

  largest <- names(fit$W)[which.max(fit$W)]
  w_loo <- if (inherits(out$fits[[largest]], "scm_fit")) {
    out$fits[[largest]]$W
  } else {
    out$fits[[largest]]
  }
  expect_equal(sum(w_loo), 1)
})
