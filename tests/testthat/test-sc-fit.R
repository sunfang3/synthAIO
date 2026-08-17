prop99_formula <- cigsale ~ lnincome + age15to24 + retprice + beer +
  cigsale(1988) + cigsale(1980) + cigsale(1975)

smoking_spec <- function() {
  sc_spec(
    formula = prop99_formula,
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    xperiod = 1980:1988
  )
}

.smoking_allopt <- new.env(parent = emptyenv())

smoking_allopt_fit <- function() {
  if (is.null(.smoking_allopt$fit)) {
    .smoking_allopt$spec <- smoking_spec()
    .smoking_allopt$fit <- sc_fit(.smoking_allopt$spec, method = "allopt")
  }
  .smoking_allopt$fit
}

test_that("smoking allopt has rmse, att, 7-row balance, ~5 nonzero weights", {
  spec <- smoking_spec()
  fit <- smoking_allopt_fit()
  expect_s3_class(fit, "scm_fit")
  expect_identical(fit$spec, spec)
  expect_identical(fit$method, "allopt")
  expect_gt(fit$rmse, 0)
  expect_lt(fit$att, 0)
  expect_s3_class(fit$balance, "data.frame")
  expect_equal(nrow(fit$balance), 7L)
  expect_equal(
    names(fit$balance),
    c(
      "predictor", "v_weight", "treated", "synthetic", "bias_pct",
      "avg_control", "avg_bias_pct"
    )
  )
  n_nz <- sum(fit$W > 0)
  expect_gte(n_nz, 3L)
  expect_lte(n_nz, 7L)
  expect_equal(unname(fit$y_synth), as.numeric(spec$Y0 %*% fit$W))
  expect_equal(unname(fit$effect), unname(fit$y_treated - fit$y_synth))
  expect_equal(names(fit$effect), names(spec$Y1))
})

test_that("effect in 1988 is near 0 relative to 2000", {
  fit <- smoking_allopt_fit()
  expect_lt(
    abs(unname(fit$effect["1988"])),
    0.5 * abs(unname(fit$effect["2000"]))
  )
})

test_that("non-sc_spec input errors", {
  expect_error(sc_fit(list()), regexp = "sc_spec")
  expect_error(sc_fit(smoking), regexp = "sc_spec")
})

test_that("r2 is in (0, 1] on smoking", {
  fit <- smoking_allopt_fit()
  expect_gt(fit$r2, 0)
  expect_lte(fit$r2, 1)
})
