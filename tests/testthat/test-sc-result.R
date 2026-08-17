prop99_formula <- cigsale ~ lnincome + age15to24 + retprice + beer +
  cigsale(1988) + cigsale(1980) + cigsale(1975)

.smoking_reg <- new.env(parent = emptyenv())

smoking_regression_fit <- function() {
  if (is.null(.smoking_reg$fit)) {
    spec <- sc_spec(
      formula = prop99_formula,
      data = smoking,
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989,
      xperiod = 1980:1988
    )
    .smoking_reg$fit <- sc_fit(spec, method = "regression")
  }
  .smoking_reg$fit
}

test_that("print/summary return object invisibly and emit RMSE", {
  fit <- smoking_regression_fit()

  printed <- expect_invisible(print(fit))
  expect_identical(printed, fit)
  print_out <- capture.output(print(fit))
  expect_true(any(grepl("RMSE", print_out, fixed = TRUE)))

  summarised <- expect_invisible(summary(fit))
  expect_identical(summarised, fit)
  summary_out <- capture.output(summary(fit))
  expect_true(any(grepl("RMSE", summary_out, fixed = TRUE)))
})

test_that("summary without placebo does not mention p-values", {
  fit <- smoking_regression_fit()
  expect_null(fit$placebo_space)
  expect_null(fit$placebo_time)
  out <- capture.output(summary(fit))
  expect_false(any(grepl("p-value|p value|pvalue|p_value", out, ignore.case = TRUE)))
})

test_that("v_weights(1) errors", {
  expect_error(v_weights(1), regexp = "scm_fit")
})

test_that("coef(fit) names match donor ids", {
  fit <- smoking_regression_fit()
  expect_identical(coef(fit), fit$W)
  expect_equal(names(coef(fit)), as.character(fit$spec$donors))
  expect_true("34" %in% names(coef(fit)))
})

test_that("capture.output(summary(fit)) contains RMSE on smoking regression fit", {
  fit <- smoking_regression_fit()
  out <- capture.output(summary(fit))
  expect_true(any(grepl("RMSE", out, fixed = TRUE)))
})

test_that("extractors return stored components without recomputing", {
  fit <- smoking_regression_fit()
  expect_identical(v_weights(fit), fit$V)
  expect_identical(balance(fit), fit$balance)
  expect_identical(effects(fit), fit$effect)

  g <- glance(fit)
  expect_s3_class(g, "data.frame")
  expect_equal(nrow(g), 1L)
  expect_equal(g$att, fit$att)
  expect_equal(g$rmse, fit$rmse)
  expect_equal(g$r2, fit$r2)
  expect_equal(g$pre_mspe, fit$pre_mspe)
  expect_equal(g$post_mspe, fit$post_mspe)
  expect_equal(g$mspe_ratio, fit$mspe_ratio)
})
