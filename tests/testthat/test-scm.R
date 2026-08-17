prop99_formula <- cigsale ~ lnincome + age15to24 + retprice + beer +
  cigsale(1988) + cigsale(1980) + cigsale(1975)

# Tiny 4-unit / 3-donor toys live in this file only (not package data).
scm_toy_panel <- function() {
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

scm_toy_call <- function(...) {
  dots <- list(...)
  args <- list(
    formula = y ~ x + y(4),
    data = scm_toy_panel(),
    unit = "id",
    time = "t",
    treat = 1,
    trperiod = 5,
    xperiod = 1:4
  )
  args[names(dots)] <- dots
  do.call(scm, args)
}

test_that("scm smoking nested=FALSE returns scm_fit with numeric att", {
  fit <- scm(
    formula = prop99_formula,
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    xperiod = 1980:1988,
    nested = FALSE
  )

  expect_s3_class(fit, "scm_fit")
  expect_true(is.numeric(fit$att))
  expect_length(fit$att, 1L)
  expect_false(is.na(fit$att))
  expect_identical(fit$method, "regression")
  expect_equal(length(fit$spec$donors), 38L)
  expect_equal(fit$spec$treat, 3)
  expect_equal(fit$spec$trperiod, 1989)
})

test_that("loo=TRUE on a tiny toy attaches $loo", {
  fit <- scm_toy_call(loo = TRUE)

  expect_s3_class(fit, "scm_fit")
  expect_false(is.null(fit$loo))
  expect_s3_class(fit$loo$band, "data.frame")
  expect_equal(names(fit$loo$band), c("time", "effect", "loo_min", "loo_max"))
  expect_true(length(fit$loo$fits) >= 1L)
})

test_that("placebo unit ids are passed through scm()", {
  fit <- scm_toy_call(placebo = list(unit = 2))

  expect_false(is.null(fit$placebo_space))
  expect_equal(sort(fit$placebo_space$table$unit), c(1, 2))
  expect_false(3 %in% fit$placebo_space$table$unit)
  expect_false(4 %in% fit$placebo_space$table$unit)
})

test_that("placebo=list() is treated as no placebo", {
  fit <- scm_toy_call(placebo = list())

  expect_s3_class(fit, "scm_fit")
  expect_null(fit$placebo_space)
  expect_null(fit$placebo_time)
})

test_that("missing trperiod errors", {
  expect_error(
    scm(
      formula = prop99_formula,
      data = smoking,
      unit = "state",
      time = "year",
      treat = 3
    ),
    regexp = "trperiod"
  )
})

test_that("allopt=TRUE with nested=FALSE warns and runs allopt", {
  fit <- NULL
  expect_warning(
    fit <- scm_toy_call(allopt = TRUE, nested = FALSE),
    regexp = "allopt"
  )
  expect_s3_class(fit, "scm_fit")
  expect_identical(fit$method, "allopt")
  expect_true(is.numeric(fit$att))
})
