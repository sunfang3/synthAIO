# 4-unit / 3-donor toys live in this file only (not package data).

plot_panel <- function() {
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

plot_fit <- function(...) {
  scm(
    formula = y ~ x + y(4),
    data = plot_panel(),
    unit = "id",
    time = "t",
    treat = 1,
    trperiod = 5,
    xperiod = 1:4,
    ...
  )
}

test_that("always-available types return a ggplot", {
  fit <- plot_fit()
  for (type in c("path", "effect", "v", "w")) {
    p <- ggplot2::autoplot(fit, type = type)
    expect_true(ggplot2::is_ggplot(p), info = type)
  }
})

test_that("type='gaps' without placebo is synthaio_plot_missing", {
  fit <- plot_fit()
  expect_null(fit$placebo_space)
  expect_error(
    ggplot2::autoplot(fit, type = "gaps"),
    class = "synthaio_plot_missing",
    regexp = "placebo_space"
  )
})

test_that("type='nope' errors", {
  fit <- plot_fit()
  expect_error(ggplot2::autoplot(fit, type = "nope"))
})

test_that("w plot omits exact zeros", {
  fit <- plot_fit()
  fit$W[[1L]] <- 0
  dropped <- names(fit$W)[1L]
  expect_true(any(fit$W == 0))

  p <- ggplot2::autoplot(fit, type = "w")
  expect_true(ggplot2::is_ggplot(p))
  expect_false(any(p$data$weight == 0))
  expect_false(dropped %in% as.character(p$data$unit))
  expect_setequal(
    as.character(p$data$unit),
    names(fit$W)[as.numeric(fit$W) > 0]
  )
})

test_that("placebo and loo types error when the slot is missing", {
  fit <- plot_fit()
  expect_null(fit$placebo_space)
  expect_null(fit$loo)
  expect_error(
    ggplot2::autoplot(fit, type = "mspe"),
    class = "synthaio_plot_missing",
    regexp = "placebo_space"
  )
  expect_error(
    ggplot2::autoplot(fit, type = "pvalue"),
    class = "synthaio_plot_missing",
    regexp = "placebo_space"
  )
  expect_error(
    ggplot2::autoplot(fit, type = "loo"),
    class = "synthaio_plot_missing",
    regexp = "loo"
  )
})

test_that("optional types return a ggplot when slots are attached", {
  fit <- plot_fit(placebo = list(unit = TRUE), loo = TRUE)
  for (type in c("gaps", "mspe", "pvalue", "loo")) {
    p <- ggplot2::autoplot(fit, type = type)
    expect_true(ggplot2::is_ggplot(p), info = type)
  }
})
