# 3-unit toys live in this file only (not package data).

jump_panel <- function() {
  data.frame(
    id = rep(1:3, each = 6L),
    t = rep(1:6, 3L),
    y = c(
      10.0, 11.0, 12.0, 13.0, 30.0, 31.0,
      10.2, 10.8, 12.1, 12.9, 14.0, 15.1,
      9.8, 11.2, 11.9, 13.1, 13.8, 14.9
    ),
    x = c(
      1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
      1.02, 1.01, 0.99, 1.00, 1.00, 1.00,
      0.98, 0.99, 1.02, 1.01, 1.00, 1.00
    )
  )
}

jump_spec <- function(trperiod = 5, xperiod = 1:4, formula = y ~ x + y(4)) {
  sc_spec(
    formula = formula,
    data = jump_panel(),
    unit = "id",
    time = "t",
    treat = 1,
    trperiod = trperiod,
    xperiod = xperiod
  )
}

wild_panel <- function() {
  data.frame(
    id = rep(1:3, each = 6L),
    t = rep(1:6, 3L),
    y = c(
      10, 11, 12, 13, 30, 31,
      10, 11, 12, 13, 14, 15,
      50, 1, 50, 1, 20, 20
    ),
    x = c(
      1, 1, 1, 1, 1, 1,
      1, 1, 1, 1, 1, 1,
      8, 0, 8, 0, 3, 3
    )
  )
}

wild_spec <- function() {
  sc_spec(
    y ~ x + y,
    data = wild_panel(),
    unit = "id",
    time = "t",
    treat = 1,
    trperiod = 5,
    xperiod = 1:4
  )
}

test_that("3-unit toy: treated has the unique post jump and the largest ratio", {
  spec <- jump_spec()
  fit <- sc_fit(spec, method = "regression")
  out <- synthaio:::sc_placebo_space(spec, fit)

  expect_s3_class(out$table, "data.frame")
  expect_equal(
    names(out$table),
    c("unit", "pre_mspe", "post_mspe", "ratio", "kept")
  )
  expect_equal(sort(out$table$unit), 1:3)
  expect_true(all(out$table$kept))

  treated <- out$table[out$table$unit == 1L, ]
  others <- out$table[out$table$unit != 1L, ]
  expect_gt(treated$ratio, max(others$ratio))
})

test_that("unfiltered Fisher p is 1/3 when treated has the largest ratio", {
  spec <- jump_spec()
  fit <- sc_fit(spec, method = "regression")
  out <- synthaio:::sc_placebo_space(spec, fit)

  expect_equal(out$p_unfiltered, 1 / 3)
  expect_equal(out$p_filtered, 1 / 3)
  n_kept <- sum(out$table$kept)
  expect_true(all(out$p_left >= 1 / n_kept, na.rm = TRUE))
  expect_true(all(out$p_two >= 1 / n_kept, na.rm = TRUE))
})

test_that("cut that drops every placebo yields NA year-wise p and synthaio_cut_empty", {
  spec <- jump_spec()
  fit <- sc_fit(spec, method = "regression")
  unf <- synthaio:::sc_placebo_space(spec, fit)
  tr_pre <- unf$table$pre_mspe[unf$table$unit == spec$treat]
  pl_pre <- unf$table$pre_mspe[unf$table$unit != spec$treat]
  cut <- min(pl_pre) / max(tr_pre, .Machine$double.eps) / 2

  out <- NULL
  expect_warning(
    out <- synthaio:::sc_placebo_space(spec, fit, cut = cut),
    class = "synthaio_cut_empty"
  )
  expect_false(any(out$table$kept[out$table$unit != spec$treat]))
  expect_true(all(is.na(out$p_left)))
  expect_true(all(is.na(out$p_two)))
})

test_that("in-time period >= trperiod errors", {
  spec <- jump_spec()
  expect_error(synthaio:::sc_placebo_time(spec, period = spec$trperiod))
  expect_error(synthaio:::sc_placebo_time(spec, period = spec$trperiod + 1))
})

test_that("cut = 2 excludes a placebo with huge pre MSPE", {
  spec <- wild_spec()
  fit <- sc_fit(spec, method = "regression")
  out <- synthaio:::sc_placebo_space(spec, fit, cut = 2)

  wild <- out$table[out$table$unit == 3L, ]
  expect_false(wild$kept)
  expect_gt(wild$pre_mspe, 2 * fit$pre_mspe)
  expect_true(all(out$table$kept[out$table$unit != 3L]))
  expect_equal(sum(out$table$kept), 2L)
})

test_that("mixed args call in-time then in-space on the shifted spec", {
  spec <- jump_spec()
  fit <- sc_fit(spec, method = "regression")
  t0 <- 4

  mixed <- synthaio:::sc_placebo(
    spec, fit, unit = TRUE, period = t0, cut = Inf
  )
  time_res <- synthaio:::sc_placebo_time(spec, period = t0, method = fit$method)
  space_res <- synthaio:::sc_placebo_space(time_res$spec, time_res$fit, cut = Inf)
  orig_space <- synthaio:::sc_placebo_space(spec, fit, cut = Inf)

  expect_equal(time_res$spec$trperiod, t0)
  expect_equal(mixed$placebo_time$spec$trperiod, t0)
  expect_equal(mixed$placebo_space$table, space_res$table)
  expect_equal(mixed$placebo_time$fit$effect, time_res$fit$effect)
  expect_false(isTRUE(all.equal(mixed$placebo_space$table$ratio, orig_space$table$ratio)))
})

test_that("unit = c(2) runs only that donor as a placebo", {
  spec <- jump_spec()
  fit <- sc_fit(spec, method = "regression")
  out <- synthaio:::sc_placebo_space(spec, fit, unit = 2)

  expect_equal(sort(out$table$unit), c(1, 2))
  expect_false(3 %in% out$table$unit)
  expect_equal(out$p_unfiltered, sum(out$table$ratio >= out$table$ratio[1]) / 2)
})

test_that("unknown placebo unit id errors", {
  spec <- jump_spec()
  fit <- sc_fit(spec, method = "regression")
  expect_error(
    synthaio:::sc_placebo_space(spec, fit, unit = 99),
    class = "synthaio_bad_placebo_unit"
  )
})

test_that("in-time clips explicit mspeperiod and preperiod to times < t0", {
  spec <- sc_spec(
    y ~ x + y(3),
    data = jump_panel(),
    unit = "id",
    time = "t",
    treat = 1,
    trperiod = 5,
    xperiod = 1:4,
    mspeperiod = 1:4,
    preperiod = 1:4,
    postperiod = 5:6
  )
  expect_equal(spec$mspe_times, 1:4)
  expect_equal(spec$pre_times, 1:4)

  out <- synthaio:::sc_placebo_time(spec, period = 4, method = "regression")
  expect_equal(out$spec$mspe_times, 1:3)
  expect_equal(out$spec$pre_times, 1:3)
  expect_true(all(out$spec$post_times >= 4))
  expect_false(any(out$spec$mspe_times >= 4))
})

test_that("in-time t0 drops y(t >= t0) and re-averages bare x on xperiod < t0", {
  xperiod <- 2:4
  t0 <- 4
  spec <- jump_spec(trperiod = 5, xperiod = xperiod, formula = y ~ x + y(4))
  expect_true("y_4" %in% names(spec$X1))
  expect_equal(spec$x_times, 2:4)

  out <- synthaio:::sc_placebo_time(spec, period = t0, method = "regression")
  shifted <- out$spec

  expect_equal(shifted$trperiod, t0)
  expect_false("y_4" %in% names(shifted$X1))
  expect_equal(names(shifted$X1), "x")
  expect_equal(shifted$x_times, 2:3)
  expect_equal(shifted$xperiod, 2:3)

  panel <- jump_panel()
  expect_equal(
    unname(shifted$X1[["x"]]),
    mean(panel$x[panel$id == 1L & panel$t %in% 2:3], na.rm = TRUE)
  )
  expect_true(all(as.numeric(names(out$effect_fake)) >= t0))
  expect_true(all(as.numeric(names(out$effect_fake)) < spec$trperiod))
  expect_equal(as.numeric(names(out$effect_post)), spec$post_times)
})
