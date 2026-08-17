# In-space placebo on the gold spec. skip_on_cran: 39 nested fits.
# Do not assert pre_mspe == 3.1668; that JSON block is the other run.

test_that("California MSPE ratio ranks first and unfiltered p is 1/n_units", {
  skip_on_cran()
  gold <- gold_smoking()
  fit <- scm(
    cigsale ~ lnincome + age15to24 + retprice + beer +
      cigsale(1988) + cigsale(1980) + cigsale(1975),
    data = smoking, unit = "state", time = "year",
    treat = gold$treat_id, trperiod = 1989,
    xperiod = 1980:1988, nested = TRUE,
    placebo = list(unit = TRUE)
  )
  plc <- fit$placebo_space
  expect_false(is.null(plc))
  n_units <- nrow(plc$table)
  expect_equal(n_units, length(fit$spec$donors) + 1L)
  expect_equal(plc$table$unit[[1L]], gold$treat_id)
  expect_equal(which.max(plc$table$ratio), 1L)
  expect_equal(plc$p_unfiltered, 1 / n_units, tolerance = 1e-9)
})
