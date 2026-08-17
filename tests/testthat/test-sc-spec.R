prop99_formula <- cigsale ~ lnincome + age15to24 + retprice + beer +
  cigsale(1988) + cigsale(1980) + cigsale(1975)

test_that("treat = 'California' resolves via state_name", {
  spec <- sc_spec(
    formula = prop99_formula,
    data = smoking,
    unit = "state",
    time = "year",
    treat = "California",
    trperiod = 1989,
    xperiod = 1980:1988
  )
  expect_equal(spec$treat, 3)
  expect_equal(sc_format_unit(spec, 3), "California (3)")
  expect_equal(sc_format_unit(spec, c(3, 34)), c("California (3)", "Utah (34)"))
})

test_that("smoking Prop 99 spec has K=7, J=38, cigsale_1988 ≈ 90.1", {
  spec <- sc_spec(
    formula = prop99_formula,
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    xperiod = 1980:1988
  )

  expect_s3_class(spec, "sc_spec")
  expect_equal(length(spec$X1), 7L)
  expect_equal(dim(spec$X0), c(7L, 38L))
  expect_equal(nrow(spec$Y0), 31L)
  expect_equal(length(spec$Y1), 31L)
  expect_equal(ncol(spec$Y0), 38L)
  expect_equal(unname(spec$X1[["cigsale_1988"]]), 90.1, tolerance = 0.05)
  expect_equal(
    names(spec$X1),
    c(
      "lnincome", "age15to24", "retprice", "beer",
      "cigsale_1988", "cigsale_1980", "cigsale_1975"
    )
  )
  expect_equal(rownames(spec$X0), names(spec$X1))
  expect_equal(spec$outcome, "cigsale")
  expect_equal(spec$treat, 3)
  expect_equal(spec$trperiod, 1989)
  expect_equal(spec$times, 1970:2000)
  expect_equal(spec$pre_times, 1970:1988)
  expect_equal(spec$post_times, 1989:2000)
  expect_equal(length(spec$donors), 38L)
  expect_false(3 %in% spec$donors)

  expect_equal(spec$formula, prop99_formula)
  expect_equal(spec$data, smoking)
  expect_equal(spec$unit, "state")
  expect_equal(spec$time, "year")
  expect_equal(spec$xperiod, 1980:1988)
  expect_null(spec$mspeperiod)
  expect_null(spec$counit)
  expect_null(spec$preperiod)
  expect_null(spec$postperiod)

  X_all <- cbind(spec$X0, spec$X1)
  sds <- apply(X_all, 1L, stats::sd)
  expect_equal(dim(spec$X0_scaled), dim(spec$X0))
  expect_equal(length(spec$X1_scaled), length(spec$X1))
  expect_equal(as.numeric(spec$X1_scaled), as.numeric(spec$X1 / sds))
  expect_equal(as.matrix(spec$X0_scaled), as.matrix(spec$X0 / sds))
})

test_that("0-row data is synthaio_not_panel", {
  expect_error(
    sc_spec(
      cigsale ~ lnincome,
      data = smoking[0, ],
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989
    ),
    class = "synthaio_not_panel"
  )
})

test_that("duplicate unit-time is synthaio_not_panel", {
  dup <- rbind(smoking[1, ], smoking)
  expect_error(
    sc_spec(
      cigsale ~ lnincome,
      data = dup,
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989
    ),
    class = "synthaio_not_panel"
  )
})

test_that("predictor time absent from panel is synthaio_bad_time", {
  expect_error(
    sc_spec(
      cigsale ~ cigsale(1960),
      data = smoking,
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989
    ),
    class = "synthaio_bad_time"
  )
})

test_that("treat id missing from the panel is synthaio_missing_treat", {
  expect_error(
    sc_spec(
      cigsale ~ lnincome,
      data = smoking,
      unit = "state",
      time = "year",
      treat = 99,
      trperiod = 1989
    ),
    class = "synthaio_missing_treat"
  )
})

test_that("fewer than two donors is synthaio_single_donor", {
  expect_error(
    sc_spec(
      cigsale ~ lnincome,
      data = smoking,
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989,
      counit = 1
    ),
    class = "synthaio_single_donor"
  )
})

test_that("missing required unit-time is synthaio_unbalanced", {
  unbalanced <- smoking[!(smoking$state == 1L & smoking$year == 1988L), ]
  expect_error(
    sc_spec(
      cigsale ~ cigsale(1988),
      data = unbalanced,
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989
    ),
    class = "synthaio_unbalanced"
  )
})

test_that("counit restricts donors and drops treated with a warning", {
  spec <- sc_spec(
    cigsale ~ lnincome,
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    counit = c(1, 2, 4, 5)
  )
  expect_equal(sort(spec$donors), c(1, 2, 4, 5))
  expect_equal(ncol(spec$X0), 4L)
  expect_equal(spec$counit, c(1, 2, 4, 5))

  spec2 <- NULL
  expect_warning(
    spec2 <- sc_spec(
      cigsale ~ lnincome,
      data = smoking,
      unit = "state",
      time = "year",
      treat = 3,
      trperiod = 1989,
      counit = c(3, 1, 2, 4)
    ),
    regexp = "counit"
  )
  expect_false(3 %in% spec2$donors)
  expect_equal(sort(spec2$donors), c(1, 2, 4))
  expect_equal(ncol(spec2$X0), 3L)
})

test_that("default xperiod is all times < trperiod", {
  spec <- sc_spec(
    cigsale ~ lnincome,
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989
  )
  expect_equal(spec$x_times, 1970:1988)
  expect_null(spec$xperiod)
  expect_equal(spec$pre_times, 1970:1988)
  expect_equal(spec$post_times, 1989:2000)
})

test_that("mspeperiod default is independent of xperiod", {
  spec <- sc_spec(
    cigsale ~ lnincome,
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    xperiod = 1980:1988,
    mspeperiod = NULL
  )
  expect_equal(spec$x_times, 1980:1988)
  expect_equal(spec$mspe_times, 1970:1988)
  expect_equal(spec$xperiod, 1980:1988)
  expect_null(spec$mspeperiod)
})

test_that("formula times name predictors and override xperiod", {
  spec <- sc_spec(
    cigsale ~ lnincome + cigsale(1988) + cigsale(1980:1988) +
      cigsale(c(1982, 1986, 1988)),
    data = smoking,
    unit = "state",
    time = "year",
    treat = 3,
    trperiod = 1989,
    xperiod = 1975:1979
  )
  expect_equal(
    names(spec$X1),
    c("lnincome", "cigsale_1988", "cigsale_1980_1988", "cigsale_1982_1986_1988")
  )
  cal <- smoking[smoking$state == 3L, ]
  expect_equal(
    unname(spec$X1[["lnincome"]]),
    mean(cal$lnincome[cal$year %in% 1975:1979], na.rm = TRUE)
  )
  expect_equal(unname(spec$X1[["cigsale_1988"]]), 90.1, tolerance = 0.05)
  expect_equal(
    unname(spec$X1[["cigsale_1980_1988"]]),
    mean(cal$cigsale[cal$year %in% 1980:1988], na.rm = TRUE)
  )
  expect_equal(
    unname(spec$X1[["cigsale_1982_1986_1988"]]),
    mean(cal$cigsale[cal$year %in% c(1982, 1986, 1988)], na.rm = TRUE)
  )
})
