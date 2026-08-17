prop99_formula <- cigsale ~ lnincome + age15to24 + retprice + beer +
  cigsale(1988) + cigsale(1980) + cigsale(1975)

toy_spec <- function() {
  panel <- data.frame(
    state = rep(c(1, 2, 3), each = 3),
    year = rep(2000:2002, 3),
    y = c(1, 2, 5, 1, 2, 3, 1.5, 2.5, 4),
    x = c(1, 1, 1, 1, 1, 1, 2, 2, 2)
  )
  sc_spec(
    y ~ x + y(2001),
    data = panel,
    unit = "state",
    time = "year",
    treat = 1,
    trperiod = 2002
  )
}

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

test_that("custom equal weights on a 2-donor toy spec is a no-op normalize", {
  spec <- toy_spec()
  custom_v <- c(0.5, 0.5)
  out <- synthaio:::sc_solve_v(spec, method = "custom", custom_v = custom_v)
  expect_equal(unname(out$V), c(0.5, 0.5))
  expect_equal(sum(out$V), 1)
  expect_equal(names(out$V), names(spec$X1))
  expect_length(out$W, 2L)
  expect_equal(names(out$W), colnames(spec$X0))
  expect_equal(sum(out$W), 1)
  expect_true(all(out$W >= 0))
})

test_that("custom with custom_v = NULL errors", {
  spec <- toy_spec()
  expect_error(
    synthaio:::sc_solve_v(spec, method = "custom", custom_v = NULL)
  )
})

test_that("custom_v with negatives errors", {
  spec <- toy_spec()
  expect_error(
    synthaio:::sc_solve_v(spec, method = "custom", custom_v = c(1, -0.2))
  )
})

test_that("nested MSPE is at most regression MSPE on smoking", {
  spec <- smoking_spec()
  reg <- synthaio:::sc_solve_v(spec, method = "regression")
  nes <- synthaio:::sc_solve_v(spec, method = "nested")
  expect_lte(nes$mspe, reg$mspe)
  expect_equal(sum(nes$V), 1)
  expect_equal(sum(reg$V), 1)
})

test_that("allopt MSPE is at most nested MSPE on smoking", {
  spec <- smoking_spec()
  nes <- synthaio:::sc_solve_v(spec, method = "nested")
  ao <- synthaio:::sc_solve_v(spec, method = "allopt", start_seed = 1)
  expect_lte(ao$mspe, nes$mspe)
})

test_that("allopt is deterministic given start_seed", {
  spec <- smoking_spec()
  a <- synthaio:::sc_solve_v(spec, method = "allopt", start_seed = 1)
  b <- synthaio:::sc_solve_v(spec, method = "allopt", start_seed = 1)
  expect_equal(a$V, b$V)
  expect_equal(a$W, b$W)
  expect_equal(a$mspe, b$mspe)
})

test_that("smoking allopt donor-set overlap is at least 3 of gold donors", {
  spec <- smoking_spec()
  ao <- synthaio:::sc_solve_v(spec, method = "allopt", start_seed = 1)
  gold_ids <- c("4", "5", "19", "21", "34")
  overlap <- gold_ids[ao$W[gold_ids] > 0]
  expect_gte(length(overlap), 3L)
})
