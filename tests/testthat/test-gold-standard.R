# Gold-standard alignment for smoking nested-allopt.
# Donor set is hard against the synth2 fixture. Individual W uses
# max(tol$w_abs, 0.05*|w|) after smoking-synth2-achieved.json (Q1).
# Do not assert pre_mspe == 3.1668 on the allopt fit.

test_that("gold donor set matches synth2 nested allopt", {
  gold <- gold_smoking()
  fit <- gold_smoking_fit()
  tol <- gold$tol$w_abs
  nz_ids <- names(fit$W)[fit$W > tol]
  got <- sort(gold_id_to_name(gold, nz_ids))
  expect_identical(got, sort(names(gold$w_nonzero)))
})

test_that("nonzero W is within relaxed gold tolerance", {
  gold <- gold_smoking()
  fit <- gold_smoking_fit()
  for (nm in names(gold$w_nonzero)) {
    id <- as.character(gold$unit_ids[[nm]])
    target <- gold$w_nonzero[[nm]]
    tol <- max(gold$tol$w_abs, 0.05 * abs(target))
    expect_lt(abs(fit$W[[id]] - target), tol + 1e-12, label = nm)
  }
})

test_that("RMSE relative error and ATT match gold tolerances", {
  gold <- gold_smoking()
  fit <- gold_smoking_fit()
  expect_lte(abs(fit$rmse - gold$rmse) / gold$rmse, gold$tol$rmse_rel)
  expect_lte(abs(fit$att - gold$att), gold$tol$att_abs)
})

test_that("selected years' effects are within 0.15 packs", {
  gold <- gold_smoking()
  fit <- gold_smoking_fit()
  for (yr in names(gold$effect)) {
    expect_lte(
      abs(unname(fit$effect[[yr]]) - gold$effect[[yr]]),
      0.15,
      label = yr
    )
  }
})

test_that("Stata V recovers Stata W through the inner QP", {
  gold <- gold_smoking()
  spec <- sc_spec(
    cigsale ~ lnincome + age15to24 + retprice + beer +
      cigsale(1988) + cigsale(1980) + cigsale(1975),
    data = smoking, unit = "state", time = "year",
    treat = gold$treat_id, trperiod = 1989,
    xperiod = 1980:1988
  )
  v <- unlist(gold$v[rownames(spec$X0_scaled)])
  out <- synthaio:::sc_solve_v(spec, method = "custom", custom_v = v)
  for (nm in names(gold$w_nonzero)) {
    id <- as.character(gold$unit_ids[[nm]])
    expect_lt(
      abs(out$W[[id]] - gold$w_nonzero[[nm]]),
      gold$tol$w_abs + 1e-12,
      label = nm
    )
  }
  expect_identical(
    sort(gold_id_to_name(gold, names(out$W)[out$W > gold$tol$w_abs])),
    sort(names(gold$w_nonzero))
  )
})

test_that("V top-2 matches the Stata nested-allopt basin", {
  gold <- gold_smoking()
  fit <- gold_smoking_fit()
  got_top <- names(sort(fit$V, decreasing = TRUE))[1:2]
  gold_top <- names(sort(unlist(gold$v), decreasing = TRUE))[1:2]
  expect_identical(sort(got_top), sort(gold_top))
})
