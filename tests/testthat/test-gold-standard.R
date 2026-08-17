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

test_that("V top-2 and weights lock to the achieved basin", {
  achieved <- gold_smoking_achieved()
  fit <- gold_smoking_fit()
  got_top <- names(sort(fit$V, decreasing = TRUE))[1:2]
  expect_identical(sort(got_top), sort(achieved$v_top2))
  # Gold V (age15to24 + cigsale_1975) is a different local min.
  expect_identical(sort(achieved$gold_v_top2), c("age15to24", "cigsale_1975"))
  for (nm in names(achieved$v)) {
    expect_equal(unname(fit$V[[nm]]), achieved$v[[nm]], tolerance = 1e-4)
  }
})
