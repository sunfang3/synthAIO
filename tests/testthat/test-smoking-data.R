test_that("smoking is the 39-state Prop 99 panel", {
  expect_s3_class(smoking, "data.frame")
  expect_equal(nrow(smoking), 1209L)
  expect_equal(length(unique(smoking$state)), 39L)
  expect_equal(range(as.integer(smoking$year)), c(1970L, 2000L))
  expect_true(all(
    c("state", "year", "cigsale", "lnincome", "age15to24", "retprice", "beer")
    %in% names(smoking)
  ))

  gold <- gold_smoking()
  states <- unique(smoking$state)
  cal_present <-
    "California" %in% as.character(states) ||
    identical(gold$treat_name, "California") && (
      gold$treat_id %in% states ||
        gold$treat_name %in% as.character(states)
    )
  expect_true(cal_present)
})

test_that("gold_smoking errors if the JSON is missing", {
  expect_error(
    gold_smoking(path = "missing-gold-fixture.json"),
    regexp = "(gold|fixture).*(not found|missing)|(not found|missing).*(gold|fixture)"
  )
})

test_that("gold fixture has required schema", {
  gold <- gold_smoking()
  expect_true(all(c("w_nonzero", "v", "tol", "_source") %in% names(gold)))
  expect_identical(gold$treat_name, "California")
  expect_false(is.null(gold$treat_id))
  if (is.numeric(smoking$state)) {
    expect_true(all(names(gold$w_nonzero) %in% names(gold$unit_ids)))
  }
})
