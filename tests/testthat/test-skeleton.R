test_that("package name is synthaio not synth2", {
  dcf <- read.dcf(file.path(find.package("synthaio"), "DESCRIPTION"))
  expect_equal(unname(dcf[1, "Package"]), "synthaio")
  expect_false(grepl("synth2", dcf[1, "Package"], fixed = TRUE))
})

test_that("NAMESPACE does not export synth2", {
  ns_lines <- readLines(file.path(find.package("synthaio"), "NAMESPACE"))
  expect_false(any(grepl("synth2", ns_lines, fixed = TRUE)))
})

test_that("package version is 0.0.0.9000", {
  expect_equal(as.character(utils::packageVersion("synthaio")), "0.0.0.9000")
})
