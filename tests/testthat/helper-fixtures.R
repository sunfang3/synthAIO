gold_smoking <- function(path = NULL) {
  if (is.null(path)) {
    path <- testthat::test_path("fixtures/smoking-synth2-nested-allopt.json")
  }
  if (!file.exists(path)) {
    stop("gold fixture not found: ", path, call. = FALSE)
  }
  jsonlite::fromJSON(path)
}
