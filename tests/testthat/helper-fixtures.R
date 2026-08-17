gold_smoking <- function(path = NULL) {
  if (is.null(path)) {
    path <- testthat::test_path("fixtures/smoking-synth2-nested-allopt.json")
  }
  if (!file.exists(path)) {
    stop("gold fixture not found: ", path, call. = FALSE)
  }
  jsonlite::fromJSON(path)
}

.gold_fit_cache <- new.env(parent = emptyenv())

gold_smoking_fit <- function() {
  if (is.null(.gold_fit_cache$fit)) {
    gold <- gold_smoking()
    .gold_fit_cache$fit <- scm(
      cigsale ~ lnincome + age15to24 + retprice + beer +
        cigsale(1988) + cigsale(1980) + cigsale(1975),
      data = smoking, unit = "state", time = "year",
      treat = gold$treat_id, trperiod = 1989,
      xperiod = 1980:1988, allopt = TRUE
    )
  }
  .gold_fit_cache$fit
}

gold_id_to_name <- function(gold, ids) {
  revmap <- stats::setNames(
    names(gold$unit_ids),
    as.character(unlist(gold$unit_ids, use.names = FALSE))
  )
  unname(revmap[as.character(ids)])
}
