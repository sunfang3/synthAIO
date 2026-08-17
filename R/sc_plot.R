SC_PLOT_TREATED <- "#B45C3D"
SC_PLOT_CONTROL <- "#8A8A8A"

#' Require an optional `scm_fit` slot for a plot type
#' @noRd
sc_plot_require <- function(object, slot) {
  val <- object[[slot]]
  if (is.null(val)) {
    synthaio_stop(
      "synthaio_plot_missing",
      "plot type requires `", slot, "`"
    )
  }
  val
}

#' Outcome path as a two-column time series
#' @noRd
sc_plot_series <- function(object) {
  times <- object$spec$times
  time_names <- as.character(times)
  data.frame(
    time = times,
    treated = as.numeric(object$y_treated[time_names]),
    synthetic = as.numeric(object$y_synth[time_names]),
    effect = as.numeric(object$effect[time_names]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' @noRd
sc_plot_theme <- function() {
  ggplot2::theme_minimal()
}

#' @noRd
sc_plot_trline <- function(object) {
  ggplot2::geom_vline(
    xintercept = object$spec$trperiod,
    linetype = "dashed",
    color = SC_PLOT_CONTROL
  )
}

#' @noRd
sc_plot_path <- function(object) {
  series <- sc_plot_series(object)
  df <- data.frame(
    time = rep(series$time, 2L),
    value = c(series$treated, series$synthetic),
    series = factor(
      rep(c("Treated", "Synthetic"), each = nrow(series)),
      levels = c("Treated", "Synthetic")
    ),
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(df, ggplot2::aes(time, value, color = series)) +
    ggplot2::geom_line() +
    sc_plot_trline(object) +
    ggplot2::scale_color_manual(
      values = c(Treated = SC_PLOT_TREATED, Synthetic = SC_PLOT_CONTROL)
    ) +
    sc_plot_theme() +
    ggplot2::labs(x = object$spec$time, y = object$spec$outcome, color = NULL)
}

#' @noRd
sc_plot_effect <- function(object) {
  df <- sc_plot_series(object)
  ggplot2::ggplot(df, ggplot2::aes(time, effect)) +
    ggplot2::geom_hline(yintercept = 0, color = SC_PLOT_CONTROL) +
    ggplot2::geom_line(color = SC_PLOT_TREATED) +
    sc_plot_trline(object) +
    sc_plot_theme() +
    ggplot2::labs(x = object$spec$time, y = "Gap")
}

#' @noRd
sc_plot_v <- function(object) {
  v <- object$V
  df <- data.frame(
    predictor = factor(names(v), levels = names(v)),
    weight = as.numeric(v),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(df, ggplot2::aes(predictor, weight)) +
    ggplot2::geom_col(fill = SC_PLOT_TREATED) +
    sc_plot_theme() +
    ggplot2::labs(x = NULL, y = "V")
}

#' @noRd
sc_plot_w <- function(object) {
  df <- sc_result_nonzero_w(object$W)
  df$unit <- factor(df$unit, levels = df$unit)
  ggplot2::ggplot(df, ggplot2::aes(unit, weight)) +
    ggplot2::geom_col(fill = SC_PLOT_TREATED) +
    sc_plot_theme() +
    ggplot2::labs(x = NULL, y = "W")
}

#' @noRd
sc_plot_gaps <- function(object) {
  plc <- sc_plot_require(object, "placebo_space")
  effect <- plc$effect
  times <- as.numeric(rownames(effect))
  units <- colnames(effect)
  treat <- as.character(object$spec$treat)
  df <- data.frame(
    time = rep(times, ncol(effect)),
    unit = rep(units, each = nrow(effect)),
    gap = as.vector(effect),
    stringsAsFactors = FALSE
  )
  df$series <- ifelse(df$unit == treat, "Treated", "Placebo")
  placebos <- df[df$series == "Placebo", , drop = FALSE]
  treated <- df[df$series == "Treated", , drop = FALSE]
  ggplot2::ggplot(df, ggplot2::aes(time, gap, group = unit)) +
    ggplot2::geom_line(
      data = placebos,
      color = SC_PLOT_CONTROL,
      alpha = 0.55
    ) +
    ggplot2::geom_line(data = treated, color = SC_PLOT_TREATED) +
    sc_plot_trline(object) +
    sc_plot_theme() +
    ggplot2::labs(x = object$spec$time, y = "Gap")
}

#' @noRd
sc_plot_mspe <- function(object) {
  plc <- sc_plot_require(object, "placebo_space")
  df <- plc$table
  treat <- object$spec$treat
  df$series <- ifelse(df$unit == treat, "Treated", "Placebo")
  ord <- order(df$ratio, decreasing = TRUE)
  df <- df[ord, , drop = FALSE]
  df$unit <- factor(df$unit, levels = df$unit)
  ggplot2::ggplot(df, ggplot2::aes(unit, ratio, fill = series)) +
    ggplot2::geom_col() +
    ggplot2::scale_fill_manual(
      values = c(Treated = SC_PLOT_TREATED, Placebo = SC_PLOT_CONTROL)
    ) +
    sc_plot_theme() +
    ggplot2::labs(x = NULL, y = "MSPE ratio", fill = NULL)
}

#' @noRd
sc_plot_pvalue <- function(object) {
  plc <- sc_plot_require(object, "placebo_space")
  p_left <- plc$p_left
  df <- data.frame(
    time = as.numeric(names(p_left)),
    p_left = as.numeric(p_left),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  ggplot2::ggplot(df, ggplot2::aes(time, p_left)) +
    ggplot2::geom_line(color = SC_PLOT_TREATED) +
    sc_plot_trline(object) +
    sc_plot_theme() +
    ggplot2::labs(x = object$spec$time, y = "Left p-value")
}

#' @noRd
sc_plot_loo <- function(object) {
  band <- sc_plot_require(object, "loo")$band
  ggplot2::ggplot(band, ggplot2::aes(time, effect)) +
    ggplot2::geom_ribbon(
      ggplot2::aes(ymin = loo_min, ymax = loo_max),
      fill = SC_PLOT_CONTROL,
      alpha = 0.35
    ) +
    ggplot2::geom_hline(yintercept = 0, color = SC_PLOT_CONTROL) +
    ggplot2::geom_line(color = SC_PLOT_TREATED) +
    sc_plot_trline(object) +
    sc_plot_theme() +
    ggplot2::labs(x = object$spec$time, y = "Gap")
}

#' Plot an `scm_fit`
#'
#' ggplot2 charts for treated versus synthetic paths, the gap series,
#' predictor and donor weights, in-space placebos, and leave-one-out
#' bands. Does not refit.
#'
#' `type` is one of `"path"`, `"effect"`, `"v"`, `"w"`, `"gaps"`,
#' `"mspe"`, `"pvalue"`, or `"loo"`. Types `"gaps"`, `"mspe"`, and
#' `"pvalue"` need `$placebo_space`. Type `"loo"` needs `$loo`.
#'
#' @param object,x An `scm_fit`.
#' @param type Plot kind. Default `"path"`.
#' @param ... Passed from `plot()` to [autoplot()].
#'
#' @return A ggplot. `plot()` prints it and returns it invisibly.
#'
#' @examples
#' fit <- scm(
#'   cigsale ~ lnincome + cigsale(1988),
#'   data = smoking, unit = "state", time = "year",
#'   treat = 3, trperiod = 1989, xperiod = 1980:1988
#' )
#' ggplot2::autoplot(fit, type = "path")
#'
#' @export
#' @importFrom ggplot2 autoplot
autoplot.scm_fit <- function(object,
                             type = c(
                               "path", "effect", "v", "w",
                               "gaps", "mspe", "pvalue", "loo"
                             ),
                             ...) {
  type <- match.arg(type)
  switch(
    type,
    path = sc_plot_path(object),
    effect = sc_plot_effect(object),
    v = sc_plot_v(object),
    w = sc_plot_w(object),
    gaps = sc_plot_gaps(object),
    mspe = sc_plot_mspe(object),
    pvalue = sc_plot_pvalue(object),
    loo = sc_plot_loo(object)
  )
}

#' @export
#' @rdname autoplot.scm_fit
plot.scm_fit <- function(x, ...) {
  print(autoplot(x, ...))
}

utils::globalVariables(c(
  "time", "value", "series", "effect", "predictor", "weight",
  "unit", "gap", "ratio", "p_left", "loo_min", "loo_max"
))
