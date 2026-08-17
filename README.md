# synthaio

Classic Abadie–Diamond–Hainmueller synthetic control, plus a synth2-class
report: nested / allopt \(V\), donor weights, predictor balance, optional
in-space and in-time placebos, leave-one-out, and ggplot2 plots. The
numerical target is Stata `synth2` on the smoking / Proposition 99
specification.

## What this is not

- Not [gsynth](https://yiqingxu.org/packages/gsynth/) (interactive fixed
  effects, multiple treated units, staggered adoption).
- Not synthetic difference-in-differences (`synthdid`).
- Not `augsynth`, `scpi`, or a thin wrapper around `Synth` / `tidysynth`.

If you need those estimators, use those packages.

## Installation

From a local checkout:

```r
# install.packages("remotes")
remotes::install_local("path/to/synthaio")
```

From GitHub:

```r
# install.packages("remotes")
remotes::install_github("sunfang3/synthAIO")
```

## Example

```r
library(synthaio)

fit <- scm(
  cigsale ~ lnincome + age15to24 + retprice + beer +
    cigsale(1988) + cigsale(1980) + cigsale(1975),
  data = smoking, unit = "state", time = "year",
  treat = "California", trperiod = 1989, xperiod = 1980:1988,
  allopt = TRUE
)
print(fit)
balance(fit)
coef(fit)[coef(fit) > 0]
```

`treat` may be `3` or `"California"`. `allopt = TRUE` is multi-start nested optimization.
Donor set, ATT, RMSE, and the two largest \(V\) names match Stata
`synth2` 2.1.0 `nested allopt` on this panel
(`data-raw/synth2_smoking_allopt.log`). Passing the Stata \(V\) as
`custom_v` recovers the Stata \(W\) to `tol$w_abs`. See the vignettes
and `tests/testthat/test-gold-standard.R`.

```r
browseVignettes("synthaio")
```
