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

From GitHub (placeholder — replace `OWNER` when the public repo exists):

```r
# remotes::install_github("OWNER/synthaio")
```

## Example

```r
library(synthaio)

fit <- scm(
  cigsale ~ lnincome + age15to24 + retprice + beer +
    cigsale(1988) + cigsale(1980) + cigsale(1975),
  data = smoking, unit = "state", time = "year",
  treat = 3, trperiod = 1989, xperiod = 1980:1988,
  allopt = TRUE
)
print(fit)
balance(fit)
coef(fit)[coef(fit) > 0]
```

California is unit `3`. `allopt = TRUE` is multi-start nested optimization.
The donor *set* matches published `synth2, nested allopt` replicas; \(V\)
may land in a different local minimum. See the vignettes for the full walk
through, and `tests/testthat/test-gold-standard.R` for locked tolerances.

```r
browseVignettes("synthaio")
```
