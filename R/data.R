#' Cigarette consumption panel (Proposition 99)
#'
#' Annual state-level cigarette sales and covariates for 39 U.S. states,
#' 1970--2000, as used by Abadie, Diamond, and Hainmueller (2010) to study
#' California's Tobacco Control Program (Proposition 99). This is the
#' standard 39-state panel, not the 8-unit `Synth::synth.data` toy.
#'
#' `state` is the integer id from the source Stata file. California is 3.
#' The five `synth2, nested allopt` nonzero-weight donors are Utah (34),
#' Nevada (21), Montana (19), Colorado (4), and Connecticut (5).
#'
#' @format A data frame with 1209 rows and 7 variables:
#' \describe{
#'   \item{state}{Integer state id, 1--39. California is 3.}
#'   \item{year}{Calendar year, 1970--2000.}
#'   \item{cigsale}{Cigarette sales per capita, in packs.}
#'   \item{lnincome}{Log state per-capita income.}
#'   \item{age15to24}{Share of the state population aged 15--24.}
#'   \item{retprice}{Average retail price of cigarettes.}
#'   \item{beer}{Beer consumption per capita.}
#' }
#'
#' @source Abadie, Diamond, and Hainmueller (2010). Vendored from QuarCS
#'   `smoking_sc.dta`
#'   <https://github.com/quarcs-lab/data-open/raw/master/isds/smoking_sc.dta>.
#'   Rebuild notes and the public-replication license comment are in
#'   `data-raw/smoking.R`.
#'
#' @references
#' Abadie, A., Diamond, A., and Hainmueller, J. (2010).
#' Synthetic Control Methods for Comparative Case Studies: Estimating the
#' Effect of California's Tobacco Control Program.
#' *Journal of the American Statistical Association*, 105(490), 493--505.
#' \doi{10.1198/jasa.2009.ap08746}
#'
#' @examples
#' data("smoking", package = "synthaio")
#' unique(smoking$state)
#' subset(smoking, state == 3L & year %in% c(1975L, 1980L, 1988L))
"smoking"
