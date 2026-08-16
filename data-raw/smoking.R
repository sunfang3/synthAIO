# Rebuild the packaged smoking panel from the committed CSV.
#
# Source: Abadie, Diamond, and Hainmueller (2010), 39-state Proposition 99
# panel. Vendored from QuarCS `smoking_sc.dta` (public dataset collection):
#   https://github.com/quarcs-lab/data-open/raw/master/isds/smoking_sc.dta
# Do not use Synth::synth.data (8-unit toy, not Prop 99).
#
# License: public academic replication extract of published statistical
# series used in ADH (2010). The QuarCS copy is redistributed as a public
# dataset. Underlying series (cigarette sales, prices, income, beer,
# age shares) are compiled public statistics as cited in that paper.
# This file is committed so CI never hits the network.
#
# One-time vendor (not run by this script):
#   curl -fsSL -o /tmp/smoking_sc.dta \
#     https://github.com/quarcs-lab/data-open/raw/master/isds/smoking_sc.dta
# Then read with haven, drop labels, reorder columns to
# state, year, cigsale, lnincome, age15to24, retprice, beer, and write
# data-raw/smoking.csv.
#
# Stata value labels on `state` (California is 3):
#   1 Alabama, 2 Arkansas, 3 California, 4 Colorado, 5 Connecticut,
#   6 Delaware, 7 Georgia, 8 Idaho, 9 Illinois, 10 Indiana, 11 Iowa,
#   12 Kansas, 13 Kentucky, 14 Louisiana, 15 Maine, 16 Minnesota,
#   17 Mississippi, 18 Missouri, 19 Montana, 20 Nebraska, 21 Nevada,
#   22 New Hampshire, 23 New Mexico, 24 North Carolina, 25 North Dakota,
#   26 Ohio, 27 Oklahoma, 28 Pennsylvania, 29 Rhode Island,
#   30 South Carolina, 31 South Dakota, 32 Tennessee, 33 Texas, 34 Utah,
#   35 Vermont, 36 Virginia, 37 West Virginia, 38 Wisconsin, 39 Wyoming.

smoking <- read.csv(
  "data-raw/smoking.csv",
  colClasses = c(
    state = "integer",
    year = "integer",
    cigsale = "numeric",
    lnincome = "numeric",
    age15to24 = "numeric",
    retprice = "numeric",
    beer = "numeric"
  ),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(smoking) == 1209L,
  length(unique(smoking$state)) == 39L,
  identical(range(smoking$year), c(1970L, 2000L)),
  3L %in% smoking$state
)

print(sort(unique(smoking$state)))

dir.create("data", showWarnings = FALSE)
save(smoking, file = "data/smoking.rda", version = 3L, compress = "bzip2")
