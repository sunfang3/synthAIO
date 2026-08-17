# synthaio 0.0.0.9000

## Gold standard

* Smoking / Prop 99 `nested allopt` numbers are locked to the author
  Stata log `data-raw/synth2_smoking_allopt.log` (synth2 2.1.0).
* Feeding the Stata \(V\) as `custom_v` recovers the Stata \(W\).
* Default `allopt` lands in the same \(V\) basin (top-2 names).

## User-facing

* `treat = "California"` works when `state_name` is on the panel.
* `print()` / `summary()` show `California (3)` and `Utah (34)`.
* `summary()` prints in-time placebo gaps and the leave-one-out band.
* `placebo = list(unit = c(4, 5))` runs only those donors.
* In-time placebos clip explicit `mspeperiod` / `preperiod`.
* `sc_effects()` is the preferred extractor; `effects()` remains an alias.

## Infrastructure

* GitHub Actions runs `R CMD check` on push and pull request.
