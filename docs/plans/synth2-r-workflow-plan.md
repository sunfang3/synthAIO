# synth2 R Workflow Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship an R package `synthaio` whose `scm()` entry point estimates classic Abadie SCM and emits a synth2-class report (weights, balance, placebos, LOO, plots) that matches Stata `synth2` on the smoking / Prop 99 specification.

**Architecture:** Long-panel data is parsed into an `sc_spec` (predictor matrix \(X_0,X_1\), outcome paths \(Y_0,Y_1\)). A solver layer computes \(W\mid V\) by constrained QP and \(V\) by regression / nested / allopt. `sc_fit()` turns one spec into an `scm_fit` object. Placebo and LOO re-call `sc_fit()` on mutated specs. `scm()` orchestrates and attaches optional diagnostics. No dependency on `Synth`, `tidysynth`, or `SCtools` for estimation.

**Tech Stack:** R ≥ 4.2, `osqp` (inner QP), `stats::optim` (outer \(V\)), `ggplot2` (plots), `testthat` + fixture JSON (gold standard), GPL-3. Working package name `synthaio`; exported function `scm()` — never `synth2` (CRAN `seewave::synth2`).

**Spec:** `docs/brainstorms/synth2-r-workflow.md` (R1–R18).

**Solver default (Q3):** implement the QP ourselves; do not wrap `tidysynth`. If smoking gold standard is missed after Unit 5, refine scaling / starts / `sigf` in Unit 5 before inventing a new estimator.

---

## Architecture (one page)

```
long data + formula + trunit/trperiod
            │
            ▼
        sc_spec()          R/sc_spec.R
            │  X0, X1, Y0, Y1, scaled copies, windows
            ▼
     sc_solve_w(V)         R/sc_solve_w.R     QP: min (X1-X0w)'V(X1-X0w)
            │                                    s.t. w≥0, 1'w=1
            ▼
     sc_solve_v()          R/sc_solve_v.R     regression | nested | allopt | custom
            │
            ▼
        sc_fit()           R/sc_fit.R         W, V, paths, ATT, RMSE, R², balance
            │
     ┌──────┼──────────┐
     ▼      ▼          ▼
 placebo  loo        scm()                 R/sc_placebo.R  R/sc_loo.R  R/scm.R
     └──────┴──────────┘
            ▼
        scm_fit            R/sc_result.R + R/sc_plot.R
```

**Scaling (must match Synth / Stata `synth`, not optional):** scale each predictor by the **sample** standard deviation (`n-1`) of that row in `cbind(X0, X1)` (treated column included). If gold \(W\) still misses in Unit 12, the first fallback is donor-only population sd. \(V\) is normalized to \(\sum_k V_{kk}=1\). Inner QP uses the scaled \(X\). Reported balance tables use the unscaled predictor means.

**Regression \(V\) (default when `nested = FALSE`):** do **not** OLS \(X_1\) on \(X_0\), and do **not** collapse \(Z\) to a scalar per unit. Port `Synth::synth` exactly (cite `synth.R` in a comment):

```r
Xall <- cbind(rep(1, 1 + J), t(cbind(X1.scaled, X0.scaled)))  # (1+J) × (1+K)
Zall <- cbind(Z1, Z0)                                        # T_mspe × (1+J)
B    <- solve(crossprod(Xall), crossprod(Xall, t(Zall)))     # (1+K) × T_mspe
B    <- B[-1, , drop = FALSE]                                # drop intercept
v0   <- diag(tcrossprod(B))                                  # V_kk ∝ sum_t β_{k,t}^2
V    <- v0 / sum(v0)
```

Do **not** floor tiny \(V_{kk}\) at ε — the gold fixture has `lnincome: 0.0`. If `sum(v0) == 0`, fall back to equal \(V = 1/K\). Nested then replaces this \(V\) by minimizing outcome MSPE over `mspeperiod`.

---

## File map

| Path | Responsibility |
|------|----------------|
| `DESCRIPTION`, `NAMESPACE`, `LICENSE` | Package metadata, GPL-3, `Roxygen: list(markdown = TRUE)` |
| `R/synthaio-package.R` | Package doc, `_PACKAGE` |
| `R/sc_spec.R` | Formula + panel → `sc_spec` |
| `R/sc_solve_w.R` | \(W\mid V\) QP + `margin`/`maxiter`/`sigf`/`bound` |
| `R/sc_solve_v.R` | regression / nested / allopt / custom \(V\) |
| `R/sc_fit.R` | One-spec estimate → core `scm_fit` fields |
| `R/sc_placebo.R` | In-space, in-time; reserved mixed args |
| `R/sc_loo.R` | Leave-one-out on nonzero-\(W\) donors |
| `R/scm.R` | User entry `scm()` |
| `R/sc_result.R` | `print` / `summary` / extractors |
| `R/sc_plot.R` | ggplot2 methods |
| `R/data.R` | `smoking` dataset documentation |
| `data/smoking.rda` | Vendored Prop 99 panel |
| `data-raw/smoking.csv`, `data-raw/smoking.R` | Rebuild script |
| `tests/testthat/fixtures/smoking-synth2-nested-allopt.json` | Gold numbers |
| `vignettes/prop99.Rmd`, `vignettes/prop99-zh.Rmd` | Teaching replica |

Do not create `R/utils.R` dumping grounds. Tiny shared helpers live next to their only caller; if a helper is used twice, put it in `R/sc_spec.R` (data) or `R/sc_fit.R` (metrics).

---

## Unit DAG

```
U1 skeleton → U2 data/fixtures → U3 spec → U4 W|V → U5 V methods
                                                      ↓
                                                    U6 fit
                                                      ↓
                                    U7 placebos    U8 LOO
                                      ↓               ↓
                                      └──── U9 scm() ─┘
                                              ↓
                                    U10 result methods
                                              ↓
                                    U11 plots
                                              ↓
                                    U12 gold-standard + check
                                              ↓
                                    U13 vignettes
```

---

### Unit 1: Package skeleton

**Goal:** Installable R package named `synthaio` that `R CMD check` can see.

**Requirements trace:** R17

**Dependencies:** none

**Files:**
- `DESCRIPTION` — package metadata
- `LICENSE` / `LICENSE.md` — GPL-3
- `NAMESPACE` — roxygen-managed; empty exports at first
- `R/synthaio-package.R` — package doc
- `.Rbuildignore` — ignore `docs/`, `data-raw/`
- `tests/testthat.R`, `tests/testthat/test-skeleton.R`

**Approach:** `usethis`-style layout by hand (do not require usethis at runtime). `Package: synthaio`. `Title: Classic Synthetic Control Workflow`. `License: GPL-3`. `Depends: R (>= 4.2.0)`. `Imports: osqp, ggplot2`. `Suggests: testthat (>= 3.0.0), knitr, rmarkdown`. `Config/testthat/edition: 3`. `RoxygenNote` set after first `devtools::document()`. Encoding UTF-8. LazyData true once Unit 2 adds data.

**Patterns:** Standard CRAN package layout. No README claims until Unit 13.

**Test scenarios:**
- [ ] Happy path: `devtools::load_all()` succeeds; `packageVersion("synthaio")` is `0.0.0.9000`
- [ ] DESCRIPTION has no `synth2` package or function name
- [ ] `R CMD check` on the empty package has no ERROR (warnings about missing docs are fixed before leaving this unit)

**Verification:** `Rscript -e 'devtools::check(document = FALSE, args = "--no-manual")'` exits 0 or only has notes we accept (e.g. new submission).

**Planning-time unknowns:** none.

- [ ] **Step 1:** Create the files above. `DESCRIPTION` must contain:

```
Package: synthaio
Type: Package
Title: Classic Synthetic Control Workflow
Version: 0.0.0.9000
Authors@R: person("synthAIO", "Authors", role = c("aut", "cre"),
                  email = "devnull@example.com")
Description: Estimates the Abadie-Diamond-Hainmueller synthetic control
    and produces a synth2-class report: nested V-W optimization, in-space
    and in-time placebos, leave-one-out robustness, balance tables, and
    plots. Numerical targets are Stata synth2 on standard examples.
License: GPL-3
Encoding: UTF-8
Language: en-US
Roxygen: list(markdown = TRUE)
Depends: R (>= 4.2.0)
Imports: osqp, ggplot2
Suggests: testthat (>= 3.0.0), knitr, rmarkdown
Config/testthat/edition: 3
```

Use a real maintainer email if the user supplies one; otherwise keep the placeholder and add a comment in DESCRIPTION that it must change before CRAN.

- [ ] **Step 2:** Add `tests/testthat/test-skeleton.R`:

```r
test_that("package name is synthaio not synth2", {
  expect_equal(unname(desc::desc_get("Package")), "synthaio")
  ns <- parseNamespaceFile(".", package = "synthaio")
  expect_false("synth2" %in% ns$exports)
})
```

If `desc` is not a dependency, read `DESCRIPTION` with `read.dcf` instead.

- [ ] **Step 3:** `devtools::test()` for that file passes. Commit:

```
feat: scaffold synthaio package skeleton
```

---

### Unit 2: Smoking data and gold-standard fixtures

**Goal:** Package contains the 39-state smoking panel and a locked JSON of published `synth2, nested allopt` numbers.

**Requirements trace:** R8, R16

**Dependencies:** Unit 1

**Files:**
- `data-raw/smoking.csv` — 1209 rows
- `data-raw/smoking.R` — writes `data/smoking.rda`
- `data/smoking.rda`
- `R/data.R` — roxygen for `smoking`
- `tests/testthat/fixtures/smoking-synth2-nested-allopt.json`
- `tests/testthat/helper-fixtures.R` — `gold_smoking()` loader
- `tests/testthat/test-smoking-data.R`

**Approach:** Vendor the standard Abadie et al. (2010) **39-state** panel (`state`, `year`, `cigsale`, `lnincome`, `age15to24`, `retprice`, `beer`). Source: QuarCS `smoking_sc.dta` or RePEc `smoking.dta`, or `tidysynth::smoking` if already installed. **Do not** vendor `Synth::synth.data` — that is an 8-unit toy, not Prop 99. Commit CSV so CI never hits the network. Document source and license in `data-raw/smoking.R` comments.

Fixture JSON is copied from Carlos Mendez's `synth2, nested allopt` replica (synth2 2.1.0) and the published MSPE table. Header comment in the JSON (`"_source"` field) records origin. Tolerances live in the JSON, not in test code.

```json
{
  "_source": "Carlos Mendez tutorial replica of Yan-Chen synth2 2.1.0 nested allopt; replace with author Stata log when available",
  "_spec": "cigsale ~ lnincome + age15to24 + retprice + beer + cigsale(1988)+cigsale(1980)+cigsale(1975); trunit=California; trperiod=1989; xperiod=1980:1988; nested allopt",
  "rmse": 1.75567,
  "r2": 0.97434,
  "att": -19.0018,
  "v": {
    "lnincome": 0.0,
    "age15to24": 0.5459,
    "retprice": 0.0174,
    "beer": 0.0031,
    "cigsale_1988": 0.0049,
    "cigsale_1980": 0.0066,
    "cigsale_1975": 0.4221
  },
  "w_nonzero": {
    "Utah": 0.3340,
    "Nevada": 0.2350,
    "Montana": 0.2020,
    "Colorado": 0.1610,
    "Connecticut": 0.0680
  },
  "effect": {
    "1989": -7.5945,
    "1990": -9.7039,
    "1993": -17.7897,
    "1997": -23.9123,
    "1999": -26.3711,
    "2000": -25.7550
  },
  "mspe_placebo_nested_sigf6": {
    "_note": "From the *placebo* run (nested, no allopt, sigf=6), NOT the allopt baseline. When windows coincide, allopt pre_mspe should be ≈ rmse^2 ≈ 3.082, not 3.1668. Unit 12 must not retune the solver to match this block.",
    "pre": 3.1668,
    "post": 391.2533,
    "ratio": 123.5490
  },
  "p_unfiltered": 0.025641025641,
  "tol": {
    "w_abs": 0.001,
    "v_abs": 0.001,
    "rmse_rel": 0.01,
    "att_abs": 0.05
  }
}
```

**Patterns:** Fixtures are data, not code. Tests load via helper.

**Test scenarios:**
- [ ] Happy path: `smoking` has 39 states, years 1970–2000, 1209 rows, California present
- [ ] Nil/empty: `gold_smoking()` errors with a clear message if the JSON is missing
- [ ] Schema: fixture has `w_nonzero`, `v`, `tol`, `_source`

**Verification:** `devtools::test(filter = "smoking-data")` passes.

**Planning-time unknowns:** Exact numeric state labels (California is 3 in the Stata file). **Resolve before Unit 3:** after vendoring, print `unique(smoking$state)` / labels and lock California's id in the fixture as both name and code.

- [ ] **Step 1:** Write `test-smoking-data.R` first (panel shape assertions) — it fails because `smoking` does not exist.
- [ ] **Step 2:** Vendor data + JSON + helper until those tests pass.
- [ ] **Step 3:** Commit `test: add smoking panel and synth2 gold fixture`

---

### Unit 3: `sc_spec()` — panel and predictor windows

**Goal:** Turn a long panel plus a synth2-like formula into matrices an estimator can consume.

**Requirements trace:** R3, R4

**Dependencies:** Unit 2

**Files:**
- `R/sc_spec.R`
- `tests/testthat/test-sc-spec.R`

**Approach:** One S3 class `sc_spec` (list). Constructor:

```r
sc_spec(
  formula,
  data,
  unit,
  time,
  treat,
  trperiod,
  counit = NULL,
  xperiod = NULL,
  mspeperiod = NULL,
  preperiod = NULL,
  postperiod = NULL
)
```

`unit`, `time`, `treat` are unquoted names or strings. `treat` is the treated unit id (numeric or character matching `data[[unit]]`).

**Formula grammar (lock this, do not invent a second syntax):**

- LHS is the outcome.
- Bare `lnincome` → mean of that variable over `xperiod` (default: all pre-treatment times).
- `cigsale(1988)` → outcome (or any panel var) at time 1988. Implemented as a call in the formula; the parser walks `terms`/`call` and treats a single numeric argument as a time point.
- `cigsale(1980:1988)` → mean over that integer sequence.
- `cigsale(c(1982, 1986, 1988))` → mean over those times (Stata `var(t1&t2&t3)`).
- Do not support Stata `t1&t2` inside R parse.

Predictor names in the spec become `lnincome`, `cigsale_1988`, `cigsale_1980_1988`, `cigsale_1982_1986_1988` as needed. These names must match fixture keys after a documented mapping.

**Matrices:**

- `X1` length-\(K\) treated predictor vector; `X0` \(K \times J\) donors.
- `Y1` full time path of treated outcome; `Y0` \(T \times J\).
- `X1_scaled`, `X0_scaled`: each predictor row divided by the **sample** sd (`n-1`) of that row in `cbind(X0, X1)` (treated included). This is the default, not the fallback. Donor-only population sd is only a Unit 12 rescue if gold \(W\) still misses.
- Slots: `donors`, `times`, `pre_times`, `post_times`, `mspe_times`, `x_times`, `outcome`, `treat`, `trperiod`.
- Persist constructor inputs on the object (`formula`, `data` or a handle, `unit`, `time`, `treat`, `xperiod`, `mspeperiod`, `counit`, `preperiod`, `postperiod`) so Units 7–8 rebuild specs by mutating those inputs, not by reverse-engineering matrices.

**Window defaults (lock; do not treat these as aliases of each other):**

| Argument | Default when `NULL` | Used for |
|----------|---------------------|----------|
| `preperiod` | all times `< trperiod` | RMSE, \(R^2\), path highlighting |
| `postperiod` | all times `>= trperiod` | ATT, post MSPE |
| `xperiod` | all times `< trperiod` | averaging bare covariates |
| `mspeperiod` | all times `< trperiod` (**not** `xperiod`) | regression \(Z\), nested/allopt objective, pre MSPE |

Gold spec: `xperiod = 1980:1988`, `mspeperiod` omitted → optimize 1970–1988. If an implementer defaults `mspeperiod` to `xperiod` once `xperiod` is set, Units 5–12 will miss gold for a spec bug, not a solver bug.

Predictor-specific times (`cigsale(1988)`) always override `xperiod` for that predictor only.

Errors (named, user-visible):

- `synthaio_not_panel` — duplicate `unit`×`time`
- `synthaio_missing_treat` — treat id not in data
- `synthaio_bad_time` — predictor time not in panel
- `synthaio_single_donor` — fewer than 2 donors
- `synthaio_unbalanced` — a required unit-time is missing (v1 requires balanced on the used times)

**Patterns:** Pure functions; no ggplot, no solver.

**Test scenarios:**
- [ ] Happy path: smoking Prop 99 formula yields \(K=7\) predictors, \(J=38\) donors, `X1` cigsale_1988 ≈ 90.1
- [ ] Nil/empty: `data` with 0 rows → `synthaio_not_panel` (or explicit empty-data error)
- [ ] Error path: `cigsale(1960)` when 1960 absent → `synthaio_bad_time`
- [ ] Edge: `counit` restricts donors; treated unit listed in `counit` is dropped with a warning
- [ ] Edge: default `xperiod` is all times `< trperiod`
- [ ] Edge: `xperiod = 1980:1988` and `mspeperiod = NULL` → `x_times` is 1980:1988 and `mspe_times` is all times `< trperiod` (1970:1988 on smoking)

**Verification:** `devtools::test(filter = "sc-spec")` passes. An engineer can `str(spec)` and see `X0`, `X1`, `Y0`, `Y1`.

**Planning-time unknowns:** How beer / age / income missingness is averaged (`na.rm = TRUE`, matching `synth`). **Deferred to implementation** — use `na.rm = TRUE` and document.

- [ ] **Step 1:** Write `test-sc-spec.R` covering the five scenarios. Fail on missing `sc_spec`.
- [ ] **Step 2:** Implement parser + matrices + class. No solver.
- [ ] **Step 3:** Commit `feat: parse long-panel SCM specifications`

---

### Unit 4: `sc_solve_w()` — inner QP

**Goal:** Given scaled \(X_0,X_1\) and diagonal \(V\), return \(W\) on the simplex.

**Requirements trace:** R1, R14

**Dependencies:** Unit 3 (only for integration tests; the solver itself is matrix-in / vector-out)

**Files:**
- `R/sc_solve_w.R`
- `tests/testthat/test-sc-solve-w.R`

**Approach:**

\[
\min_w\ (X_1-X_0 w)^\top V (X_1-X_0 w)
\quad w\ge 0,\ \mathbf{1}^\top w=1
\]

QP: \(P = X_0^\top V X_0\) (symmetrized), \(q = -X_0^\top V X_1\). Solve with `osqp` (`verbose = FALSE`, `eps_abs`/`eps_rel` derived from `sigf`: `10^(-sigf)`). `bound` clips \(w_j\) to \([0, \texttt{bound}]\) then the simplex equality still holds — for classic SCM `bound` defaults to 10 and is inert because \(w_j\le 1\). `margin` is the extra primal-feasibility slack passed to osqp as `eps_prim_inf` / a relaxed equality tolerance of `margin * 10^(-sigf)` — do not leave it unused and do not invent a third meaning in Unit 12. After solve, zero out weights below `10^(-sigf)` and renormalize if the remainder sums to > 0.

Function signature:

```r
sc_solve_w(X0, X1, V, margin = 0.05, maxiter = 1000, sigf = 7, bound = 10)
```

`V` may be a length-\(K\) vector or \(K\times K\) diagonal matrix. Return a named numeric vector of length \(J\) (names from `colnames(X0)`).

Do not call `Synth::synth`.

**Patterns:** No data frames inside the solver.

**Test scenarios:**
- [ ] Happy path: two donors, \(X_1\) exactly equal to donor 1 → \(w \approx (1,0)\)
- [ ] Happy path: \(X_1\) midpoint of two identical-dimension donors → \(w \approx (0.5, 0.5)\)
- [ ] Nil/empty: `X0` with 0 columns → error `synthaio_single_donor`
- [ ] Error path: `V` length ≠ `nrow(X0)` → error
- [ ] Edge: `sum(w) == 1`, `all(w >= -10^(-sigf))`
- [ ] Edge: `custom V` that zeros a predictor does not use that row (weight on that discrepancy is 0)

**Verification:** unit tests pass; a smoke call on `sc_spec(smoking, ...)` with equal \(V\) returns a simplex vector of length 38.

**Planning-time unknowns:** osqp vs Stata interior-point residual. **Deferred to Unit 5/12** — if gold \(W\) set is wrong, try tightening `sigf`, adding a simplex projection polish, or a fallback `quadprog::solve.QP` and keep the lower objective.

- [ ] **Step 1:** Write the six tests with a tiny handmade \(X_0,X_1\).
- [ ] **Step 2:** Implement `sc_solve_w()`.
- [ ] **Step 3:** Commit `feat: solve donor weights given V`

---

### Unit 5: `sc_solve_v()` — regression, nested, allopt, custom

**Goal:** Produce diagonal \(V\) by each synth2 method and the \(W\) that goes with it.

**Requirements trace:** R2, R8 (partial), R14

**Dependencies:** Unit 4

**Files:**
- `R/sc_solve_v.R`
- `tests/testthat/test-sc-solve-v.R`

**Approach:**

```r
sc_solve_v(spec, method = c("regression", "nested", "allopt", "custom"),
           custom_v = NULL, ...)
```

Returns `list(V = named numeric length K, W = named numeric length J, mspe = scalar, starts = ...)`.

1. **custom:** normalize `custom_v` to sum 1; `W = sc_solve_w(...)`.
2. **regression:** starting \(V\) from the locked `Xall`/`Zall`/`diag(BB')` block in Architecture (copy that code; do not re-derive). `W = sc_solve_w`. Default when `nested = FALSE`.
3. **nested:** parameterize \(\theta \in \mathbb{R}^{K-1}\) (or \(\mathbb{R}^K\) with softmax) so \(V_k = e^{\theta_k}/\sum e^{\theta_j}\). Objective = outcome MSPE over `mspeperiod` using `W(V)`. Optimize with `stats::optim(..., method = "BFGS")`, start from regression \(\theta\). Guard `maxiter`.
4. **allopt:** run nested from three starts — (a) regression \(V\), (b) equal \(V=1/K\), (c) `optim` with `fnscale` and a jittered start (`rnorm` on \(\theta\), seed documented as `start_seed = 1` for reproducibility). Keep the lowest MSPE.

Always return \(V\) summing to 1.

**Patterns:** Outer optimizer may only call `sc_solve_w` + a small `sc_mspe(spec, W)` helper in `R/sc_fit.R` if that helper already exists; until Unit 6, put `sc_pre_mspe(Y1, Y0, W, mspe_index)` in `R/sc_solve_v.R`.

**Test scenarios:**
- [ ] Happy path: `custom_v` equal weights on a 2-donor toy spec is a no-op normalize
- [ ] Happy path: `nested` MSPE ≤ `regression` MSPE on smoking (allow equality)
- [ ] Happy path: `allopt` MSPE ≤ `nested` MSPE
- [ ] Nil/empty: `custom` with `custom_v = NULL` → error
- [ ] Error path: `custom_v` with negatives → error
- [ ] Edge: `allopt` is deterministic given `start_seed`

**Verification:** On smoking, `allopt` produces at least 3 of the 5 gold nonzero donors (Utah, Nevada, Montana, Colorado, Connecticut). If fewer than 3, **stop and fix the solver** before Unit 6. Do not proceed with a wrong donor set.

**Planning-time unknowns:** None. Copy the Architecture `Xall`/`Zall`/`B`/`diag(BB')` block. Cite `Synth::synth` file/line. Do not OLS \(X_1\) on \(X_0\). Do not average \(Z\) first.

- [ ] **Step 1:** Toy tests for custom + simplex.
- [ ] **Step 2:** Implement custom + regression + nested.
- [ ] **Step 3:** Add smoking smoke test (donor-set overlap ≥ 3) and `allopt`.
- [ ] **Step 4:** Commit `feat: estimate V by regression, nested, and allopt`

---

### Unit 6: `sc_fit()` — paths, ATT, RMSE, R², balance

**Goal:** One estimated spec becomes a complete point-estimate object (no placebos yet).

**Requirements trace:** R1, R9

**Dependencies:** Unit 5

**Files:**
- `R/sc_fit.R`
- `tests/testthat/test-sc-fit.R`

**Approach:**

```r
sc_fit(spec, method = "regression", custom_v = NULL, ...)
```

Calls `sc_solve_v`, then:

- `y_synth = Y0 %*% W`
- `effect = Y1 - y_synth` (named by time)
- `att = mean(effect[post])`
- `rmse = sqrt(mean(effect[pre]^2))`
- `r2 = 1 - ss_res / ss_tot` on **pre-treatment** path (`ss_tot` = sum of squared deviations of \(Y_1^{pre}\) from its mean). Ignore the synth2 help typo that says “posttreatment”.
- `pre_mspe` = mean squared effect over `mspe_times` (nested objective; gold `mspe.pre`)
- `rmse` / `r2` = over `pre_times` (metric window; may equal `mspe_times` under defaults)
- `post_mspe` = mean squared effect over `post_times`; `mspe_ratio = post_mspe / pre_mspe`
- `balance`: data.frame with columns `predictor`, `v_weight`, `treated`, `synthetic`, `bias_pct`, `avg_control`, `avg_bias_pct`
- class `scm_fit`

`avg_control` is the unweighted mean of donors (not \(W\)-weighted). `bias_pct = 100 * (synthetic - treated) / treated` when treated ≠ 0.

This is the **only** new public-ish abstraction in this unit (`scm_fit`). `sc_fit` is an internal/exported low-level function; users normally call `scm()`.

**Patterns:** Metrics stay in `R/sc_fit.R`. No plotting.

**Test scenarios:**
- [ ] Happy path: smoking + `allopt` has `rmse > 0`, `att < 0`, 7-row balance, 5-ish nonzero weights
- [ ] Happy path: `effect["1988"]` near 0 relative to `effect["2000"]`
- [ ] Nil/empty: not applicable (spec already validated)
- [ ] Error path: passing a non-`sc_spec` → error
- [ ] Edge: `r2` in \((0,1]\) on smoking

**Verification:** `summary` is not required yet; `str(fit)` shows `W`, `V`, `att`, `balance`.

**Planning-time unknowns:** R² definition. **Resolved above** (pre-treatment). If gold `r2 = 0.97434` is far off, re-read SJ paper and adjust one definition — do not invent a second `r2` field.

- [ ] **Step 1:** Tests on smoking for `att < 0`, balance rows, class.
- [ ] **Step 2:** Implement `sc_fit`.
- [ ] **Step 3:** Commit `feat: build scm_fit point estimates and balance`

---

### Unit 7: Placebos (in-space, in-time; mixed reserved)

**Goal:** Automate fake-unit and fake-time re-estimation and Fisher / MSPE tables.

**Requirements trace:** R5, R6, R13 (interface only)

**Dependencies:** Unit 6

**Files:**
- `R/sc_placebo.R`
- `tests/testthat/test-sc-placebo.R`

**Approach:**

```r
sc_placebo_space(spec, fit, cut = Inf, ...)
sc_placebo_time(spec, period, ...)
```

**In-space:** for each donor \(j\), rebuild a spec with `treat = j` and the remaining units (including the original treated) as donors, `sc_fit` with the **same** `method`. Collect pre/post MSPE and ratio. **Locked Fisher p:** \(p = \#\{\text{ratio} \ge \text{ratio}_{tr}\} / N_{\text{units}}\) including the treated unit, matching synth2's `1/39 = 0.0256`. Do not use a mid-p or `+1` in the denominator. Document that this is the unfiltered p.

`cut = 2` drops placebos whose **pre** MSPE \(> 2 \times\) treated pre MSPE before ranking and before year-wise p. Year-wise p uses the **same** include-treated convention as Fisher p, so the floor is \(1/n_{\mathrm{kept}}\) (Mendez `cut(2)` → \(1/20 = 0.05\)). Left-sided p at time \(t\): share of retained units (treated + kept placebos) with `effect[t] <= effect_tr[t]`. Two-sided: share with `|effect[t]| >= |effect_tr[t]|`.

**In-time:** user supplies fake `period` \(t_0 < trperiod\). Rebuild spec with `trperiod = t_0`. Point predictors whose time is \(\ge t_0\) (e.g. `cigsale(1988)` when \(t_0=1985\)) are **dropped**. Bare covariates whose `xperiod` crosses \(t_0\) are **re-averaged on times in `xperiod` that are `< t_0`** (not dropped). If `xperiod` was `NULL`, it already defaults to times `< trperiod`, so after the shift it becomes times `< t_0`. If no predictors remain, error `synthaio_no_predictors`. Fit. Return fake-window effects and the original post window evaluated at the new weights.

**Mixed (R13):** `scm(..., placebo = list(unit = TRUE, period = 1985, cut = 10))` is parsed in Unit 9. If both `unit` and `period` are set, v1 runs **in-time first** (shift `trperiod`) then in-space on that shifted spec — this matches the synth2 mixed example that sets `trperiod(1985)` and `placebo(unit ...)`. Implement that composition; do not invent a third algorithm. Add a test that the composition is invoked, not a full smoking mixed run (too slow for default testthat unless we skip on CRAN with `skip_on_cran` for the long one).

Store results as a list on the object: `placebo_space`, `placebo_time` with tidy data frames (`unit`, `pre_mspe`, `post_mspe`, `ratio`, `kept`).

**Patterns:** Reuse `sc_spec` + `sc_fit`. No copy-paste of QP.

**Test scenarios:**
- [ ] Happy path: 3-unit toy panel, treated has the unique post jump → treated has the largest ratio
- [ ] Happy path: unfiltered \(p = 1/3\)
- [ ] Nil/empty: `cut` drops every placebo → year-wise p is `NA` with warning `synthaio_cut_empty`
- [ ] Error path: in-time `period >= trperiod` → error
- [ ] Edge: `cut = 2` on toy data where one donor has huge pre MSPE excludes it
- [ ] Edge: mixed args call in-time-then-space (spy via a tiny toy, not smoking)
- [ ] Edge: in-time \(t_0\) drops a point predictor `y(t>=t0)` and re-averages a bare covariate on `xperiod < t0`

**Verification:** Toy p-values match the locked formula. Smoking full placebo is **not** required in this unit (too slow / depends on solver quality); it belongs in Unit 12 with `skip_on_cran()`.

**Planning-time unknowns:** Exact Fisher formula (mid-p vs +1). **Resolved above** to match \(1/39\). If a later Stata log uses a different convention, change one function and the fixture together.

- [ ] **Step 1:** Build a 3×T toy inside the test file (do not add package data).
- [ ] **Step 2:** Implement space + time + mixed composition.
- [ ] **Step 3:** Commit `feat: automate in-space and in-time placebos`

---

### Unit 8: Leave-one-out

**Goal:** Drop each nonzero-weight donor once and record effect bands.

**Requirements trace:** R7

**Dependencies:** Unit 6

**Files:**
- `R/sc_loo.R`
- `tests/testthat/test-sc-loo.R`

**Approach:**

```r
sc_loo(spec, fit, ...)
```

Let \(S = \{j : W_j > 10^{-sigf}\}\). For each \(j \in S\), rebuild spec with `counit = donors \ {j}`, `sc_fit` same method. Return a data frame `time, effect, loo_min, loo_max` and a list of per-dropped-unit `W`.

If \(|S|=0\) (should not happen after a valid fit), error `synthaio_loo_empty`.

**Patterns:** Same spec mutation style as placebos.

**Test scenarios:**
- [ ] Happy path: 3-donor toy with two positive weights → two LOO fits; band contains the original effect
- [ ] Nil/empty: force `W` all below threshold → `synthaio_loo_empty`
- [ ] Error path: non-`scm_fit` → error
- [ ] Edge: dropping the largest donor still produces \(\sum W=1\)

**Verification:** Toy tests pass.

**Planning-time unknowns:** none.

- [ ] **Step 1:** Toy tests.
- [ ] **Step 2:** Implement `sc_loo`.
- [ ] **Step 3:** Commit `feat: leave-one-out donor robustness`

---

### Unit 9: `scm()` single entry

**Goal:** One function runs estimate and optional placebo / LOO.

**Requirements trace:** R11, R13 (flags), R14 (optimizer args)

**Dependencies:** Units 6, 7, 8

**Files:**
- `R/scm.R`
- `tests/testthat/test-scm.R`

**Approach:**

```r
scm(
  formula, data, unit, time, treat, trperiod,
  counit = NULL, xperiod = NULL, mspeperiod = NULL,
  preperiod = NULL, postperiod = NULL,
  nested = FALSE, allopt = FALSE, custom_v = NULL,
  placebo = NULL,
  loo = FALSE,
  margin = 0.05, maxiter = 1000, sigf = 7, bound = 10,
  start_seed = 1
)
```

`method` mapping: `allopt = TRUE` → `"allopt"`; else `nested = TRUE` → `"nested"`; else if `custom_v` → `"custom"`; else `"regression"`. `allopt` implies nested (ignore a contradictory `nested = FALSE` with a warning).

`placebo` is `NULL` or a list:

- `unit = TRUE` or a vector of ids
- `period = NULL` or a time
- `cut = Inf`

Orchestration: `spec → sc_fit → [placebo_time] → [placebo_space] → [loo]`. Attach pieces to the same `scm_fit`.

**Patterns:** `scm()` contains no math.

**Test scenarios:**
- [ ] Happy path: `scm(... smoking ..., nested = FALSE)` returns `scm_fit` with `att` numeric
- [ ] Happy path: `loo = TRUE` on a tiny toy attaches `$loo`
- [ ] Nil/empty: `placebo = list()` is treated as no placebo
- [ ] Error path: missing `trperiod` → error
- [ ] Edge: `allopt = TRUE, nested = FALSE` warns and runs allopt

**Verification:** A 15-line smoking example in the test (regression method, no placebo) finishes quickly and has 38 donors.

**Planning-time unknowns:** none.

- [ ] **Step 1:** Tests for mapping of flags and attachment of optional slots.
- [ ] **Step 2:** Implement `scm()`.
- [ ] **Step 3:** Commit `feat: add scm() workflow entry point`

---

### Unit 10: print, summary, extractors

**Goal:** Users can read the object without `str()`.

**Requirements trace:** R12, R9

**Dependencies:** Unit 9

**Files:**
- `R/sc_result.R`
- `tests/testthat/test-sc-result.R`

**Approach:** S3, not S7 (one class is enough).

- `print.scm_fit` — treated, time, method, RMSE, R², ATT, n donors
- `summary.scm_fit` — print + balance table + nonzero \(W\) + MSPE table if placebo present
- `coef.scm_fit` — \(W\)
- `v_weights(x)`, `balance(x)`, `effects(x)`, `glance(x)` (one-row tibble/data.frame of scalars)

No `broom` dependency. `glance` is our function.

**Patterns:** Methods do not recompute estimates.

**Test scenarios:**
- [ ] Happy path: `print` / `summary` return the object invisibly and emit RMSE
- [ ] Nil/empty: `summary` without placebo does not mention p-values
- [ ] Error path: `v_weights(1)` → error
- [ ] Edge: `coef(fit)` names match donor ids

**Verification:** `capture.output(summary(fit))` contains `"RMSE"` and `"Utah"` (or the numeric id mapped to Utah) on a smoking regression fit.

**Planning-time unknowns:** none.

- [ ] **Step 1:** Tests with `capture.output`.
- [ ] **Step 2:** Implement methods.
- [ ] **Step 3:** Commit `feat: print, summary, and extractors for scm_fit`

---

### Unit 11: Plots

**Goal:** ggplot2 charts covering the synth2 figure set the user actually needs.

**Requirements trace:** R10

**Dependencies:** Unit 10 (can depend on Unit 6 only; prefer Unit 9 so placebo plots have data)

**Files:**
- `R/sc_plot.R`
- `tests/testthat/test-sc-plot.R`

**Approach:**

```r
autoplot.scm_fit(object, type = c("path", "effect", "v", "w",
                                  "gaps", "mspe", "pvalue", "loo"), ...)
plot.scm_fit <- function(x, ...) print(autoplot(x, ...))
```

| type | content | requires |
|------|---------|----------|
| `path` | treated vs synthetic, vertical line at `trperiod` | always |
| `effect` | gap series | always |
| `v` | bar of \(V\) | always |
| `w` | bar of nonzero \(W\) | always |
| `gaps` | spaghetti of placebo gaps + treated | `placebo_space` |
| `mspe` | bar of ratios, treated highlighted | `placebo_space` |
| `pvalue` | year-wise left p | `placebo_space` |
| `loo` | original effect + min/max ribbon | `loo` |

Missing optional data → error `synthaio_plot_missing` naming the slot. Theme: `ggplot2::theme_minimal()` plus a single treated color `#B45C3D` and control grey `#8A8A8A`. No custom ggproto.

**Patterns:** Plot functions do not refit.

**Test scenarios:**
- [ ] Happy path: each always-available type returns a ggplot
- [ ] Nil/empty: `type = "gaps"` without placebo → `synthaio_plot_missing`
- [ ] Error path: `type = "nope"` → error
- [ ] Edge: `w` plot omits exact zeros

**Verification:** Tests use `ggplot2::is_ggplot()` (or `inherits(p, "ggplot")`). Do not `ggsave` in tests.

**Planning-time unknowns:** none.

- [ ] **Step 1:** Tests for class and error.
- [ ] **Step 2:** Implement eight types (optional types still implemented, fail if slot missing).
- [ ] **Step 3:** Commit `feat: ggplot2 autoplot methods for scm_fit`

---

### Unit 12: Gold-standard alignment and `R CMD check`

**Goal:** Prove R8 on smoking and leave the package CRAN-check clean.

**Requirements trace:** R8, R16, R17

**Dependencies:** Units 5–11

**Files:**
- `tests/testthat/test-gold-standard.R`
- `tests/testthat/test-gold-placebo.R` (skip_on_cran)
- `.github/workflows/R-CMD-check.yaml` (if we add CI; optional but recommended)
- `R/synthaio-package.R` / NAMESPACE — export the public API: `scm`, `sc_spec`, `sc_fit`, `v_weights`, `balance`, `effects`, `glance`, and the S3 methods `print.scm_fit`, `summary.scm_fit`, `coef.scm_fit`, `plot.scm_fit`, `autoplot.scm_fit` (each needs `@export` or `print()` / `plot()` will not dispatch)

**Approach:** Load `gold_smoking()`. Run the **full fixture `_spec`**, not a shorthand:

```r
scm(cigsale ~ lnincome + age15to24 + retprice + beer +
      cigsale(1988) + cigsale(1980) + cigsale(1975),
    data = smoking, unit = "state", time = "year",
    treat = california_id, trperiod = 1989,
    xperiod = 1980:1988, allopt = TRUE)
```

Cache that fit in a helper. Assert:

1. `sort(names(w[w > tol]))` equals `sort(names(gold$w_nonzero))` after mapping ids ↔ names.
2. Each nonzero \(W\) within `tol$w_abs`.
3. Top-2 \(V\) names match; each \(V\) within `tol$v_abs`.
4. RMSE relative error ≤ `tol$rmse_rel`; ATT within `tol$att_abs`.
5. Selected years' effects within `0.15` packs (looser than ATT; nested local minima).
6. Placebo test (skip_on_cran): California ratio ranks first; unfiltered p within `1e-9` of `1/n_units`. Do **not** assert `pre_mspe == 3.1668` on the allopt fit — that JSON block is the other run.

If (1) fails: **do not weaken the test**. Go back to Units 4–5 (scaling, starts, `sigf`, optional `quadprog` polish). Weakening tolerances requires editing the fixture `_source` and a plan amendment, not a silent `abs(x-y) < 1`.

**Patterns:** One source of truth: the JSON.

**Test scenarios:**
- [ ] Happy path: gold assertions above
- [ ] Error path: none (this is an integration test)
- [ ] Edge: `skip_on_cran()` on the 38-fit placebo loop
- [ ] `R CMD check --as-cran --no-manual` : 0 ERROR, 0 WARNING

**Verification:** `devtools::test(filter = "gold-standard")` passes locally. `devtools::check()` clean.

**Planning-time unknowns:** Whether allopt on this machine hits the same local min as Stata. **If blocked:** keep the donor-set assertion hard; relax individual \(W\) to `max(0.001, 0.05*|w|)` only after writing the achieved numbers into a `*_achieved.json` **and** opening a comment in the fixture that Stata log is still missing (Q1). Do not delete the gold file.

- [ ] **Step 1:** Write `test-gold-standard.R` from the JSON (fails if solver is off).
- [ ] **Step 2:** Tune Units 4–5 until it passes or the blocked protocol above is followed.
- [ ] **Step 3:** `devtools::document()`; export surface; `devtools::check()`.
- [ ] **Step 4:** Commit `test: lock smoking nested-allopt against synth2 fixtures`

---

### Unit 13: Vignettes

**Goal:** English and Chinese Prop 99 walkthroughs that a 连享会 reader can follow.

**Requirements trace:** R15

**Dependencies:** Unit 12 (vignettes must not contradict passing gold tests)

**Files:**
- `vignettes/prop99.Rmd`
- `vignettes/prop99-zh.Rmd`
- `README.md`
- `DESCRIPTION` — `VignetteBuilder: knitr`

**Approach:** Each vignette: load `smoking`, `scm(..., allopt = TRUE)`, print balance and nonzero \(W\), `autoplot` path + effect. Do **not** run full in-space placebo in the vignette by default (precompute or `eval=FALSE` with saved figures, or run with `cut` and `skip` on CRAN via `\donttest` style — for vignettes, use `allopt = TRUE` for the point estimate only, and show placebo code in a chunk with `eval = FALSE` plus a sentence that `tests/testthat/test-gold-placebo.R` is the checked run).

README: what the package is, what it is not (not gsynth), install, 10-line example.

**Patterns:** Numbers in prose must come from `fit` (`sprintf`), not hand-typed, so they cannot drift.

**Test scenarios:**
- [ ] Happy path: `devtools::build_vignettes()` succeeds
- [ ] Edge: Chinese vignette is valid UTF-8
- [ ] README example is copy-pasteable

**Verification:** `R CMD check` still 0 WARNING with vignettes (use `--no-build-vignettes` on CI if needed, but local build must work).

**Planning-time unknowns:** none.

- [ ] **Step 1:** Write both vignettes + README.
- [ ] **Step 2:** Build vignettes.
- [ ] **Step 3:** Commit `docs: add Prop 99 vignettes and README`

---

## Public API (do not grow in v1)

```r
scm()
sc_spec()
sc_fit()
coef.scm_fit()      # W
v_weights()
balance()
effects()
glance()
print.scm_fit()
summary.scm_fit()
autoplot.scm_fit()
plot.scm_fit()
```

Internals (`sc_solve_w`, `sc_solve_v`, `sc_placebo_*`, `sc_loo`) may be exported later; v1 keeps them unexported unless a test needs them (`::: ` is fine in tests).

---

## Quality bar

- [x] Every unit traces to requirement IDs
- [x] Dependencies are a DAG
- [x] Every unit has ≥ 3 test scenarios
- [x] No unit lists > 8 files
- [x] No unit introduces > 2 new abstractions (`sc_spec`, `scm_fit` are the only product types; solvers are functions)
- [x] Unknowns classified (Q1 Stata log deferred; regression-\(V\) formula locked to Synth OLS-on-predictors-with-intercept; `mspeperiod` default locked independent of `xperiod`; QP vs LOQO deferred to Unit 12 protocol)
- [x] Handoff: an engineer should not invent product behavior — formula grammar, Fisher p, R² window, scaling, mixed composition, and gold numbers are specified

---

## Out of scope (do not implement while executing this plan)

gsynth, SDID, scpi, augsynth, multi-treated joint estimation, Shiny, wrapping `tidysynth`.

---

## Execution notes

- TDD each unit: failing test → implement → pass → commit (messages given).
- Skill for implementation: `@superpowers-ruby:test-driven-development` and `@superpowers-ruby:verification-before-completion`.
- If Unit 5 donor-set smoke fails, do not start Unit 6.
- If Unit 12 gold fails, stay in Units 4–5.
