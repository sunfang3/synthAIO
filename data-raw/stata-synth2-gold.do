* Gold-standard dump for synthaio (Abadie Prop 99 / smoking).
* Run this in Stata on the Windows host, not in WSL.
*
* First time only (Stata command window):
*   ssc install synth, replace
*   ssc install synth2, replace
*
* Working directory MUST be writable. Change the next line if needed.
* Example:  cd "C:\Users\YOURNAME\Desktop"

clear all
set more off

* Text log (readable in WSL). SMCL is Stata's default and harder to grep.
capture log close
log using "synth2_smoking_allopt.log", replace text

di as txt "==== which synth2 ===="
which synth2
which synth

* smoking.dta is NOT a Stata official sysuse file. Do not use: sysuse smoking
*
* Easiest: copy the repo CSV to the Desktop first (from WSL):
*   cp /home/fang/Project/synthAIO/data-raw/smoking.csv /mnt/c/Users/你的用户名/Desktop/
*
* Then set this path (edit the username):
local csv `"C:\Users\你的用户名\Desktop\smoking.csv"'
*
* If you already copied it next to this .do, you can instead:
* local csv `"smoking.csv"'

capture confirm file `"`csv'"'
if _rc {
    di as err "CSV not found: `csv'"
    di as err "Copy data-raw/smoking.csv to the Desktop and edit the local csv path."
    exit 601
}

import delimited "`csv'", clear varnames(1)
destring state year cigsale lnincome age15to24 retprice beer, replace force
xtset state year
describe
summarize
count
tab state if year == 1988

di as txt "==== BASELINE: nested allopt ===="
synth2 cigsale lnincome age15to24 retprice beer ///
    cigsale(1988) cigsale(1980) cigsale(1975), ///
    trunit(3) trperiod(1989) xperiod(1980(1)1988) ///
    nested allopt

di as txt "==== ereturn scalars ===="
ereturn list

di as txt "==== V weights ===="
matrix list e(V_wt), format(%12.8f)

di as txt "==== unit weights ===="
matrix list e(U_wt), format(%12.8f)

di as txt "==== balance ===="
matrix list e(bal), format(%12.6f)

di as txt "==== key scalars ===="
di "rmse = " e(rmse)
di "r2   = " e(r2)
di "att  = " e(att)

log close

di as txt "Log written to:"
pwd
di as txt "synth2_smoking_allopt.log"

* Optional, slow (~minutes). Uncomment if you also want placebo p-values.
* log using "synth2_smoking_placebo.log", replace text
* synth2 cigsale lnincome age15to24 retprice beer ///
*     cigsale(1988) cigsale(1980) cigsale(1975), ///
*     trunit(3) trperiod(1989) xperiod(1980(1)1988) ///
*     nested placebo(unit cut(2)) sigf(6)
* matrix list e(mspe), format(%12.6f)
* matrix list e(pval), format(%12.6f)
* log close
