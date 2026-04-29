*! Capture golden-baseline outputs from KRLS Stata 1.01 (frozen reference).
*!
*! Runs a deterministic example through the source files in the repo root
*! and captures the key returns to dev/baseline_1.01.txt. Phase 2 fixes
*! must reproduce this byte-for-byte.

clear all
set more off
set seed 20260429

adopath ++ "`c(pwd)'"
which krls
which kpredict

use "growthdata.dta", clear

log using "dev/baseline_1.01.log", replace text

* canonical Barro-style growth example: regress growth on initial GDP,
* trade share, schooling, and assassinations.
krls growth rgdp60 tradeshare yearsschool assassinations, deriv vcov

log close

file open out using "dev/baseline_1.01.txt", write replace text
file write out "KRLS Stata baseline -- frozen 1.01 reference" _n _n

file write out "---- e(R2) ----" _n
file write out %18.10f (e(R2)) _n _n

file write out "---- e(lambda) ----" _n
file write out %18.10f (e(lambda)) _n _n

file write out "---- e(sigma) ----" _n
file write out %18.10f (e(sigma)) _n _n

file write out "---- e(Looloss) ----" _n
file write out %18.10f (e(Looloss)) _n _n

file write out "---- e(Eff_Degrees) ----" _n
file write out %18.10f (e(Eff_Degrees)) _n _n

file write out "---- e(Output) — pointwise derivatives summary ----" _n
mat O = e(Output)
forv i = 1/`=rowsof(O)' {
  forv j = 1/`=colsof(O)' {
    file write out %18.10f (O[`i',`j']) " "
  }
  file write out _n
}

file close out

display "Wrote dev/baseline_1.01.txt"
