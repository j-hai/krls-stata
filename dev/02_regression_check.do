*! Regression check: current source vs. frozen 1.01 baseline.
*!
*! Re-runs the canonical Barro growth example and writes the same matrix
*! dump to dev/current.txt. Then `diff dev/baseline_1.01.txt dev/current.txt`
*! is the truth test — clean = byte-identical.

clear all
set more off
set seed 20260429

adopath ++ "`c(pwd)'"

use "growthdata.dta", clear

log using "dev/current.log", replace text
krls growth rgdp60 tradeshare yearsschool assassinations, deriv vcov
log close

file open out using "dev/current.txt", write replace text
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
display "Wrote dev/current.txt"
