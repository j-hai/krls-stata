*! Integration tests for krls — broader coverage than baseline.
*!
*! Exercises the main option paths so a refactor that breaks one of them
*! is caught even if the basic growth example still reproduces:
*!
*!   t1: krls default (deriv + vcov, default sigma/lambda)
*!   t2: krls suppress (no derivatives reported)
*!   t3: krls custom lambda (skip optimization)
*!   t4: krls custom sigma
*!   t5: kpredict on the fitted sample (round-trip)
*!
*! How to run from repo root:
*!   /Applications/Stata/StataMP.app/Contents/MacOS/stata-mp -b do dev/03_integration_tests.do

clear all
set more off
set seed 20260429

adopath ++ "`c(pwd)'"

log using "dev/integration_tests.log", replace text

file open out using "dev/integration_outputs.txt", write replace text
file write out "KRLS Stata integration tests — current source" _n _n

capture program drop dump
program dump
  args label
  file write out "==== `label' ====" _n
  file write out "R2:        " %14.10f (e(R2)) _n
  file write out "lambda:    " %14.10f (e(lambda)) _n
  file write out "sigma:     " %14.10f (e(sigma)) _n
  file write out "Looloss:   " %14.10f (e(Looloss)) _n
  file write out "Eff_Deg:   " %14.10f (e(Eff_Degrees)) _n _n
end

* -------------------- t1: krls default --------------------
display _newline as txt "=== t1: krls default ==="
use "growthdata.dta", clear
krls growth rgdp60 tradeshare yearsschool assassinations, deriv vcov
dump "t1 default"

* -------------------- t2: krls suppress --------------------
display _newline as txt "=== t2: krls suppress ==="
use "growthdata.dta", clear
krls growth rgdp60 tradeshare yearsschool assassinations, suppress
dump "t2 suppress"

* -------------------- t3: krls user lambda --------------------
display _newline as txt "=== t3: krls user lambda=0.5 ==="
use "growthdata.dta", clear
krls growth rgdp60 tradeshare yearsschool assassinations, lambda(0.5) deriv vcov
dump "t3 user lambda"

* -------------------- t4: krls user sigma --------------------
display _newline as txt "=== t4: krls user sigma=2 ==="
use "growthdata.dta", clear
krls growth rgdp60 tradeshare yearsschool assassinations, sigma(2) deriv vcov
dump "t4 user sigma"

* -------------------- t5: kpredict round-trip --------------------
display _newline as txt "=== t5: kpredict round-trip ==="
use "growthdata.dta", clear
krls growth rgdp60 tradeshare yearsschool assassinations, deriv vcov
kpredict yhat
sum yhat growth
file write out "==== t5 kpredict ====" _n
qui sum yhat
file write out "yhat mean:    " %14.10f (r(mean)) _n
qui sum growth
file write out "growth mean:  " %14.10f (r(mean)) _n
qui corr yhat growth
file write out "cor(yhat, growth):  " %14.10f (r(rho)) _n _n

file close out
log close

display "Wrote dev/integration_outputs.txt and dev/integration_tests.log"
