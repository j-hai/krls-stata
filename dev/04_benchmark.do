*! Quick wall-clock benchmark for krls at varying n.

clear all
set more off
adopath ++ "`c(pwd)'"

local seed 20260429

quietly {
  set obs 1000
  set seed `seed'
  gen x1 = rnormal()
  gen x2 = rnormal()
  gen x3 = rnormal()
  gen x4 = rnormal()
  gen y  = x1 + 0.5*x2 + 0.3*x3^2 + rnormal(0, 0.3)
}

foreach N in 100 200 400 600 {
  preserve
    qui keep in 1/`N'
    timer clear 1
    timer on 1
    qui krls y x1 x2 x3 x4, deriv vcov
    timer off 1
    qui timer list 1
    di as txt "N=`N': " as res %9.3f r(t1) as txt " seconds"
  restore
}
