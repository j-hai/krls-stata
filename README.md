# KRLS for Stata

Stata implementation of **Kernel-Based Regularized Least Squares (KRLS)**
— a machine-learning method for regression and classification that fits
a flexible function `y = f(x)` without assuming linearity or additivity.
Pure Stata/Mata implementation (no compiled plugin).

The method is described in:

> Hainmueller, J., & Hazlett, C. (2014). "Kernel Regularized Least Squares: Reducing Misspecification Bias with a Flexible and Interpretable Machine Learning Approach." *Political Analysis*, 22(2), 143–168.
>
> Ferwerda, J., Hainmueller, J., & Hazlett, C. (2017). "Kernel-Based Regularized Least Squares in R (KRLS) and Stata (krls)." *Journal of Statistical Software*, 79(3).

## Installation

```stata
* From SSC (stable):
ssc install krls, replace

* Development version from GitHub:
net install krls, from(https://raw.githubusercontent.com/j-hai/krls-stata/main/k/) replace
```

When net-installing directly from GitHub, Stata installs the command
and help files, but not the bundled example dataset. To download
`growthdata.dta` into your current working directory, run:

```stata
net get krls, from(https://raw.githubusercontent.com/j-hai/krls-stata/main/k/)
```

## Quick start

```stata
use https://raw.githubusercontent.com/j-hai/krls-stata/main/g/growthdata.dta, clear

* Fit
krls growth rgdp60 tradeshare yearsschool assassinations, deriv vcov

* Predict on the same sample
kpredict yhat

* Predicted values' standard errors
kpredict yhat_se, se

* Residuals
kpredict resid, residuals
```

After `krls`, the standard ereturned values are available:

* `e(R2)` — coefficient of determination
* `e(lambda)` — chosen regularization parameter
* `e(sigma)` — kernel bandwidth (default = number of predictors)
* `e(Looloss)` — leave-one-out loss
* `e(Eff_Degrees)` — effective degrees of freedom
* `e(Output)` — table of average and quantile pointwise derivatives
* `e(Varcovfit)` — variance-covariance of fitted values (with `vcov` option)

## What's new in 1.02

* Several bug fixes — syntax glitch in `krls.ado` (stray apostrophe in
  the `svcov` branch), comma-separated options in `kpredict.ado`'s
  syntax line, dead `[DOUBLE]` option in `krls.ado`'s syntax,
  "Supressed" typo replaced with a real `printf`.
* `version` declaration bumped 11 → 13.
* No numerical changes — verified byte-for-byte against the 1.01
  baseline on the canonical Barro growth example.

See [`NEWS.md`](NEWS.md) for the full change log.

## License

GPL (>= 2).
