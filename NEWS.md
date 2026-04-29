# KRLS Stata 1.02

## Bug fixes

* **`krls.ado`**: removed a stray trailing apostrophe on the
  `local vcov = "vcov"` line (was `local vcov = "vcov"'` — orphaned
  closing-tick from a deleted macro reference). The `svcov()` option
  path is now syntactically clean.
* **`kpredict.ado`**: fixed the `syntax` declaration to space-separate
  options (`Se Residuals`) instead of comma-separate
  (`Se, Fitted, Residuals`). Removed the unused `Fitted` option that
  was declared but never referenced in the body.
* **`krls.ado`**: removed the unused `[DOUBLE]` option that was
  declared on the `syntax` line but never referenced.
* **`krls.ado`**: the "Derivatives suppressed" status string was a
  bare Mata expression that evaluated and discarded; it now actually
  prints via `printf` (and the typo "Supressed" is fixed).

## Stata version declaration

* Bumped `version 11` → `version 13` in both `krls.ado` and
  `kpredict.ado` (in both the program body and the Mata block).
  Stata 11 was 2009; Stata 13 (2013) is well below the user-base
  floor and removes a 16-year-old compatibility layer.

## Verified

The canonical Barro growth example reproduces byte-for-byte against
the 1.01 baseline:

* `e(R2)` = 0.5237911523
* `e(lambda)` = 0.4805161997
* `e(sigma)` = 4.0
* `e(Looloss)` = 97.4977706922
* `e(Eff_Degrees)` = 16.1742104324
* `e(Output)` (pointwise derivatives matrix): unchanged

5 integration tests in `dev/03_integration_tests.do` cover default,
`suppress`, user-supplied `lambda`, user-supplied `sigma`, and a
`kpredict` round-trip — all pass.
