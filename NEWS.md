# KRLS Stata 1.03

## Performance

* Vectorized the three pairwise-distance helpers in `krls.ado` and
  `kpredict.ado`:
  - `m_euclidian_distance` (used to build the kernel matrix `K`)
  - `m_distance` (used to build pointwise derivatives)
  - the identical helper inside `kpredict.ado`

  Each was a manual `for (i; for (j ...))` Mata loop costing
  `O(n^2 * d)` interpreted iterations. The new code uses the BLAS
  identity `||x_i - x_j||^2 = ||x_i||^2 + ||x_j||^2 - 2 x_i' x_j`,
  collapsing the inner work into a single `X * X'` matrix
  multiplication. Outer-difference for derivatives is similarly one
  broadcast.

* Skipped a redundant `sqrt(...)^2` round-trip — both call sites of
  `m_euclidian_distance` only needed the squared distance to feed
  `exp(-d^2 / sigma)`, but the old code computed `sqrt`, returned,
  then squared back. The new internal helper
  `m_euclidian_distance_sq` short-circuits this. The original
  `m_euclidian_distance` is preserved as a wrapper for any external
  caller.

Benchmark (Apple Silicon Stata 17 MP, simulated 4-predictor data):

  N  =  100:  27 ms ->  21 ms   (1.29x)
  N  =  200:  86 ms ->  65 ms   (1.32x)
  N  =  400: 371 ms -> 295 ms   (1.26x)
  N  =  600:1168 ms -> 993 ms   (1.18x)

Numerical results unchanged — verified byte-for-byte against the
1.01 baseline on the canonical Barro growth example.

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
