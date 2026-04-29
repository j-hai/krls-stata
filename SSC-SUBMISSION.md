# Submitting krls 1.03 to SSC

The Stata `krls` package is distributed via SSC. To update the SSC copy
with this 1.03 release, email Kit Baum at Boston College — but **don't
attach the zip** (Gmail/Google Workspace tends to block emails carrying
binaries; that bounced for the synth submission). Instead point Kit at
the GitHub release URL.

## Files in the SSC bundle

```
krls.ado
krls.sthlp
kpredict.ado
kpredict.sthlp
krls.pkg
```

(Don't ship the `.dta` example datasets to SSC — those are repo
illustrations only.)

## Email template

To: `kit.baum@bc.edu`
Subject: `krls update — version 1.03`

```
Dear Kit,

Please find an updated version (1.03) of the krls package for SSC,
downloadable from:

  https://github.com/j-hai/krls-stata/releases/download/v1.03/krls-1.03-ssc.zip

This rolls together two release notches relative to the 1.01 (Dec 2013)
copy currently on SSC:

1.02 — Bug fixes (no numerical change):
  * krls.ado: removed a stray trailing apostrophe on
    'local vcov = "vcov"' (line 54) — orphaned closing-tick from a
    deleted macro reference. svcov() option path now syntactically clean.
  * kpredict.ado: syntax declaration's options list was comma-separated
    (Se, Fitted, Residuals) instead of space-separated. Now (Se Residuals);
    the unused 'Fitted' option that was declared but never referenced
    is removed.
  * krls.ado: removed unused [DOUBLE] option declared on the syntax
    line but never referenced.
  * 'Derivatives suppressed' status string was a bare Mata expression;
    now uses printf, also fixes 'Supressed' typo.
  * version declaration bumped 11 -> 13 in both .ado files.

1.03 — Performance:
  * Vectorized the three pairwise-distance helpers via the BLAS identity
    ||x_i - x_j||^2 = ||x_i||^2 + ||x_j||^2 - 2 x_i' x_j. Each was a
    manual nested-loop in Mata; the rewrite collapses them to one
    matrix multiplication per call.
  * Skipped a redundant sqrt(...)^2 round-trip — both call sites only
    needed squared distance to feed exp(-d^2 / sigma).
  * Wall-clock measurements (Apple Silicon Stata 17 MP, 4-predictor
    synthetic data): 1.18-1.32x faster across n=100..600.

Numerical results unchanged — byte-for-byte against the 1.01 baseline
on the canonical Barro growth example (R2 = 0.5237911523,
lambda = 0.4805161997).

Source repository (with full change log, test harness, and benchmark
script):
  https://github.com/j-hai/krls-stata

Best,
Jens Hainmueller
jhain@stanford.edu
```

## How to (re-)package the zip

```sh
cd /Users/jhainmueller/Documents/GitHub/krls-stata
zip ~/Desktop/krls-1.03-ssc.zip \
  krls.ado krls.sthlp kpredict.ado kpredict.sthlp krls.pkg
```

The release on GitHub also has the zip attached at the URL above —
that's the canonical download for Kit and end users.

## Test users can do before SSC is live

```stata
net install krls, from(https://raw.githubusercontent.com/j-hai/krls-stata/main/) replace
```
