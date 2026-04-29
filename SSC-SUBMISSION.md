# Submitting krls 1.02 to SSC

The Stata `krls` package is distributed via SSC. To update the SSC copy
with this 1.02 release, email Kit Baum at Boston College.

## Files to send

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
Subject: `krls update — version 1.02`

```
Dear Kit,

Please find attached an updated version (1.02) of the krls package
for SSC. This is a small bug-fix release relative to 1.01 (Dec 2013):

* krls.ado: removed a stray trailing apostrophe on
  'local vcov = "vcov"' (line 54) — orphaned closing-tick from a
  deleted macro reference. The svcov() option path is now
  syntactically clean.
* kpredict.ado: fixed the syntax declaration to space-separate
  options (Se Residuals) instead of comma-separate
  (Se, Fitted, Residuals), and dropped the unused 'Fitted' option
  that was declared but never referenced.
* krls.ado: removed the unused [DOUBLE] option declared on the
  syntax line but never referenced.
* The 'Derivatives suppressed' status string was a bare Mata
  expression that evaluated and discarded; it now actually prints
  via printf (and the typo "Supressed" is fixed).
* version declaration bumped 11 -> 13 in both .ado files.

No numerical changes — verified byte-for-byte against the 1.01
baseline on the canonical Barro growth example
(R2 = 0.5237911523, lambda = 0.4805161997, etc.).

Source repository (with full change log and test harness):
https://github.com/j-hai/krls-stata

Best,
Jens Hainmueller
jhain@stanford.edu
```

## How to package the zip

```sh
cd /Users/jhainmueller/Documents/GitHub/krls-stata
zip ~/Desktop/krls-1.02-ssc.zip \
  krls.ado krls.sthlp kpredict.ado kpredict.sthlp krls.pkg
```

## Test users can do before SSC is live

```stata
net install krls, from(https://raw.githubusercontent.com/j-hai/krls-stata/main/) replace
```
