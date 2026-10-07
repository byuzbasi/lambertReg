# Validation and contributing

## What is checked

The package’s local validation compares the public implementation
against the preserved reference engine on small scaled, correlated,
high-dimensional, and budget-failure examples. Checks cover coefficient
paths, fold losses, selected penalties, support counts, both starting
points, preprocessing, numerical status, and original-scale predictions.
Plot checks cover stored-fit displays and gaps for unavailable fits.

The latest local validation for version 0.1.10 passed 142 targeted
checks and `R CMD check` with status OK. These are development checks,
not CRAN acceptance or a guarantee of a global minimum. The
documentation examples are also executed during the local site build.

## Interpret numerical status

The objective is nonconvex. A returned fit is a stationary solution that
meets the package’s numerical acceptance rules. The solver considers
warm and zero starts, checks normalized stationarity, and reports
active-curvature diagnostics. These checks do not certify global
optimality.

- `ok` describes the selected fit’s availability and numerical
  acceptance.
- `partial_search` records unavailable CV candidates.
- `full_path_complete` describes the stored full-training path.
- `cv_diagnostics` and path `diagnostics` provide candidate-level
  details.

An unavailable candidate is ineligible for CV selection. Predictions and
coefficient extraction reject an unverified selected fit. Plotting
presents unavailable path positions as gaps.

## Scope of the interface

The package supports dense finite matrices and Gaussian responses. It
does not currently support observation weights, automatic missing-value
handling, or a high-dimensional oracle guarantee. The shape is fixed and
early stopping is not applied. Consult the function reference for
defaults and controls.

## Report a problem

Open an issue at the [source
repository](https://github.com/byuzbasi/lambertReg/issues). Include a
minimal synthetic or public-data example, seed and folds when relevant,
package version,
[`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html), controls,
warnings, and fit status. For plots, include the device size. Do not
attach confidential data.

The repository’s [contribution
guide](https://github.com/byuzbasi/lambertReg/blob/main/CONTRIBUTING.md)
describes local checks and the review of methodological or numerical
changes.
