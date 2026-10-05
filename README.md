# lambertReg

Sparse Gaussian regression with the **Lambert** penalty. Shape is fixed at
c = 1; cross-validation selects only the relative penalty-path position.
This initial package uses the manuscript's reference R/C++ implementation.
It is a development package hosted in a private GitHub repository; it is not
a CRAN release. Repository access requires authorization.

## Installation

Install the source archive with `R CMD INSTALL lambertReg_0.1.9.tar.gz`.
R (>= 4.1.0), Rcpp, and an R-compatible C++ toolchain are required.
The fitting functions do not require glmnet, ncvreg, mltools, knitr or
kableExtra. Building the vignette additionally uses knitr and rmarkdown
under Suggests. A source installation requires Rtools on Windows or the
corresponding R-compatible developer tools on macOS/Linux.

Authorized GitHub users can install from a local clone:

```sh
gh repo clone byuzbasi/lambertReg
R CMD INSTALL lambertReg
```

Use your existing GitHub authentication; do not place access tokens in scripts.
To build the installed vignette from a clone, with knitr, rmarkdown and Pandoc
available, run `R CMD build lambertReg` and install the resulting archive.
The checked source archive already includes the built vignette.

## Example

```r
library(lambertReg)
set.seed(42)
x <- matrix(rnorm(100 * 8), 100, 8)
y <- 2 * x[, 1] - x[, 2] + rnorm(100)
foldid <- sample(rep(1:5, length.out = nrow(x)))
fit <- cv.lambert(x, y, foldid = foldid)
print(fit)
coef(fit)
predict(fit, newx = x[1:3, , drop = FALSE])
plot(fit)
```

These predictions only illustrate the API. Evaluate performance on separate
observations or with outer cross-validation. CV tuning error is not external
test error.

## Matching the manuscript

- Fixed c = 1, 60 descending penalty fractions from 1 to 0.001.
- Training-only predictor centering/RMS scaling and response centering.
- Unpenalized intercept and predictions returned on the original scale.
- Five-fold CV by default; arithmetic mean of fold MSEs, even for unequal folds.
- Each fold uses its own entry score. The grid aligns fractions rather than
  absolute lambda values across folds. `lambda.min` is the selected fraction
  times the full-training entry score.
- Both warm and zero starts; objective ties favor warm starts, CV ties favor
  larger penalty fractions. Relative tie tolerance is 1e-10.
- 10,000 sweeps per start and normalized KKT tolerance 1e-7 by default.
- Selected-fit active-curvature checks follow the reference implementation.
- No response variance scaling, post-selection refit, shape tuning, or
  coefficient interpolation. No optional early pruning in this initial version.

The package contains only the Lambert estimation components, not the study's
simulation runner, comparator methods, or numerical datasets. Source hashes
and extraction notes are in `system.file("PROVENANCE.json", package="lambertReg")`.

## Failure handling and scope

Inspect `fit$ok`, `fit$status`, `fit$partial_search`, and
`fit$cv_diagnostics`. An incomplete CV search differs from a failed selected
fit. Predict/coef refuse unverified fits. A stationary solution is not a
certificate of a global optimum. At curvature joins the Hessian diagnostic
can be unavailable; this is reported explicitly.

Inputs are finite dense numeric matrices. Missing values, constant training
columns, zero entry scores, and incompatible prediction columns cause clear
errors. No observations or predictors are silently removed. Preprocessing
inside each CV training fold can reject a column that varies only in its
held-out observations.

CV uses the caller's random-number state when `foldid` is absent. For
reproducibility, set a seed and save `foldid`. Supplied folds do not consume
random numbers. The package does not impose the simulation's master seed on
new applications.

## Citation and license

Use `citation("lambertReg")`. The associated manuscript is currently
unpublished: *The Lambert Penalty: Logarithmic Shrinkage for Sparse Regression*
(Bahadir Yuzbasi, 2026). No journal or DOI is claimed. Source repository:
<https://github.com/byuzbasi/lambertReg> (private). Report issues there when
authorized, or contact the maintainer by email.

GPL-3. Maintainer: Bahadir Yuzbasi <b.yzb@hotmail.com>.

## CV error plot (0.1.9)

```r
plot(fit)                    # detailed CV plot, error bars and selected-fit summary
plot(fit, xaxis = "lambda")  # log full-training lambda
plot(fit, style = "paper")    # compact CV plot for publication
plot(fit, show_nzero = FALSE)
```

Blue open circles and a connecting line show arithmetic mean validation-fold
MSE, with pale error bars and subtle horizontal guides. An orange marker and
dashed line identify the CV-selected candidate. The logarithmic horizontal
axis increases from left to right, from weaker to stronger regularization.
The default detailed view adds a title and a compact summary of the selected
lambda, CV MSE and variable count. Use `style = "paper"` to omit these annotations.
The upper axis counts
variables in each verified full-training model, excluding the intercept and
using the study threshold `abs(beta) > 1e-8` on the original coefficient scale.
Missing/unverified counts appear as `--`. Counts are not averages over CV folds.

`cv.lambert()` now computes and stores the complete full-training path. This
adds work after the selected position but leaves the CV scores, selection rule,
and selected coefficients unchanged. `plot()` never refits models. Saved
0.1.0 objects can still be plotted; unavailable counts are clearly marked.
For the default 60-point grid, use a wide plot (about 9 inches) to read every
upper label. Colors, title, axis labels and upper-label size can be customized.
Error bars show the mean plus/minus one fold-based standard error:
`sd(fold MSE) / sqrt(K)`. This is a descriptive display because training folds
overlap; it is not an independent-sample confidence interval. Every fold loss
must be finite to draw a bar. Unavailable bars are omitted, and the existing
minimum-CV selection rule is unchanged (no one-SE rule).

Upper counts are unrotated and sized to avoid overlap. Error bars are drawn
from stored fold losses, without refitting. Use `plot(fit, error_bars = FALSE)`
to hide the bars. The default y limits include their full extent.


## Coefficient paths

```r
plot(fit, type = "coefficients")
plot(fit, type = "coefficients", labels = TRUE, xaxis = "lambda")
plot(fit, type = "coefficients", style = "paper")
plot(fit$lambert.fit)         # also works for objects returned by lambert()
```

Each colored curve follows one original-scale coefficient, excluding the
intercept. Stronger regularization is on the left. The zero reference line
shows where coefficients enter or leave the model. From a CV object, a
vertical dashed line marks the selected penalty; a plain path object has
no automatic selection. `labels = TRUE` adds a predictor legend outside the
panel. Counts and coefficients do not describe causal importance.

Unverified fits produce gaps, never invented values or interpolated lines.
Only stored path positions are plotted; older saved objects may have shorter
paths. Both plotting styles use the same numerical values and never refit
models. Predictions, coefficient extraction, CV scores and tuning rules are
unchanged in this release.

## Help and worked guide

```r
help(package = "lambertReg")
?lambert
?cv.lambert
?predict.lambert
?lambert_penalty
vignette("introduction", package = "lambertReg")
citation("lambertReg")
```

The vignette uses a small reproducible synthetic example, a separate test set,
and saved CV folds. It describes every preprocessing operation, returned-object
layout, numerical status and plotting convention. It is a usage demonstration,
not a benchmark or an empirical result from the manuscript.

## Local checking

With the declared dependencies available, run `R CMD build lambertReg`, then
`R CMD check --as-cran lambertReg_0.1.9.tar.gz`. This builds and checks the help,
examples, tests and vignette. A private GitHub URL is inaccessible to anonymous
URL validators. Any such findings must be reported when preparing a public
submission. CRAN acceptance and cross-platform checks are separate steps.
No CI jobs or external check services are required by this repository.
