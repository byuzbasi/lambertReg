# Training-Only Cross-Validation for Lambert Regression

Choose a relative penalty-path position by the arithmetic mean of
validation-fold squared errors, with fixed Lambert shape and
fold-specific preprocessing.

## Usage

``` r
cv.lambert(x, y, nfolds = 5L, foldid = NULL,
           lambda_fraction = 10^seq(0, -3, length.out = 60L),
           max_sweeps = 10000L, kkt_tol = 1e-7)
```

## Arguments

- x, y:

  Training predictors and response, as in [`lambert`](lambert.md).

- nfolds:

  Number of folds, default five. Each training subset must contain at
  least three observations. Inferred from `foldid` when that argument is
  supplied and `nfolds` is omitted.

- foldid:

  Optional integer fold label for every row, with consecutive labels
  starting at one and no empty folds. Store and reuse this vector for
  matched comparisons. Without it, balanced labels are randomly permuted
  using the caller's R random-number state. Use
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for
  reproducibility.

- lambda_fraction:

  Strictly decreasing fractions in \\(0,1\]\\; default 60
  logarithmically spaced positions from 1 to 0.001.

- max_sweeps, kkt_tol:

  Solver controls described in [`lambert`](lambert.md).

## Details

Each training fold estimates its own predictor centers, RMS scales,
response mean and entry score. Candidate \\k\\ in fold \\f\\ uses
\\\lambda\_{k,f}=\lambda\_{\max,f}\tau_k\\, where \\\tau_k\\ is the
supplied fraction. The candidate is aligned by fraction, not by a common
absolute lambda. The validation observations are never used to estimate
these values.

For validation sets \\I_f\\ and training-fold predictions \\\hat y_i\\,
the score is \$\$\mathrm{CV}(k)=\frac{1}{K}\sum\_{f=1}^K
\frac{1}{\|I_f\|}\sum\_{i\in I_f}(y_i-\hat y_i)^2.\$\$ This is an
unweighted mean of the fold MSEs, including when fold sizes differ. Only
candidates with finite loss in every fold are eligible. Ties within
\\10^{-10}(1+\mathrm{CV}\_{\min})\\ favor the larger penalty fraction.
`lambda.min` is the selected fraction multiplied by the full training
entry score. A full-training path is fitted at every supplied fraction.
Later positions provide variable counts for plotting and do not affect
CV selection. The selected full and fold fits receive the reference
active-curvature check. A rejected selected fit is reported explicitly;
another candidate is not silently substituted.

An incomplete search is flagged by `partial_search`; it does not imply
that the selected fit failed. Solver convergence and curvature status
must be inspected separately. No one-standard-error rule or optional
early pruning is implemented in this initial package. CV error is a
tuning score, not an unbiased external performance estimate; use
independent test data or outer cross-validation to evaluate prediction.
Changing the number of folds or numerical controls is supported but
departs from the study defaults.

## Value

A list of class `cv.lambert`. Write \\L\\ for the number of fractions,
\\K\\ for the number of folds, and \\n\\ for the number of input rows.
Components are:

- cvm:

  Length-\\L\\ arithmetic means of fold MSEs. An ineligible candidate
  has value `Inf`.

- fold_loss:

  Numeric \\L\times K\\ validation-MSE matrix, with `NA` where a fold
  candidate is unavailable.

- foldid, nfolds:

  Length-\\n\\ fold labels and fold count.

- index:

  Selected path index, or `NA_integer_` when none is eligible.

- lambda_fraction, lambda:

  Length-\\L\\ fractions and corresponding full-training penalty levels.
  These absolute levels need not equal fold levels.

- fraction.min, lambda.min:

  Selected fraction and full-training penalty; `NA_real_` if no
  candidate is eligible.

- lambert.fit:

  Full-training `lambert` path, including when no CV candidate is
  eligible. Its fields are documented in [`lambert`](lambert.md).

- nzero, support_tol:

  Length-\\L\\ selected-variable counts and their threshold. Counts
  exclude the intercept and use original-scale coefficient magnitudes
  greater than \\10^{-8}\\. Unverified counts are `NA`.

- full_path_complete:

  Whether every full-training position is verified; this is separate
  from fold-search completeness and selected-fit validity.

- ok, status:

  Logical selected-fit acceptance and character status. `no_valid_cv`
  means no candidate had finite loss in all folds. Other status values
  follow [`lambert`](lambert.md). A negative-curvature selected fold can
  set `status = "negative_curvature"`.

- partial_search:

  Whether any candidate has a nonfinite CV score. A partial search can
  still have `ok = TRUE`.

- cv_diagnostics:

  Fold-stacked solver diagnostics, with additional `fold` and validation
  `loss` columns.

- attempts:

  Fold-stacked records for both starts, with a `fold` column. The
  `executed` and `reused` flags are described in
  [`lambert`](lambert.md). Full-training attempts are in
  `lambert.fit$attempts`.

- selected_curvature:

  Length-\\K\\ list of selected-fold Hessian assessments; `NULL` if no
  CV candidate is eligible. Full-training assessments are in
  `lambert.fit$curvature`.

- preprocessing:

  Length-\\K+1\\ list: each training fold, then the full training
  sample. Each entry contains `mx`, `my`, `sx` and `lambda0`, the
  predictor means, response mean, RMS scales and entry score.

- control, shape, call:

  Numerical settings, fixed shape (one), and call.

Completing the full-training path can cost more than fitting only its
selected prefix. Plotting reads stored values and never fits models.

## See also

[`lambert`](lambert.md), [`predict.lambert`](methods.md)

## Examples

``` r
set.seed(12)
x <- matrix(rnorm(60 * 6), 60, 6)
y <- 2 * x[, 1] - x[, 2] + rnorm(60)
foldid <- sample(rep(1:5, length.out = 60))
fit <- cv.lambert(x, y, foldid = foldid,
                  lambda_fraction = c(1, 0.5, 0.2))
print(fit)
#> Cross-validated Lambert regression (5 folds)
#> Selected index: 3; full-training lambda: 0.3597978
#> Status: verified_stationary; partial search: FALSE
coef(fit)
#>                     s3
#> (Intercept)  0.1123142
#> V1           1.9992149
#> V2          -1.1984651
#> V3           0.0000000
#> V4           0.0000000
#> V5           0.0000000
#> V6           0.0000000
predict(fit, x[1:3, , drop = FALSE])
#>             s3
#> [1,] -3.334891
#> [2,]  2.073636
#> [3,] -2.826032
```
