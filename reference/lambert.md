# Gaussian Regression Paths with the Fixed Lambert Penalty

Fit a regularization path with the bounded Lambert penalty of shape
\\c=1\\. Exact cyclic scalar updates are computed by a compiled solver
that preserves the reference fitting rule. No shape selection or
post-selection refit is performed.

## Usage

``` r
lambert(x, y, lambda = NULL,
        lambda_fraction = 10^seq(0, -3, length.out = 60L),
        max_sweeps = 10000L, kkt_tol = 1e-7)
```

## Arguments

- x:

  Finite numeric predictor matrix, with at least three rows and one
  column. Columns must vary in every training sample. Optional column
  names must be nonempty and unique. Sparse matrices, missing values and
  automatic factor encoding are not supported.

- y:

  Finite numeric response vector. The response is centered, not scaled
  to unit variance. A positive entry score is required.

- lambda:

  Optional strictly decreasing vector of positive penalty levels on the
  standardized-predictor scale. Cannot be supplied together with
  `lambda_fraction`. Zero is not supported.

- lambda_fraction:

  Strictly decreasing fractions in \\(0,1\]\\, multiplied by the
  training entry score \\\lambda\_{\max}\\. Defaults to 60 positions
  from 1 to 0.001.

- max_sweeps:

  Positive integer sweep budget for each start and penalty level.

- kkt_tol:

  Positive normalized KKT acceptance tolerance. The internal solver uses
  half this tolerance, followed by an independent residual check.

## Details

For \\n\\ training observations, each predictor is centered and divided
by its root-mean-square centered value so that the standardized column
has squared norm \\n\\. The intercept is unpenalized. All centering and
scaling statistics are estimated from the supplied training
observations. The objective is \$\$\\y_c-Z\beta\\\_2^2/(2n)+\sum_j
p\_\lambda(\|\beta_j\|).\$\$ Here \\y_c\\ is the centered response,
\\Z\\ the standardized design and \\p\_\lambda\\ the fixed penalty
documented in `lambert_penalty`. The entry score is
\\\lambda\_{\max}=\max_j\|Z_j^T y_c\|/n\\.

Each path position uses both the previous verified solution and the zero
vector as starts. Among numerically accepted candidates the smaller
objective is chosen, with relative tolerance \\10^{-10}\\ and the warm
start favored in a tie. A failed position resets the next warm start to
zero. Diagnostics include the normalized KKT residual, objective,
iterations and both starts. When the two initial vectors are exactly
identical and the warm-start result passes the solver and independent
KKT checks, that result is reused for the zero-start record. Distinct
starts and failed warm starts are solved separately. Curvature on the
active coordinates is also recorded; a negative local curvature rejects
verification. At a curvature join this check is reported as unavailable
rather than as a positive-curvature certificate. These checks do not
certify a global minimum of a nonconvex objective.

The numerical controls preserve the reference study's defaults:
objective increase tolerance \\10^{-9}\\, scaling tolerance
\\10^{-12}\\, curvature join tolerance \\10^{-8}\\ and
negative-curvature tolerance \\10^{-8}\\. No automatic column removal,
imputation, response transformation, or relaxed convergence is applied.
Constant or nearly constant training columns and a zero entry score
produce an error. Unverified solutions are retained for diagnosis but
refused by the prediction and coefficient methods.

## Value

A list of class `lambert`. With \\p\\ predictors and \\L\\ penalty
levels, its components are:

- beta:

  Numeric \\p\times L\\ original-scale coefficient matrix; rows are
  predictors and columns are path positions, named `s1`, etc.

- intercept:

  Length-\\L\\ original-scale intercept vector.

- beta_standardized:

  Numeric \\p\times L\\ matrix in the internally standardized predictor
  coordinates.

- lambda, lambda_fraction:

  Length-\\L\\ penalty levels and fractions of the full-training entry
  score.

- center, scale:

  Length-\\p\\ predictor means and centered RMS scales.

- response_mean, lambda_max:

  Response mean and entry score.

- nobs, nvars, xnames:

  Training dimensions and supplied predictor names (or `NULL` if
  unnamed). Coefficient row names are generated if needed.

- diagnostics:

  One row per path position. Columns include `index`, `lambda`, `ok`,
  `status`, normalized `kkt`, `objective`, `iterations`,
  `selected_start`, `elapsed` (seconds), `max_increase`, `warning` and
  `error`. The `solver_lambda` and `solver_shape` columns record solver
  inputs. These are solver/KKT records; consult `verified` for the
  additional active-curvature assessment.

- attempts:

  The same diagnostic columns for both logical starts at each position,
  plus logical `executed` and `reused` columns. A reused record has
  `executed = FALSE`, `reused = TRUE`, zero `elapsed` and the iteration
  count of the verified warm-start solve it shares.

- curvature:

  Length-\\L\\ list with `status`, `min_eigenvalue` and `condition` for
  each active-set Hessian.

- verified, status:

  Length-\\L\\ logical acceptance indicators and character status values
  after the active-curvature check.

- control, shape, call:

  Numerical settings, fixed shape (one), and call.

- elapsed:

  Total solver path-fitting elapsed time in seconds; subsequent
  wrapper-level curvature checks are not included.

## Status codes

`verified_stationary` means the solver and KKT checks passed.
`budget_exhausted` means a start did not converge within its sweep
budget; `kkt_not_met` means the independent KKT check failed;
`objective_increased` means the objective monotonicity check failed; and
`solver_error` records a caught solver error. `negative_curvature`
rejects an otherwise accepted path solution. Curvature status may also
be `nonnegative_local`, `empty_active`, `join_unavailable`, or
`unavailable`. A join or empty active set is not a positive-curvature
certificate. Unavailable solutions may contain `NA` coefficients and
must not be used for prediction.

## See also

[`cv.lambert`](https://byuzbasi.github.io/lambertReg/reference/cv.lambert.md),
[`lambert_penalty`](https://byuzbasi.github.io/lambertReg/reference/lambert_penalty.md),
[`predict.lambert`](https://byuzbasi.github.io/lambertReg/reference/methods.md)

## Examples

``` r
set.seed(7)
x <- matrix(rnorm(60 * 6), 60, 6)
y <- 2 * x[, 1] - x[, 2] + rnorm(60)
fit <- lambert(x, y, lambda_fraction = c(1, 0.5, 0.2))
coef(fit, index = 2)
#>                     s2
#> (Intercept) 0.03569433
#> V1          1.44891258
#> V2          0.00000000
#> V3          0.00000000
#> V4          0.00000000
#> V5          0.00000000
#> V6          0.00000000
predict(fit, x[1:3, , drop = FALSE], index = 2)
#>              s2
#> [1,]  3.3497155
#> [2,] -1.6983232
#> [3,] -0.9702748
```
