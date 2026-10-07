# lambertReg

Sparse Gaussian regression with the **Lambert penalty**. The package provides
regularization paths, cross-validation, prediction, and plots of CV error and
coefficient paths. The shape is fixed at c = 1; cross-validation selects the
penalty strength.

## Installation

R (>= 4.1.0), Rcpp, and an R-compatible C++ compiler are required.
Install from a local clone:

```sh
git clone https://github.com/byuzbasi/lambertReg.git
R CMD INSTALL lambertReg
```

To include the worked vignette, install the suggested packages knitr and
rmarkdown and make Pandoc available, then build before installing:

```sh
R CMD build lambertReg
R CMD INSTALL lambertReg_0.1.10.tar.gz
```

## Quick start

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
```

The default is five-fold CV over 60 penalty fractions. Predictors are centered
and scaled within each training fold; coefficients and predictions are
returned on the original scale. `fit$lambda.min` is the CV-selected penalty
for the full training data. Use held-out observations to assess prediction
performance.

## Plots

```r
plot(fit)                                  # CV error, error bars and variable counts
plot(fit, style = "paper")                 # compact publication style
plot(fit, type = "coefficients")           # coefficient paths
plot(fit, type = "coefficients", labels = TRUE)
```

The CV plot marks the selected penalty and displays the number of selected
variables above the curve. Error bars show one fold-based standard error.
Plots use stored results and do not refit models.

## Documentation

The [introductory vignette](vignettes/introduction.Rmd) works through model
fitting, tuning, prediction, plots, and numerical diagnostics.

```r
vignette("introduction", package = "lambertReg")
?lambert
?cv.lambert
?lambert_penalty
```

Inputs are finite dense numeric matrices with columns that vary in each
training sample. Inspect `fit$ok`, `fit$partial_search`, and
`fit$cv_diagnostics` when checking a fitted model. Detailed controls and
status definitions are documented in the help pages.

## Citation and support

Use `citation("lambertReg")` for the package citation. Report issues through
[GitHub Issues](https://github.com/byuzbasi/lambertReg/issues).

GPL-3. Maintainer: Bahadir Yuzbasi <b.yzb@hotmail.com>.
