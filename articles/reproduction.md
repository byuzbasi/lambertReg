# Data and reproducibility

## Know your inputs

[`lambert()`](https://byuzbasi.github.io/lambertReg/reference/lambert.md)
and
[`cv.lambert()`](https://byuzbasi.github.io/lambertReg/reference/cv.lambert.md)
take a finite dense numeric matrix `x` and a numeric response vector
`y`. Rows represent observations; columns represent predictors.
Predictor names, when present, must be unique. Prediction data must use
the same names and column order. Missing values are not imputed, and
columns must vary in every training sample.

The package centers and RMS-scales predictors separately in each
training sample, centers the response, and leaves the intercept
unpenalized. The returned coefficients and predictions are in the
original units. The
[introduction](https://byuzbasi.github.io/lambertReg/articles/introduction.html#penalty-scale-and-preprocessing)
gives the precise transformation and penalty scale.

## Demonstration data dictionary

The example generator uses `set.seed(42)` and generates 100 independent
observations. Its response is `2*x1 - x2 + 0.5*x3 + noise`, with
independent standard-normal noise. The eight predictors are independent
standard-normal variables. The demonstration has no physical units.

| File | Fields | Meaning |
|:---|:---|:---|
| `predictors.csv` | `x1`–`x8` | Eight synthetic predictors; 100 rows in observation order. |
| `observations.csv` | `id`, `y`, `sample`, `fold` | Response, training/test split, and CV fold; fold is missing for held-out rows. |
| `configuration.rds` | `seed`, `train`, `foldid`, `lambda_fraction`, `controls`, `package_version` | Frozen example settings and fitting controls. |
| `cv_path.csv` | `lambda`, `fraction`, `cv_mse`, `selected` | Full-training penalty scale, mean fold MSE, and selected-variable count. |
| `coefficients.csv` | `predictor`, `s1`–`s60` | Original-scale coefficient paths, excluding the intercept. |
| `summary.csv` | `lambda`, `selected`, `cv_mse`, `test_rmse`, `test_mae` | Selected model and held-out errors for this demonstration. |
| `fit.rds` | A `cv.lambert` object | Stored paths, CV losses, preprocessing, and numerical diagnostics. |
| `versions.txt`, `manifest.csv` | Versions; `file`, `bytes`, `sha256` | Software record and output integrity checks. |

The first 80 observations form the training sample and the final 20 form
a held-out test sample. Training rows receive an explicit balanced
five-fold partition. The test rows do not affect scaling, tuning, or
selection. No observations are removed and no missing values are
generated in predictors or responses.

## Reproduce the gallery

The [example
archive](https://byuzbasi.github.io/lambertReg/examples/lambert-example-v1.zip)
contains the synthetic inputs, stored model, plots, configuration, and
checksum manifest listed above. The
[generator](https://byuzbasi.github.io/lambertReg/examples/showcase.R)
recreates these files from the frozen seed.

Install lambertReg and make the documentation dependency `digest`
available. Run the following in R, using a new output directory:

``` r

source(system.file("examples", "showcase.R", package = "lambertReg"))
example <- lambert_showcase("lambert-example-v1")
example$summary
```

The generator refuses an existing directory. It saves the input data,
folds, controls, fit, CV curve, coefficient paths, two PNG plots,
software versions, and a SHA-256 manifest. It verifies the selected fit
and the completeness of the demonstration’s CV and full-training paths.
This is one small example; it does not launch a Monte Carlo study.

To verify saved files:

``` r

out <- "lambert-example-v1"
manifest <- read.csv(file.path(out, "manifest.csv"))
paths <- file.path(out, manifest$file)
stopifnot(all(file.info(paths)$size == manifest$bytes))
hashes <- vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), "")
stopifnot(identical(unname(hashes), manifest$sha256))
```

## Reproduce an independent analysis

Record the input source and any preprocessing outside the package, the
training/test split, `foldid`, penalty grid, fitting controls, and
software versions. Preserve raw data separately. Use an independent test
sample or outer CV for evaluation; the tuning minimum is not an external
prediction error. This guide reproduces the package demonstration. The
manuscript’s simulation and QSAR study use their separately recorded
designs and inputs.

## Cite the software

``` r

citation("lambertReg")
```

    ## To cite package 'lambertReg' in publications use:
    ## 
    ##   Yuzbasi B (2026). lambertReg: Sparse Gaussian Regression with the
    ##   Lambert Penalty. R package version 0.1.10
    ## 
    ##   Yuzbasi B (2026). The Lambert Penalty: Logarithmic Shrinkage for
    ##   Sparse Regression. arXiv preprint arXiv:2610.09627 [stat.ME].
    ##   https://arxiv.org/abs/2610.09627
    ## 
    ## To see these entries in BibTeX format, use 'print(<citation>,
    ## bibtex=TRUE)', 'toBibtex(.)', or set
    ## 'options(citation.bibtex.max=999)'.

The source repository provides `CITATION.cff` with the associated paper
as the preferred citation. The paper is [*The Lambert Penalty:
Logarithmic Shrinkage for Sparse
Regression*](https://arxiv.org/abs/2610.09627), by Bahadir Yuzbasi
(2026), arXiv preprint **arXiv:2610.09627** \[stat.ME\], first posted on
7 October 2026. `citation("lambertReg")` returns both the software and
paper citations.
