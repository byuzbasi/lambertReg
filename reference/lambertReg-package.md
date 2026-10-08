# Sparse Gaussian Regression with the Lambert Penalty

An R/C++ implementation of fixed-shape Lambert regression with
regularization paths, five-fold cross-validation by default, and
explicit numerical diagnostics.

## Details

The package implements the main estimator in the arXiv preprint *The
Lambert Penalty: Logarithmic Shrinkage for Sparse Regression* by Bahadir
Yuzbasi (2026), arXiv:2610.09627 \[stat.ME\], available at
<https://arxiv.org/abs/2610.09627>. It does not implement the extended
shape family, Lambert-Min post-selection procedures, or comparator
estimators. Start with `cv.lambert` to select a penalty, then use
`coef`, `predict` and `plot`. Use `lambert` for a path with specified
penalties. All predictions and reported coefficients use the original
input scale. Training-only predictor centering and RMS scaling are
automatic. The worked guide is available through
[`vignette("introduction", package = "lambertReg")`](https://byuzbasi.github.io/lambertReg/articles/introduction.md).
Source code and issue reporting are available at
<https://github.com/byuzbasi/lambertReg>.

## Author

Bahadir Yuzbasi

## See also

[`lambert`](https://byuzbasi.github.io/lambertReg/reference/lambert.md),
[`cv.lambert`](https://byuzbasi.github.io/lambertReg/reference/cv.lambert.md)
