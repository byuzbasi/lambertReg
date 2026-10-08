\name{lambertReg-package}
\alias{lambertReg-package}
\alias{lambertReg}
\docType{package}
\title{Sparse Gaussian Regression with the Lambert Penalty}
\description{An R/C++ implementation of fixed-shape Lambert regression with
regularization paths, five-fold cross-validation by default, and explicit
numerical diagnostics.}
\details{The package implements the main estimator in the arXiv preprint
\emph{The Lambert Penalty: Logarithmic Shrinkage for Sparse Regression}
by Bahadir Yuzbasi (2026), arXiv:2610.09627 [stat.ME], available at
\url{https://arxiv.org/abs/2610.09627}. It does not implement the extended shape family,
Lambert-Min post-selection procedures, or comparator estimators.
Start with \code{cv.lambert} to select a penalty, then use \code{coef},
\code{predict} and \code{plot}. Use \code{lambert} for a path with specified
penalties. All predictions and reported coefficients use the original input
scale. Training-only predictor centering and RMS scaling are automatic.
The worked guide is available through
\code{vignette("introduction", package = "lambertReg")}.
Source code and issue reporting are available at
\url{https://github.com/byuzbasi/lambertReg}.}
\author{Bahadir Yuzbasi}
\seealso{\code{\link{lambert}}, \code{\link{cv.lambert}}}
