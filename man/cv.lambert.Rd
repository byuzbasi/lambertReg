\name{cv.lambert}
\alias{cv.lambert}
\title{Training-Only Cross-Validation for Lambert Regression}
\description{Choose a relative penalty-path position by the arithmetic mean
of validation-fold squared errors, with fixed Lambert shape and fold-specific
preprocessing.}
\usage{
cv.lambert(x, y, nfolds = 5L, foldid = NULL,
           lambda_fraction = 10^seq(0, -3, length.out = 60L),
           max_sweeps = 10000L, kkt_tol = 1e-7)
}
\arguments{
\item{x, y}{Training predictors and response, as in \code{\link{lambert}}.}
\item{nfolds}{Number of folds, default five. Each training subset must contain
at least three observations. Inferred from \code{foldid} when that argument
is supplied and \code{nfolds} is omitted.}
\item{foldid}{Optional integer fold label for every row, with consecutive labels
starting at one and no empty folds. Store and reuse this vector for matched
comparisons. Without it, balanced labels are randomly permuted using the
caller's R random-number state. Use \code{set.seed()} for reproducibility.}
\item{lambda_fraction}{Strictly decreasing fractions in \eqn{(0,1]}; default
60 logarithmically spaced positions from 1 to 0.001.}
\item{max_sweeps, kkt_tol}{Solver controls described in \code{\link{lambert}}.}
}
\details{
Each training fold estimates its own predictor centers, RMS scales, response
mean and entry score. Candidate \eqn{k} in fold \eqn{f} uses
\eqn{\lambda_{k,f}=\lambda_{\max,f}\tau_k}, where \eqn{\tau_k} is the supplied
fraction. The candidate is aligned by fraction, not by a common absolute
lambda. The validation observations are never used to estimate these values.

For validation sets \eqn{I_f} and training-fold predictions \eqn{\hat y_i},
the score is
\deqn{\mathrm{CV}(k)=\frac{1}{K}\sum_{f=1}^K
       \frac{1}{|I_f|}\sum_{i\in I_f}(y_i-\hat y_i)^2.}
This is an unweighted mean of the fold MSEs, including when fold sizes differ.
Only candidates with finite loss in every fold are eligible. Ties within
\eqn{10^{-10}(1+\mathrm{CV}_{\min})} favor the larger penalty fraction.
\code{lambda.min} is the selected fraction multiplied by the full training
entry score. A full-training path is fitted at every supplied fraction. Later positions
provide variable counts for plotting and do not affect CV selection. The
selected full and fold fits receive the reference active-curvature check.
A rejected selected fit is reported explicitly; another candidate is not
silently substituted.

An incomplete search is flagged by \code{partial_search}; it does not imply
that the selected fit failed. Solver convergence and curvature status must
be inspected separately. No one-standard-error rule or optional early
pruning is implemented in this initial package. CV error is a tuning score,
not an unbiased external performance estimate; use independent test data or
outer cross-validation to evaluate prediction. Changing the number of folds
or numerical controls is supported but departs from the study defaults.
}
\value{A list of class \code{cv.lambert}. Write \eqn{L} for the
number of fractions, \eqn{K} for the number of folds, and \eqn{n} for the
number of input rows. Components are:
\describe{
\item{cvm}{Length-\eqn{L} arithmetic means of fold MSEs. An ineligible
candidate has value \code{Inf}.}
\item{fold_loss}{Numeric \eqn{L\times K} validation-MSE matrix, with
\code{NA} where a fold candidate is unavailable.}
\item{foldid, nfolds}{Length-\eqn{n} fold labels and fold count.}
\item{index}{Selected path index, or \code{NA_integer_} when none is eligible.}
\item{lambda_fraction, lambda}{Length-\eqn{L} fractions and corresponding
full-training penalty levels. These absolute levels need not equal fold levels.}
\item{fraction.min, lambda.min}{Selected fraction and full-training penalty;
\code{NA_real_} if no candidate is eligible.}
\item{lambert.fit}{Full-training \code{lambert} path, including when no CV
candidate is eligible. Its fields are documented in \code{\link{lambert}}.}
\item{nzero, support_tol}{Length-\eqn{L} selected-variable counts and their
threshold. Counts exclude the intercept and use original-scale coefficient
magnitudes greater than \eqn{10^{-8}}. Unverified counts are \code{NA}.}
\item{full_path_complete}{Whether every full-training position is verified;
this is separate from fold-search completeness and selected-fit validity.}
\item{ok, status}{Logical selected-fit acceptance and character status.
\code{no_valid_cv} means no candidate had finite loss in all folds.
Other status values follow \code{\link{lambert}}. A negative-curvature
selected fold can set \code{status = "negative_curvature"}.}
\item{partial_search}{Whether any candidate has a nonfinite CV score.
A partial search can still have \code{ok = TRUE}.}
\item{cv_diagnostics}{Fold-stacked solver diagnostics, with additional
\code{fold} and validation \code{loss} columns.}
\item{attempts}{Fold-stacked records for both starts, with a \code{fold} column.
Full-training attempts are in \code{lambert.fit$attempts}.}
\item{selected_curvature}{Length-\eqn{K} list of selected-fold Hessian
assessments; \code{NULL} if no CV candidate is eligible. Full-training
assessments are in \code{lambert.fit$curvature}.}
\item{preprocessing}{Length-\eqn{K+1} list: each training fold, then the full
training sample. Each entry contains \code{mx}, \code{my}, \code{sx} and
\code{lambda0}, the predictor means, response mean, RMS scales and entry score.}
\item{control, shape, call}{Numerical settings, fixed shape (one), and call.}
}
Completing the full-training path can cost more than fitting only its selected
prefix. Plotting reads stored values and never fits models.
}
\seealso{\code{\link{lambert}}, \code{\link{predict.lambert}}}
\examples{
set.seed(12)
x <- matrix(rnorm(60 * 6), 60, 6)
y <- 2 * x[, 1] - x[, 2] + rnorm(60)
foldid <- sample(rep(1:5, length.out = 60))
fit <- cv.lambert(x, y, foldid = foldid,
                  lambda_fraction = c(1, 0.5, 0.2))
print(fit)
coef(fit)
predict(fit, x[1:3, , drop = FALSE])
}
