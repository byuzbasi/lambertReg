\name{lambert}
\alias{lambert}
\title{Gaussian Regression Paths with the Fixed Lambert Penalty}
\description{Fit a regularization path with the bounded Lambert penalty of shape
\eqn{c=1}. Exact cyclic scalar updates are computed by a compiled
solver that preserves the reference fitting rule. No shape selection or post-selection refit is performed.}
\usage{
lambert(x, y, lambda = NULL,
        lambda_fraction = 10^seq(0, -3, length.out = 60L),
        max_sweeps = 10000L, kkt_tol = 1e-7)
}
\arguments{
\item{x}{Finite numeric predictor matrix, with at least three rows and one
column. Columns must vary in every training sample. Optional column names must
be nonempty and unique. Sparse matrices, missing values and automatic factor
encoding are not supported.}
\item{y}{Finite numeric response vector. The response is centered, not scaled
to unit variance. A positive entry score is required.}
\item{lambda}{Optional strictly decreasing vector of positive penalty levels
on the standardized-predictor scale. Cannot be supplied together with
\code{lambda_fraction}. Zero is not supported.}
\item{lambda_fraction}{Strictly decreasing fractions in \eqn{(0,1]}, multiplied
by the training entry score \eqn{\lambda_{\max}}. Defaults to 60 positions
from 1 to 0.001.}
\item{max_sweeps}{Positive integer sweep budget for each start and penalty level.}
\item{kkt_tol}{Positive normalized KKT acceptance tolerance. The internal solver
uses half this tolerance, followed by an independent residual check.}
}
\details{
For \eqn{n} training observations, each predictor is centered and divided by
its root-mean-square centered value so that the standardized column has squared
norm \eqn{n}. The intercept is unpenalized. All centering and scaling statistics
are estimated from the supplied training observations. The objective is
\deqn{\|y_c-Z\beta\|_2^2/(2n)+\sum_j p_\lambda(|\beta_j|).}
Here \eqn{y_c} is the centered response, \eqn{Z} the standardized design and
\eqn{p_\lambda} the fixed penalty documented in \code{lambert_penalty}.
The entry score is \eqn{\lambda_{\max}=\max_j|Z_j^T y_c|/n}.

Each path position uses both the previous solver/KKT-accepted solution and the
zero vector as starts. Among solver/KKT-accepted candidates the smaller objective
is chosen, with relative tolerance \eqn{10^{-10}} and the warm start favored
in a tie. A failed position resets the next warm start to zero. Diagnostics
include the normalized KKT residual, objective, iterations and both starts.
When the two initial vectors are exactly identical and the warm-start result
passes the solver and independent KKT checks, that result is reused for the
zero-start record. Distinct starts and failed warm starts are solved separately.
The compiled solver caches immutable training-data checks within each path.
A conservative score bound can omit an update of a zero coefficient when that
update is guaranteed to remain zero; uncertain coordinates are updated in the
original order. Full KKT and objective checks retain the reference cadence:
the first three sweeps, every 25th sweep, and the final budgeted sweep.
Every full KKT check includes every predictor. RcppArmadillo supplies working
storage; coordinate inner products retain the reference scalar summation order.
Curvature on the active coordinates is then recorded for the retained path.
Verification requires
an assessed nonnegative active-set Hessian within tolerance, or an empty
active set. Negative curvature and unavailable curvature both reject
verification. At a curvature join, \code{status = "join_unavailable"} and
\code{verified = FALSE}; coefficients remain available in the returned object
for diagnosis. These are local numerical checks, not a global-minimum certificate.

The numerical controls preserve the reference study's defaults: objective
increase tolerance \eqn{10^{-9}}, scaling tolerance \eqn{10^{-12}}, curvature
join tolerance \eqn{10^{-8}} and negative-curvature tolerance \eqn{10^{-8}}.
No automatic column removal, imputation, response transformation, or relaxed
convergence is applied. Constant or nearly constant training columns and a
zero entry score produce an error. Unverified solutions are retained for
diagnosis but refused by the prediction and coefficient methods.
}
\value{A list of class \code{lambert}. With \eqn{p} predictors and
\eqn{L} penalty levels, its components are:
\describe{
\item{beta}{Numeric \eqn{p\times L} original-scale coefficient matrix;
rows are predictors and columns are path positions, named \code{s1}, etc.}
\item{intercept}{Length-\eqn{L} original-scale intercept vector.}
\item{beta_standardized}{Numeric \eqn{p\times L} matrix in the internally
standardized predictor coordinates.}
\item{lambda, lambda_fraction}{Length-\eqn{L} penalty levels and fractions
of the full-training entry score.}
\item{center, scale}{Length-\eqn{p} predictor means and centered RMS scales.}
\item{response_mean, lambda_max}{Response mean and entry score.}
\item{nobs, nvars, xnames}{Training dimensions and supplied predictor names
(or \code{NULL} if unnamed). Coefficient row names are generated if needed.}
\item{diagnostics}{One row per path position. Columns include \code{index},
\code{lambda}, \code{ok}, \code{status}, normalized \code{kkt},
\code{objective}, \code{iterations}, \code{selected_start}, \code{elapsed}
(seconds), \code{max_increase}, \code{warning} and \code{error}.
The \code{solver_lambda} and \code{solver_shape} columns record solver inputs.
These are solver/KKT records; consult \code{verified} for the additional
active-curvature assessment.}
\item{attempts}{The same diagnostic columns for both logical starts at each
position, plus logical \code{executed} and \code{reused} columns. A reused
record has \code{executed = FALSE}, \code{reused = TRUE}, zero \code{elapsed}
and the iteration count of the verified warm-start solve it shares.}
\item{curvature}{Length-\eqn{L} list with \code{status},
\code{min_eigenvalue} and \code{condition} for each active-set Hessian.}
\item{verified, status}{Length-\eqn{L} logical acceptance indicators and
character status values after the active-curvature check.}
\item{control, shape, call}{Numerical settings, fixed shape (one), and call.}
\item{elapsed}{Total solver path-fitting elapsed time in seconds;
subsequent wrapper-level curvature checks are not included.}
}}
\section{Status codes}{
In \code{diagnostics}, \code{verified_stationary} means the solver and KKT
checks passed. The top-level \code{status} and \code{verified} additionally
require the active-curvature acceptance rule above.
\code{budget_exhausted} means a start did not converge within its sweep budget;
\code{kkt_not_met} means the independent KKT check failed;
\code{objective_increased} means the objective monotonicity check failed;
and \code{solver_error} records a caught solver error.
\code{negative_curvature} rejects an otherwise accepted path solution.
Curvature status may also be \code{nonnegative_local}, \code{empty_active},
\code{join_unavailable}, or \code{unavailable}. A join yields top-level
\code{join_unavailable}; other missing or unknown assessments yield
\code{curvature_unavailable}. Both set \code{verified = FALSE}. An empty
active set passes the acceptance rule without asserting strictly positive
curvature. Unverified solutions may have finite coefficients or \code{NA}
values and are refused by the prediction and coefficient methods.
}
\seealso{\code{\link{cv.lambert}}, \code{\link{lambert_penalty}},
\code{\link{predict.lambert}}}
\examples{
set.seed(7)
x <- matrix(rnorm(60 * 6), 60, 6)
y <- 2 * x[, 1] - x[, 2] + rnorm(60)
fit <- lambert(x, y, lambda_fraction = c(1, 0.5, 0.2))
coef(fit, index = 2)
predict(fit, x[1:3, , drop = FALSE], index = 2)
}
