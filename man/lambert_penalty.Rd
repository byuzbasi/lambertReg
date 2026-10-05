\name{lambert_penalty}
\alias{lambert_penalty}
\alias{lambert_threshold}
\title{The Fixed Lambert Penalty and Scalar Update}
\description{Evaluate the bounded Lambert penalty or its exact scalar
least-squares update with fixed shape \eqn{c=1}.}
\usage{
lambert_penalty(t, lambda)
lambert_threshold(z, lambda)
}
\arguments{
\item{t}{Finite, nonnegative coefficient magnitudes.}
\item{z}{Finite signed scalar scores.}
\item{lambda}{A finite positive penalty level.}
}
\details{
Write \eqn{W_0} for the principal real Lambert function, characterized by
\eqn{W_0(u)\exp(W_0(u))=u} for \eqn{u\geq0}. For a magnitude \eqn{t}, let
\eqn{w=W_0(t/\lambda)}. The penalty is
\deqn{p_\lambda(t)=\frac{\lambda^2}{4}
 [e^{2w}(1+2w-2w^2)-1],\quad 0\leq t\leq e\lambda.}
For \eqn{t>e\lambda}, it is constant at \eqn{\lambda^2(e^2-1)/4}.
The expression is evaluated with cancellation-resistant arithmetic near zero.

The scalar update minimizes \eqn{(z-b)^2/2+p_\lambda(|b|)} over \eqn{b\in R}.
It is zero when \eqn{|z|\leq\lambda}, equals
\eqn{z\log(|z|/\lambda)} for \eqn{\lambda<|z|<e\lambda}, and equals
\eqn{z} when \eqn{|z|\geq e\lambda}. This normalization requires unit
coordinate curvature, obtained by training-only RMS scaling in \code{lambert}.
}
\value{An unnamed numeric vector of the same length as \code{t} or \code{z},
containing penalty values or scalar updates, respectively. A scalar
\code{lambda} is shared by all elements; elementwise penalty-level vectors
are not supported. Zero-length input vectors return zero-length vectors.}
\seealso{\code{\link{lambert}}}
\examples{
lambert_penalty(c(0, 1, exp(1), 4), lambda = 1)
lambert_threshold(c(-4, -1, 0, 1, 2, 4), lambda = 1)
}
