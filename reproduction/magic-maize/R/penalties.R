# Scalar penalties for normalized quadratic loss (coordinate curvature = 1).
# Existing dimensionless plot interfaces are retained; shape defaults to e.
lambert_w_unit <- function(x) {
  cutoff <- exp(1)
  stopifnot(all(is.finite(x)), all(x >= 0), all(x <= cutoff))
  vapply(x, function(value) {
    if (value == 0) return(0)
    if (value == cutoff) return(1)
    uniroot(function(w) w * exp(w) - value, c(0, 1), tol = 1e-14)$root
  }, numeric(1))
}

penalty_lasso <- function(u) u
penalty_mcp <- function(u, a = exp(1)) {
  stopifnot(a > 1)
  v <- pmin(u, a)
  v - v^2 / (2 * a)
}
penalty_candidate <- function(u) {
  w <- lambert_w_unit(pmin(u, exp(1)))
  (expm1(2 * w) + exp(2 * w) * (2 * w - 2 * w^2)) / 4
}
penalty_scad <- function(u, a = exp(1)) {
  stopifnot(a > 2)
  ifelse(u <= 1, u, ifelse(u <= a,
    (-u^2 + 2 * a * u - 1) / (2 * (a - 1)), (a + 1) / 2))
}
shrink_lasso <- function(s) sign(s) * pmax(abs(s) - 1, 0)
shrink_mcp <- function(s, a = exp(1)) {
  stopifnot(a > 1)
  sign(s) * pmin(abs(s), a / (a - 1) * pmax(abs(s) - 1, 0))
}
shrink_candidate <- function(s) {
  value <- abs(s)
  result <- numeric(length(s))
  middle <- value > 1 & value < exp(1)
  result[middle] <- value[middle] * log(value[middle])
  result[value >= exp(1)] <- value[value >= exp(1)]
  sign(s) * result
}
shrink_scad <- function(s, a = exp(1)) {
  stopifnot(a > 2)
  value <- abs(s)
  sign(s) * ifelse(value <= 2, pmax(value - 1, 0),
    ifelse(value <= a, ((a - 1) * value - a) / (a - 2), value))
}

penalty_value <- function(beta, lambda, method, a = exp(1)) {
  stopifnot(lambda > 0)
  u <- abs(beta) / lambda
  lambda^2 * switch(method, lasso = penalty_lasso(u), mcp = penalty_mcp(u, a),
    scad = penalty_scad(u, a), candidate = penalty_candidate(u), stop("Unknown penalty"))
}
penalty_slope <- function(t, lambda, method, a = exp(1)) {
  stopifnot(all(t >= 0), lambda > 0)
  u <- t / lambda
  switch(method,
    lasso = rep(lambda, length(t)),
    mcp = lambda * pmax(1 - u / a, 0),
    scad = lambda * ifelse(u <= 1, 1, pmax((a - u) / (a - 1), 0)),
    candidate = {
      w <- lambert_w_unit(pmin(u, exp(1)))
      lambda * exp(w) * (1 - w)
    }, stop("Unknown penalty"))
}
penalty_curvature <- function(t, lambda, method, a = exp(1)) {
  u <- t / lambda
  switch(method,
    lasso = rep(0, length(t)),
    mcp = ifelse(u < a, -1 / a, 0),
    scad = ifelse(u > 1 & u < a, -1 / (a - 1), 0),
    candidate = {
      w <- lambert_w_unit(pmin(u, exp(1)))
      ifelse(u < exp(1), -w / (1 + w), 0)
    }, stop("Unknown penalty"))
}
threshold_value <- function(z, lambda, method, a = exp(1)) {
  s <- z / lambda
  lambda * switch(method, lasso = shrink_lasso(s), mcp = shrink_mcp(s, a),
    scad = shrink_scad(s, a), candidate = shrink_candidate(s), stop("Unknown penalty"))
}
