# Scalar Lambert family for loss 0.5 * (z - b)^2; see protocol/lambert_family_v1.md.
# shape >= 1 is explicit: shape = 1 recovers the original Candidate B.
# t denotes a nonnegative coefficient magnitude, not a signed coefficient.
lambert_family_constants <- function(lambda, shape) {
  if (!is.numeric(lambda) || length(lambda) != 1L ||
      !is.finite(lambda) || lambda <= 0) stop("lambda must be finite and positive.")
  if (!is.numeric(shape) || length(shape) != 1L ||
      !is.finite(shape) || shape < 1) stop("shape must be finite and at least 1.")
  cutoff <- lambda * exp(shape)
  plateau <- lambda^2 * expm1(2 * shape) / (4 * shape)
  if (!is.finite(cutoff) || !is.finite(plateau) || plateau <= 0) {
    stop("lambda and shape exceed the representable numerical range.")
  }
  list(cutoff = cutoff, plateau = plateau)
}

lambert_family_arguments <- function(x, nonnegative = TRUE) {
  if (!is.numeric(x) || any(!is.finite(x)) ||
      (nonnegative && any(x < 0))) {
    stop(if (nonnegative) "t must contain finite nonnegative magnitudes."
         else "z must contain finite numeric scores.")
  }
  invisible(TRUE)
}

# Returns min(W_0(shape * t / lambda), shape), avoiding unconstrained branches.
lambert_family_w <- function(t, lambda, shape) {
  constants <- lambert_family_constants(lambda, shape)
  lambert_family_arguments(t)
  result <- rep(shape, length(t))
  result[t == 0] <- 0
  interior <- t > 0 & t < constants$cutoff
  target <- shape * (t[interior] / lambda)
  result[interior] <- vapply(target, function(value) {
    uniroot(function(w) w * exp(w) - value, c(0, shape),
            tol = 1e-14, check.conv = TRUE)$root
  }, numeric(1))
  result
}

lambert_penalty <- function(t, lambda, shape) {
  constants <- lambert_family_constants(lambda, shape)
  w <- lambert_family_w(t, lambda, shape)
  # expm1 avoids subtraction of two nearly equal constants near t = 0.
  result <- lambda^2 * (expm1(2 * w) / (4 * shape) +
    exp(2 * w) * w * (shape - w) / (2 * shape^2))
  result[t >= constants$cutoff] <- constants$plateau
  result
}

lambert_slope <- function(t, lambda, shape) {
  constants <- lambert_family_constants(lambda, shape)
  w <- lambert_family_w(t, lambda, shape)
  result <- lambda * exp(w) * (1 - w / shape)
  result[t >= constants$cutoff] <- 0
  result
}

# At zero return p''(0+). At the exact identity cutoff p'' does not exist.
# This is the curvature in magnitude t; the signed penalty has a cusp at zero.
lambert_curvature <- function(t, lambda, shape) {
  constants <- lambert_family_constants(lambda, shape)
  w <- lambert_family_w(t, lambda, shape)
  result <- shape / (1 + w) - 1
  result[t > constants$cutoff] <- 0
  result[t == constants$cutoff] <- NA_real_
  result
}

lambert_threshold <- function(z, lambda, shape) {
  constants <- lambert_family_constants(lambda, shape)
  lambert_family_arguments(z, nonnegative = FALSE)
  magnitude <- abs(z)
  result <- numeric(length(z))
  interior <- magnitude > lambda & magnitude < constants$cutoff
  result[interior] <- magnitude[interior] / shape *
    log1p((magnitude[interior] - lambda) / lambda)
  tail <- magnitude >= constants$cutoff
  result[tail] <- magnitude[tail]
  sign(z) * result
}
