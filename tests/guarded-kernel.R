# Bounded regressions against the previous compiled scalar solver.
library(lambertReg)
ns <- asNamespace("lambertReg")
old <- new.env(parent = ns)
for (name in c("lambert", "cv.lambert", "perf_path")) {
  fn <- get(name, ns)
  environment(fn) <- old
  old[[name]] <- fn
}
old$lambert_cd_guarded_arma_v2 <- function(X, y, init, lambda, shape,
    max_sweeps, kkt_tol, pulse, profile = FALSE, cache = NULL)
  get("lambert_cd_cpp", ns)(X, y, init, lambda, shape, max_sweeps, kkt_tol, pulse)
strip <- function(x) {
  if (!is.list(x)) return(x)
  if (!is.null(names(x))) x <- x[!names(x) %in% c("elapsed", "call")]
  for (i in seq_along(x)) x[i] <- list(strip(x[[i]]))
  x
}
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
grid <- c(1, .5, .2, .06, .01)
for (case in 1:4) {
  set.seed(2026102000L + case)
  n <- if (case == 2L) 32L else 60L
  p <- if (case == 2L) 80L else 12L
  x <- matrix(rnorm(n*p), n, p)
  if (case >= 3L) {
    rho <- if (case == 3L) .95 else .99999
    x[, 2] <- rho*x[, 1] + sqrt(1-rho^2)*x[, 2]
  }
  y <- 1.5*x[, 1] - .8*x[, 3] + rnorm(n, sd = .5)
  foldid <- rep(1:5, length.out = n)
  saved_rng <- .Random.seed
  for (budget in c(10000L, 1L)) {
    a <- suppressWarnings(old$lambert(x, y, lambda_fraction = grid, max_sweeps = budget))
    b <- suppressWarnings(lambert(x, y, lambda_fraction = grid, max_sweeps = budget))
    stopifnot(identical(strip(a), strip(b)),
      identical(abs(a$beta) > 1e-8, abs(b$beta) > 1e-8))
    k <- which(a$verified)
    if (length(k)) stopifnot(identical(coef(a, index = k), coef(b, index = k)),
      identical(predict(a, x, index = k), predict(b, x, index = k)))
  }
  a <- suppressWarnings(old$cv.lambert(x, y, foldid = foldid, lambda_fraction = grid))
  b <- suppressWarnings(cv.lambert(x, y, foldid = foldid, lambda_fraction = grid))
  stopifnot(identical(strip(a), strip(b)), identical(saved_rng, .Random.seed))
  if (a$ok) stopifnot(identical(coef(a), coef(b)), identical(predict(a, x), predict(b, x)))
  # Cached and standalone calls, including a context surviving collection.
  cfg <- get(".lr_config", ns)(10000L, 1e-7)
  pp <- get(".lr_prepare", ns)(x, y, cfg)
  context <- get("lambert_context_guarded_arma_v2", ns)(pp$Z, pp$y)
  invisible(gc())
  solve <- function(cache = NULL, profile = FALSE)
    get("lambert_cd_guarded_arma_v2", ns)(pp$Z, pp$y, numeric(p),
      pp$lambda0*.2, 1, 10000L, 5e-8, function(...) NULL, profile, cache)
  fields <- c("beta", "iterations", "converged", "objective", "kkt",
    "max_increase", "increase_failed")
  f <- solve(context, TRUE)
  stopifnot(identical(f[fields], solve()[fields]),
    f$full_scans == min(f$iterations, 3L) + floor(f$iterations/25) +
      as.integer(f$iterations == 10000L && f$iterations > 3 && f$iterations %% 25 != 0),
    fails(get("lambert_cd_guarded_arma_v2", ns)(pp$Z + 0, pp$y, numeric(p),
      .1, 1, 10L, 5e-8, function(...) NULL, FALSE, context)))
}
# Scalar paths and explicit lambda values retain the same public contract.
x <- matrix(seq(-2, 2, length.out = 30L), ncol = 1L)
y <- 2*x[, 1] + sin(seq_len(30L))
a <- old$lambert(x, y, lambda = c(3, 1, .2))
b <- lambert(x, y, lambda = c(3, 1, .2))
stopifnot(identical(strip(a), strip(b)))
# A public fit must contain no nonportable context or process-global cache.
has_pointer <- function(x) typeof(x) == "externalptr" ||
  (is.list(x) && any(vapply(x, has_pointer, logical(1))))
stopifnot(!has_pointer(b), identical(b, unserialize(serialize(b, NULL))))
stopifnot(fails(get("lambert_context_guarded_arma_v2", ns)(matrix(NA_real_, 3, 2), 1:3)))
cat("PASS: ordered-kernel paths, CV, failures, predictions, context and serialization regressions.\n")
