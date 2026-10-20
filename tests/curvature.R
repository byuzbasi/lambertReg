# Small deterministic regression checks; no production study.
library(lambertReg)
ns <- asNamespace("lambertReg")
gate <- get("perf_curvature_status", ns)
record <- function(s) list(status = s)
stopifnot(
  gate(list(record("nonnegative_local"))) == "verified_stationary",
  gate(list(record("empty_active"))) == "verified_stationary",
  gate(list(record("nonnegative_local"), record("empty_active"))) == "verified_stationary",
  gate(list(record("join_unavailable"))) == "join_unavailable",
  gate(list(record("nonnegative_local"), record("join_unavailable"))) == "join_unavailable",
  gate(list(record("negative_curvature"), record("join_unavailable"))) == "negative_curvature",
  gate(list(record("unavailable"))) == "curvature_unavailable",
  gate(list(record("new_unknown_status"))) == "curvature_unavailable",
  gate(list(record(NA_character_))) == "curvature_unavailable",
  gate(list(NULL)) == "curvature_unavailable",
  gate(list()) == "curvature_unavailable")

# A coordinate-stationary join can have a joint descent direction.
n <- 3L; rho <- .9; lambda <- 1
q1 <- c(1, -1, 0)/sqrt(2); q2 <- c(1, 1, -2)/sqrt(6)
x <- sqrt(n) * cbind(q1, rho*q1 + sqrt(1-rho^2)*q2)
beta <- c(exp(1), exp(1)); y <- drop(x %*% beta)
cfg <- get(".lr_config", ns)(10000L, 1e-7, lambda = lambda)
prep <- list(Z = x, y = y, mx = c(0, 0), my = 0, sx = c(1, 1),
  score = drop(crossprod(x, y)/n), lambda0 = max(abs(crossprod(x, y)/n)))
setting <- list(family = "lambert", shape = 1)
h <- get("perf_hessian", ns)(beta, prep, lambda, setting, cfg)
diagnosis <- get("perf_diagnose", ns)(beta, prep, lambda, setting, cfg)
objective <- function(b) mean((y - drop(x %*% b))^2)/2 +
  sum(lambert_penalty(abs(b), lambda))
stopifnot(diagnosis$ok, diagnosis$kkt < 1e-14,
  h$status == "join_unavailable",
  objective(beta + 1e-3*c(-1, .9)) < objective(beta) - 1e-8)
path <- list(beta = matrix(beta, 2L), elapsed = 0, attempts = data.frame(),
  diagnostics = data.frame(lambda = lambda, ok = TRUE, status = "verified_stationary"))
fit <- get(".lr_object", ns)(path, prep, cfg, NULL, quote(lambert(x, y)))
stopifnot(!fit$verified, fit$status == "join_unavailable",
  identical(as.numeric(fit$beta), beta), fit$diagnostics$ok)
fails <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
stopifnot(fails(coef(fit, index = 1L)), fails(predict(fit, newx = x, index = 1L)))

# Even a legitimate scalar minimum at the join remains unassessed, not
# mislabeled as negative curvature. The public interfaces enforce that status.
x <- matrix(c(-5, -3, -1, 1, 3, 5), ncol = 1L)
y <- 2 * x[, 1]
fraction <- exp(-1)
fit <- suppressWarnings(lambert(x, y, lambda_fraction = fraction))
stopifnot(all(fit$diagnostics$ok), !fit$verified,
  fit$status == "join_unavailable", all(is.finite(fit$beta)))
cv <- suppressWarnings(cv.lambert(x, y, foldid = rep(1:3, each = 2L),
  lambda_fraction = fraction))
stopifnot(cv$index == 1L, all(is.finite(cv$fold_loss)), !cv$partial_search,
  !cv$ok, cv$status == "join_unavailable", is.na(cv$nzero),
  all(vapply(cv$selected_curvature, function(h) h$status == "join_unavailable", logical(1))),
  fails(coef(cv)), fails(predict(cv, newx = x)))

# Away from the join, the same one-variable design remains usable.
plain <- cv.lambert(x, y, foldid = rep(1:3, each = 2L), lambda_fraction = .2)
stopifnot(plain$ok, plain$status == "verified_stationary",
  all(is.finite(predict(plain, newx = x))))
cat("PASS: strict curvature acceptance, joint-descent counterexample, path/CV joins, prediction refusal and ordinary fits.\n")
