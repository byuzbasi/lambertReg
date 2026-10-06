library(lambertReg)
rejects <- function(expr) inherits(tryCatch(force(expr), error = identity), "error")
near <- function(x, y, tol = 1e-10) isTRUE(all.equal(x, y, tolerance = tol, check.attributes = FALSE))
# Scalar rule, symmetry, continuity and global scalar minimization.
z <- c(-4, -exp(1), -2, -1, 0, 1, 2, exp(1), 4)
t <- lambert_threshold(z, 1)
stopifnot(near(t, -rev(t)), t[4] == 0, t[6] == 0,
  near(t[7], 2 * log(2)), near(t[8:9], z[8:9]),
  near(lambert_penalty(c(exp(1), 4, 10), 1), rep(expm1(2)/4, 3)),
  lambert_penalty(0, 1) == 0)
for (s in c(0, 0.5, 1, 1.01, 2, exp(1), 4)) {
  objective <- function(b) (s - b)^2/2 + lambert_penalty(abs(b), 1)
  numerical <- optimize(objective, c(-5, 5), tol = 1e-10)
  stopifnot(objective(lambert_threshold(s, 1)) <= numerical$objective + 1e-9)
}
stopifnot(rejects(lambert_penalty(-1, 1)), rejects(lambert_threshold(1, 0)),
  rejects(lambert_threshold(Inf, 1)), rejects(lambert_penalty(1, 1e300)))
# Dense, nonstandard-scale inputs and original-scale prediction.
set.seed(104)
x <- sweep(matrix(rnorm(60 * 6), 60, 6), 2, seq_len(6), "*")
x <- sweep(x, 2, 11:16, "+"); colnames(x) <- paste0("x", 1:6)
y <- 4 + 2 * x[, 1] - x[, 2] + rnorm(60)
foldid <- sample(rep(1:5, each = 12))
grid <- 10^seq(0, -3, length.out = 5)
f <- lambert(x, y, lambda_fraction = grid)
stopifnot(all(f$verified), all(f$beta[, 1] == 0),
  near(predict(f, x), cbind(1, x) %*% coef(f)),
  near(f$scale, sqrt(colMeans(sweep(x, 2, colMeans(x), "-")^2))))
# Both logical starts remain; reuse requires identical, verified initial states.
audit <- f$attempts
stopifnot(nrow(audit) == 2L * length(grid), any(audit$reused),
  all(audit$executed != audit$reused),
  all(audit$selected_start[audit$reused] == "zero"),
  all(audit$ok[audit$reused]), all(audit$elapsed[audit$reused] == 0))
for (k in which(colSums(abs(f$beta_standardized)) > 0) + 1L) {
  if (k <= length(grid)) stopifnot(all(audit$executed[audit$index == k]))
}
rng <- .Random.seed
cv <- cv.lambert(x, y, foldid = foldid, lambda_fraction = grid)
stopifnot(identical(rng, .Random.seed), cv$ok, !cv$partial_search,
  near(cv$cvm, rowMeans(cv$fold_loss)),
  near(predict(cv, x), cbind(1, x) %*% coef(cv)))
for (i in 1:5) stopifnot(near(cv$preprocessing[[i]]$mx, colMeans(x[foldid != i, ])))
# A held-out perturbation cannot affect that fold's centering or scale.
x2 <- x; x2[foldid == 1, 1] <- x2[foldid == 1, 1] + 1
cv2 <- cv.lambert(x2, y, foldid = foldid, lambda_fraction = grid)
stopifnot(identical(cv$preprocessing[[1]], cv2$preprocessing[[1]]))
# Supplied absolute lambda path and fractional path are equivalent.
a <- lambert(x, y, lambda = f$lambda)
stopifnot(near(a$beta, f$beta), near(a$intercept, f$intercept))
# Reproducible generated folds; unequal folds remain arithmetic fold averages.
set.seed(18); c1 <- cv.lambert(x, y, lambda_fraction = c(1, .3))
set.seed(18); c2 <- cv.lambert(x, y, lambda_fraction = c(1, .3))
stopifnot(identical(c1$foldid, c2$foldid), identical(c1$cvm, c2$cvm))
u <- cv.lambert(x, y, foldid = rep(1:5, c(8, 10, 12, 14, 16)), lambda_fraction = c(1, .3))
stopifnot(near(u$cvm, rowMeans(u$fold_loss)))
# Invalid data must fail before silent recoding or fitting.
stopifnot(rejects(lambert(cbind(x, 1), y)), rejects(lambert(x, rep(0, 60))),
  rejects(lambert(x, replace(y, 1, NA))), rejects(lambert(x, y, lambda = c(1, 2))),
  rejects(lambert(x, y, lambda = 1, lambda_fraction = .5)),
  rejects(cv.lambert(x, y, foldid = rep(1, 60))),
  rejects(cv.lambert(x, y, foldid = replace(foldid, 1, 8))),
  rejects(cv.lambert(x, y, foldid = foldid, nfolds = 4)),
  rejects(predict(cv, x[, 6:1])), rejects(coef(f, index = 0)),
  rejects(coef(cv, unsupported = TRUE)))
# Fold-constant column: no automatic removal or use of validation scale.
bad <- x; bad[, 6] <- 0; bad[foldid == 1, 6] <- 1:12
stopifnot(rejects(cv.lambert(bad, y, foldid = foldid, lambda_fraction = c(1, .3))))
# Failed fits remain diagnostic artifacts, and cannot be predicted silently.
budget <- suppressWarnings(lambert(x, y, lambda_fraction = .05, max_sweeps = 1))
stopifnot(all(budget$attempts$executed), !any(budget$attempts$reused),
  !budget$verified[1], rejects(predict(budget, x)), rejects(coef(budget)))
budget_cv <- suppressWarnings(cv.lambert(x, y, foldid = foldid,
  lambda_fraction = .05, max_sweeps = 1))
stopifnot(!budget_cv$ok, budget_cv$partial_search,
  rejects(predict(budget_cv, x)), rejects(coef(budget_cv)))
# p > n smoke: no global-optimum assertion, but accepted fits pass diagnostics.
set.seed(19); hx <- matrix(rnorm(30 * 35), 30, 35); hy <- hx[, 1] + rnorm(30)
h <- lambert(hx, hy, lambda_fraction = c(1, .5))
stopifnot(all(h$verified), all(h$diagnostics$kkt <= 1e-7))
pdf_file <- tempfile(fileext = ".pdf")
grDevices::pdf(pdf_file)
plot(f); plot(cv)
grDevices::dev.off()
stopifnot(file.info(pdf_file)$size > 0)
cat("PASS: scalar, prediction, CV, preprocessing, RNG, rejection, failure and high-dimensional smoke checks.\n")
